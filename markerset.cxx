#include <cstring>
#include <iostream> // for std::cerr
#include <time.h> // for time()
#include <assert.h>
#include <unordered_map>
#include <set>
#include <numeric>
#include <random>
#ifdef USE_NPROF
#include <nvToolsExt.h>
#endif

#include "constants.hpp"
#include "parameters.hpp"
#include "barycentric-fn.hpp"
#include "binaryio.hpp"
#include "markerset.hpp"
#include "mesh.hpp"
#include "geometry.hpp"
#include "utils.hpp"
#include "matprops.hpp"


namespace {

    const int DEBUG = 0;
    const double over_alloc_ratio = 2.0;  // how many extra space to allocate for future expansion

}


MarkerSet::MarkerSet(const std::string& name) :
    _name(name)
{
    _last_id = _nmarkers = 0;
    allocate_markerdata(4 * 1024); // pre-allocate a small amount of markers
}


MarkerSet::MarkerSet(const Param& param, Variables& var, const std::string& name) :
    _name(name)
{
    _last_id = _nmarkers = 0;

    switch ( param.markers.init_marker_option ) {
    case 1:
        random_markers(param, var);
        break;
    case 2:
        regularly_spaced_markers(param, var);
        break;
    default:
        std::cerr << "Error: unknown init_marker_option: " << param.markers.init_marker_option << ". The only valid option is '1'.\n";
        std::exit(1);
        break;
    }

    for( int e = 0; e < var.nelem; e++ ) {
        int num_markers_in_elem = 0;
        for( int i = 0; i < param.mat.nmat; i++ )
            num_markers_in_elem += (*(var.elemmarkers))[e][i];

        if (num_markers_in_elem <= 0) {
            std::cerr << "Error: no marker in element #" << e
                      << ". Please increase the number of markers.\n";
            std::exit(1);
        }
    }
}


MarkerSet::MarkerSet(const Param& param, Variables& var, BinaryInput& bin, const std::string& name) :
    _name(name)
{
    // init from checkpoint file
    read_chkpt_file(var, bin);
}


void MarkerSet::allocate_markerdata( const int max_markers )
{
    _reserved_space = max_markers;
    _eta = new shapefn( max_markers );
    _elem = new int_vec( max_markers );
    _mattype = new int_vec( max_markers );
    _id = new int_vec( max_markers );
    _time = new double_vec( max_markers );
    _z = new double_vec( max_markers );
    _distance = new double_vec( max_markers );
    _slope = new double_vec( max_markers );
    _tmp = new double_vec( max_markers );
}


void MarkerSet::random_eta( double *eta )
{
    // eta for randomly scattered markers within an element.
    // An alternative would be to fix barycentric coordinates and add random perturbations.
    //
    while (1) {
        // sum(eta) == 1 && all components of eta are greater than zero
        double sum = 0;
        for( int n = 0; n < NDIMS; n++ ) {
            eta[n] = (rand()/(double)RAND_MAX);
            sum += eta[n];
        }
        if (sum < 1) {
            eta[NODES_PER_ELEM - 1] = 1 - sum;
            break;
        }
    }
}

void MarkerSet::random_eta_seed( double *eta, int seed )
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double sum = 0;
    for(int n = 0; n < NODES_PER_ELEM; n++) {
        eta[n] = dist(gen);
        sum += eta[n];
    }
    for(int n = 0; n < NODES_PER_ELEM; n++)
        eta[n] /= sum;
}


void MarkerSet::append_markers(AMD_vec &md)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    int nmarker = md.size();
    // Ensure sufficient array size
    while ( _nmarkers + nmarker > _reserved_space ) {
        // Resize the marker-related arrays if necessary.
        const int newsize = _nmarkers * over_alloc_ratio;
        resize( newsize );
    }

    int_vec last_ids(nmarker);
    std::iota(last_ids.begin(), last_ids.end(), _last_id);

    int_vec idxs(nmarker);
    std::iota(idxs.begin(), idxs.end(), _nmarkers);

    #pragma omp parallel for default(none) shared(md, idxs, last_ids, nmarker)
    for (int i = 0; i < nmarker; ++i)
        append_marker_at_i(md[i], idxs[i], last_ids[i]);

    _nmarkers += nmarker;
    _last_id += nmarker;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::append_marker_at_i(AppendMarkerData &md, int idx, int last_id)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    int m = idx;
    std::memcpy((*_eta)[m], md.eta.data(), NODES_PER_ELEM*sizeof(double));
    (*_elem)[m] = md.elem;
    (*_mattype)[m] = md.mattype;
    (*_id)[m] = last_id;
    (*_time)[m] = md.time;
    (*_z)[m] = md.depth;
    (*_distance)[m] = md.distance;
    (*_slope)[m] = md.slope;

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void MarkerSet::append_marker( const double *eta, int el, int mt , double ti, double z, double distance, double slope)
{
    // Ensure sufficient array size
    if( _nmarkers == _reserved_space ) {
        // Resize the marker-related arrays if necessary.
        const int newsize = _nmarkers * over_alloc_ratio;
        resize( newsize );
    }

    int m = _nmarkers;
    std::memcpy((*_eta)[m], eta, NODES_PER_ELEM*sizeof(double));
    (*_elem)[m] = el;
    (*_mattype)[m] = mt;
    (*_id)[m] = _last_id;
    (*_time)[m] = ti;
    (*_z)[m] = z;
    (*_distance)[m] = distance;
    (*_slope)[m] = slope;

    if(DEBUG > 1) {
        std::cout << el << " " << m << " "
                  << _nmarkers << " " << (*_mattype)[_nmarkers] << " "
                  << eta[0] << "+" << eta[1] << "+" << eta[2];
#ifdef THREED
        std::cout << "+" << eta[3]
                  << "=" << (eta[0]+eta[1]+eta[2]+eta[3]) << "\n";
#else
        std::cout << "=" << (eta[0]+eta[1]+eta[2]) << "\n";
#endif
    }

    ++_nmarkers;
    ++_last_id;
}


/*
void MarkerSet::create_marker_in_elem(Variables& var)
{

    for (int i=0; i<_nmarkers; ++i)
        (*var.marker_in_elem)[(*_elem)[i]].push_back(i);

    if (DEBUG) {
        std::cout << "Markers in elements";
        for (int i=0; i<var.nelem; i++) {
            std::cout <<"Element: "<< i  << "\n" << "Markers: ";
            for (size_t j=0; j<(*var.marker_in_elem)[i].size();j++)
                std::cout << (*var.marker_in_elem)[i][j] << " ";
            std::cout << "\n";
        }
    }
}

void MarkerSet::update_marker_in_elem(Variables& var)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    delete var.marker_in_elem;
    var.marker_in_elem = new int_vec2D(var.nelem);

#ifdef USE_NPROF
    nvtxRangePushA("loop all marker");
#endif
    for (int i=0; i<_nmarkers; ++i)
        (*var.marker_in_elem)[(*_elem)[i]].push_back(i);
#ifdef USE_NPROF
    nvtxRangePop();
#endif

    if (DEBUG) {
        std::cout << "Markers in elements";
        for (int i=0; i<var.nelem; i++) {
            std::cout <<"Element: "<< i  << "\n" << "Markers: ";
            for (size_t j=0; j<(*var.marker_in_elem)[i].size();j++)
                std::cout << (*var.marker_in_elem)[i][j] << " ";
            std::cout << "\n";
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}
*/

void MarkerSet::create_melt_markers(const int mat, int_vec& melt_markers)
{
    melt_markers.clear();

    for (int i=0;i<_nmarkers;i++)
        if ((*_mattype)[i] == mat)
                melt_markers.push_back(i);

}

// set marker  generated by deposit
void MarkerSet::set_surface_marker(const Param& param,const Variables& var, const double smallest_size, \
            const int mattype, double_vec& edvacc, int_vec2D& elemmarkers, int_vec2D& markers_in_elem)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    double_vec &src_locs = *var.surfinfo.src_locs;

//    #pragma omp parallel for default(none) shared(var,mattype,edvacc,elemmarkers,src_locs,std::cout)
    //XXX has bug of interval empty sediment marker element
    for (size_t i=0; i<(*var.surfinfo.top_facet_elems).size(); i++) {
        int e = (*var.surfinfo.top_facet_elems)[i];

        if (edvacc[e] < 0.) {
            edvacc[e] = 0.;
            continue;
        }

        // if the ratio of increased Vol. and element Vol. up to certain level,
        // set a marker in the area between edvacc.
        int_vec &a = (*var.elemmarkers)[e];
        double nmarkers;
        if (param.mesh.remeshing_option == 13 && param.mesh.meshing_elem_shape > 0) {
            nmarkers = param.markers.markers_per_element;
        } else {
            nmarkers = std::accumulate(a.begin(), a.end(), 0);
        }

        if ((*var.volume)[e] > nmarkers * edvacc[e]) continue;

        double eta[NDIMS], mcoord[NDIMS] = {}, sum = 0.;

        // random horizontal coordinate (2D:x or 3D:x-y) 
        for (int j = 0; j < NDIMS; j++) {
            eta[j] = (rand()/(double)RAND_MAX);
            sum += eta[j];
        }
        for (int j = 0; j < NDIMS; j++)
            eta[j] /= sum;

        int_vec n(NDIMS);
        for (int j=0; j<NDIMS; j++)
            n[j] = (*var.surfinfo.top_nodes)[(*var.surfinfo.elem_and_nodes)[i][j]];

        for (size_t j=0; j<NDIMS-1; j++)
            for (int k=0; k<NDIMS; k++)
                mcoord[j] += (*var.coord)[ n[k] ][j] * eta[k];

#ifdef THREED
        double base = (((*var.coord)[n[1]][0] - (*var.coord)[n[0]][0])*((*var.coord)[n[2]][1] - (*var.coord)[n[0]][1]) \
            - ((*var.coord)[n[2]][0] - (*var.coord)[n[0]][0])*((*var.coord)[n[1]][1] - (*var.coord)[n[0]][1]))/2.0;
#else
        double base = (*var.coord)[n[0]][0] - (*var.coord)[n[1]][0];
#endif
        base = (base > 0.0) ? base : -base;

        // height of marker
        double dv_apply = (*var.volume)[e] / nmarkers;
        double edh = dv_apply / base;
        double dhr = 0.8;
        for (int j=0; j<NDIMS; j++)
            mcoord[NDIMS-1] += ( (*var.coord)[ n[j] ][NDIMS-1] - (edh * dhr) ) * eta[j];

        const double *coord1[NODES_PER_ELEM];

        for (int j=0; j<NODES_PER_ELEM;j++)
            coord1[j] = (*var.coord)[(*var.connectivity)[e][j]];

        Barycentric_transformation bary(coord1, (*var.volume)[e]);

        double eta0[NDIMS], eta1[NODES_PER_ELEM];
        eta1[NODES_PER_ELEM-1] = 1;
        bary.transform(mcoord,0,eta0);

        edvacc[e] -= dv_apply;

        // the surface slope for marker
#ifdef THREED
        double slope, aspect_deg;
        const double *p0 = (*var.coord)[n[0]];
        const double *p1 = (*var.coord)[n[1]];
        const double *p2 = (*var.coord)[n[2]];
        double v1x = p1[0] - p0[0], v1y = p1[1] - p0[1], v1z = p1[2] - p0[2];
        double v2x = p2[0] - p0[0], v2y = p2[1] - p0[1], v2z = p2[2] - p0[2];

        // norm n = v1 x v2
        double nx = v1y * v2z - v1z * v2y;
        double ny = v1z * v2x - v1x * v2z;
        double nz = v1x * v2y - v1y * v2x;

        double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (norm == 0.0) {
            slope = 0.0;
            aspect_deg = -1.0;
        } else {
            // slope: normal vector and z-axis angle
            slope = std::acos(nz / norm);
            // asmith: direction of maximum slope (relative to north, clockwise)
            aspect_deg = std::atan2(ny, nx) * 180.0 / M_PI;
            aspect_deg = 90.0 - aspect_deg;
            if (aspect_deg < 0) aspect_deg += 360.0;
        }
#else
        double slope = ( (*var.coord)[n[0]][NDIMS-1] - (*var.coord)[n[1]][NDIMS-1] ) / ( (*var.coord)[n[0]][0] - (*var.coord)[n[1]][0] );
#endif
        double water_depth = var.surfinfo.base_level - mcoord[NDIMS-1];

        double distance = 0.;
        if (mcoord[0] >= src_locs[0] && mcoord[0] <= src_locs[1]) {
            double dx0 = mcoord[0] - src_locs[0];
            double dx1 = src_locs[1] - mcoord[0];
            if (slope < 0)
                distance = dx0;
            else if (slope > 0)
                distance = dx1;
            else
                distance =  std::min(dx0, dx1);
        }

        if (bary.is_inside(eta0)) {
            for (int j=0; j<NDIMS; j++) {
                eta1[j] = eta0[j];
                eta1[NDIMS] -= eta0[j];
            }

//            #pragma omp critical
            {
                append_marker(eta1, e, mattype, var.time / YEAR2SEC, water_depth, distance, slope);
                ++elemmarkers[e][mattype];
                markers_in_elem[e].push_back(_nmarkers-1);
            }
//            #pragma omp critical(deposit)
//            std::cout << "   A mattype " << mattype << " marker generated by surface processes in element " << e << ".\n";
        }
        else {
            int new_elem, inc;

//            #pragma omp critical(deposit)
            std::cout << "   A generated  marker (mat=" << mattype << ") in element " << e << " is trying to remap in element ";

            remap_marker(var,mcoord,e,new_elem,eta0,inc);

            if (inc) {
                for (int j=0; j<NDIMS; j++) {
                    eta1[j] = eta0[j];
                    eta1[NDIMS] -= eta0[j];
                }

//                #pragma omp critical(add_marker)
                {
                    append_marker(eta1, new_elem, mattype, var.time / YEAR2SEC, water_depth, distance, slope);
                    ++elemmarkers[new_elem][mattype];
                    markers_in_elem[new_elem].push_back(_nmarkers-1);
                }
//                #pragma omp critical
                std::cout << "... Success!\n";
            }
            else {
                std::cout << "... Surface marker generated fail!\n";

//                #pragma omp critical(deposit_err)
                {
                    std::cout << "Coordinate: ";
                    for (int j=0; j<NDIMS; j++)
                        std::cout << mcoord[j] << " ";
                    std::cout << "\neta: ";
                    for (int j=0; j<NDIMS; j++)
                        std::cout << j << " " << eta0[j] << " ";
                    std::cout << "\n";
                }
            }
        }
    }
//    printf("\n");
#ifdef USE_NPROF
    nvtxRangePop();
#endif

}


void MarkerSet::remap_marker(const Variables &var, const double *m_coord, \
                    const int e, int& new_elem, double *new_eta, int& inc)
{
    const double *coord[NODES_PER_ELEM];
    double eta[NODES_PER_ELEM];
    int_vec nodes((*var.connectivity)[e],(*var.connectivity)[e]+NODES_PER_ELEM);
    int_vec searched(var.nelem,0);
    searched[e]=1;

//    std::cout << "Try to remap in ";

    for (auto n = nodes.begin(); n<nodes.end(); n++) {
        for(auto ee = (*var.support)[*n].begin(); ee < (*var.support)[*n].end(); ++ee) {
            if (searched[*ee]) continue;
            searched[*ee]=1;

//            std::cout << *ee << " ";

            for (int j=0; j<NODES_PER_ELEM; j++)
                coord[j] = (*var.coord)[(*var.connectivity)[*ee][j]];

            double volume = compute_volume(coord);
            Barycentric_transformation bary(coord,volume);
            bary.transform(m_coord,0,eta);

            if (bary.is_inside(eta)) {
                for (int j=0; j<NDIMS; j++)
                    new_eta[j] = eta[j];

                new_elem = *ee;
                inc = 1;
                return;
            }
            for (int j=0; j<NODES_PER_ELEM; j++)   
                coord[j] = NULL;
        }
    }
    inc = 0;
}

void MarkerSet::append_random_marker_in_elem( int el, int mt)
{
    double eta[NODES_PER_ELEM];
    random_eta(eta);
    append_marker(eta, el, mt, 0., 0., 0., 0.);
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void MarkerSet::append_random_marker_in_elem( int el, int mt, double time)
{
    double eta[NODES_PER_ELEM];
    random_eta(eta);
    append_marker(eta, el, mt, time, 0., 0., 0.);
}

void MarkerSet::random_markers( const Param& param, Variables &var )
{
    const int ne = var.nelem;
    const int mpe = param.markers.markers_per_element;
    const int num_markers = ne*mpe;
    const int max_markers = num_markers * over_alloc_ratio;

    // allocate memory for data members.
    allocate_markerdata( max_markers );
    
    // initialize random seed:
    if (param.markers.random_seed)
        srand(param.markers.random_seed);
    else
        srand(time(NULL));

    // Generate particles in each element.
    for( int e = 0; e < ne; e++ )
        for( int m = 0; m < mpe; m++ ) {
            // random barycentric coordinate
            double eta[NODES_PER_ELEM];
            random_eta(eta);

            // decide the mattype of markers
            int mt = initial_mattype(param, var, e, eta);
            append_marker(eta, e, mt, 0., 0., 0., 0.);
            ++(*var.elemmarkers)[e][mt];
            (*var.markers_in_elem)[e].push_back(_nmarkers-1);
        }
}


void MarkerSet::regularly_spaced_markers( const Param& param, Variables &var )
{
    const int d = param.markers.init_marker_spacing * param.mesh.resolution;

    double domain_min[NDIMS], domain_max[NDIMS];
    {
        for (int d=0; d<NDIMS; d++)
            domain_min[d] = domain_max[d] = (*var.coord)[0][d];

        for (int i=1; i<var.nnode; i++) {
            for (int d=0; d<NDIMS; d++) {
                domain_min[d] = std::min(domain_min[d], (*var.coord)[i][d]);
                domain_max[d] = std::max(domain_max[d], (*var.coord)[i][d]);
            }
        }
        // print(std::cout, domain_min, NDIMS);
        // print(std::cout, domain_max, NDIMS);
    }

    const double xlength = domain_max[0] - domain_min[0];
    const int nx = xlength / d + 1;
    const double x0 = domain_min[0] + 0.5 * (xlength - (nx-1)*d);
    const double zlength = domain_max[NDIMS-1] - domain_min[NDIMS-1];
    const int nz = zlength / d + 1;
    const double z0 = domain_min[NDIMS-1] + 0.5 * (zlength - (nz-1)*d);
#ifdef THREED
    const double ylength = domain_max[1] - domain_min[1];
    const int ny = ylength / d + 1;
    const double y0 = domain_min[1] + 0.5 * (ylength - (ny-1)*d);
#else
    const int ny = 1;
#endif

    const int num_markers = nx * ny * nz;
    const int max_markers = num_markers * over_alloc_ratio;

    allocate_markerdata( max_markers );

    // nearest-neighbor search structure
    array_t *centroid = elem_center(*var.coord, *var.connectivity); // centroid of elements

    PointCloud cloud(*centroid);
    NANOKDTree kdtree(NDIMS, cloud);
    kdtree.buildIndex();

    const int k = std::min(20, var.nelem);  // how many nearest neighbors to search?

    double_vec new_volume( var.nelem );
    compute_volume( *var.coord, *var.connectivity, new_volume );
    Barycentric_transformation bary(*var.coord, *var.connectivity, new_volume);

    #pragma acc wait

    for (int n=0; n< num_markers; ++n) {
        int ix = n % nx;
        int iy = (n / nx) % ny;
        int iz = n / (nx * ny);

        // Physical coordinate of new marker
        double x[NDIMS] = { x0 + ix*d,
#ifdef THREED
                            y0 + iy*d,
#endif
                            z0 + iz*d };
        bool found = false;

        // Look for nearby elements.
        size_t_vec nn_idx(k);
        double_vec dd(k);
        KNNResultSet resultSet(k);
        resultSet.init(nn_idx.data(), dd.data());

        kdtree.findNeighbors(resultSet, x);

        for( int j=0; j<k; j++ ) {
            int e = nn_idx[j];
            double eta[NODES_PER_ELEM];

            bary.transform(x, e, eta);

            // Compute the last component of eta, with the constraint sum(eta)==1
            double tmp = 1;
            for( int d=0; d<NDIMS; ++d) {
                tmp -= eta[d];
            }
            eta[NDIMS] = tmp;

            if (bary.is_inside(eta)) {
                int mt = initial_mattype(param, var, e, eta, x);
                append_marker(eta, e, mt, 0., 0., 0., 0.);
                ++(*var.elemmarkers)[e][mt];
                (*var.markers_in_elem)[e].push_back(_nmarkers-1);
                found = true;
                break;
            }
        }

        if (! found) {
            // x is outside the domain (ex: the domain is not rectangular)
            continue;
        }
    }

    delete centroid;
}


int MarkerSet::initial_mattype( const Param& param, const Variables &var,
                                int elem, const double eta[NODES_PER_ELEM],
                                const double *x)
{
    int mt;
    if (param.ic.mattype_option == 0) {
        mt = (*var.regattr)[elem][0]; // mattype should take a reginal attribute assigned during meshing.
    }
    else {
        double p[NDIMS] = {0};
        if (x) {
            std::memcpy(p, x, NDIMS*sizeof(double));
        }
        else {
            const int *conn = (*var.connectivity)[elem];
            for(int i=0; i<NDIMS; i++) {
                for(int j=0; j<NODES_PER_ELEM; j++)
                    p[i] += (*var.coord)[ conn[j] ][i] * eta[j];
            }
        }
        // modify mt according to the marker coordinate p
        switch (param.ic.mattype_option) {
        case 1:
            mt = layered_initial_mattype(param, var, elem, eta, p);
            break;
        case 101:
            mt = custom_initial_mattype(param, var, elem, eta, p);
            break;
        default:
            std::cerr << "Error: unknown ic.mattype_option: " << param.ic.mattype_option << '\n';
            std::exit(1);
        }
    }
    return mt;
}


int MarkerSet::layered_initial_mattype( const Param& param, const Variables &var,
                                        int elem, const double eta[NODES_PER_ELEM],
                                        const double *x)
{
    int mt = param.ic.layer_mattypes[param.ic.layer_mattypes.size() - 1];
    const double_vec &layers = param.ic.mattype_layer_depths;
    for (std::size_t i=0; i<layers.size(); ++i) {
        if (x[NDIMS-1] >= -param.mesh.zlength * layers[i]) {
            mt = param.ic.layer_mattypes[i];
            break;
        }
    }
    return mt;
}


int MarkerSet::custom_initial_mattype( const Param& param, const Variables &var,
                                       int elem, const double eta[NODES_PER_ELEM],
                                       const double *x )
{
    /* User defined function */
    int mt = 0;

    return mt;
}

void MarkerSet::remove_markers(int_vec& markers, int_vec2D& markers_in_elem)
{
    if (markers.empty()) return;

#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    std::sort(markers.begin(),markers.end());

    int n = markers.size();
    int_vec replaced_markers(n), a_out, b_out;
    a_out.reserve(n);
    b_out.reserve(n);

    #pragma omp parallel default(none) shared(markers_in_elem,replaced_markers, markers, n, a_out, b_out)
    {
        #pragma omp for
        for (int i=0; i<n; i++)
            replaced_markers[i] = _nmarkers - n + i;

        #pragma omp sections
        {
            #pragma omp section
            std::set_difference(markers.begin(), markers.end(),
                        replaced_markers.begin(), replaced_markers.end(),
                        std::back_inserter(a_out));

            #pragma omp section
            std::set_difference(replaced_markers.begin(), replaced_markers.end(),
                        markers.begin(), markers.end(),
                        std::back_inserter(b_out));
        }

#ifndef ACC
        #pragma omp for
#else
        #pragma omp single
        #pragma acc parallel loop async
#endif
        for (int i = 0; i < a_out.size(); i++) {
            int_vec& emarkers = markers_in_elem[(*_elem)[b_out[i]]];
            auto it = std::find(emarkers.begin(), emarkers.end(), b_out[i]);
            emarkers[it - emarkers.begin()] = a_out[i];
            remove_marker_data(a_out[i],b_out[i]);
        }
    }

    printf("    Removed %d markers from markerset.\n", n);

    #pragma acc wait

    _nmarkers -= n;

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::remove_marker(int i)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // Replace marker i by the last marker.
    --_nmarkers;
    remove_marker_data(i, _nmarkers);
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::remove_marker_data(int is, int ie)
{
    // Replace marker i by the target marker ie.
    std::memcpy( (*_eta)[is], (*_eta)[ie], sizeof(double)*(NODES_PER_ELEM) );
    (*_id)[is] = (*_id)[ie];
    (*_elem)[is] = (*_elem)[ie];
    (*_mattype)[is] = (*_mattype)[ie];
    (*_time)[is] = (*_time)[ie];
    (*_z)[is] = (*_z)[ie];
    (*_distance)[is] = (*_distance)[ie];
    (*_slope)[is] = (*_slope)[ie];
    (*_tmp)[is] = (*_tmp)[ie];
}


void MarkerSet::set_eta( const int i, const double r[NDIMS] ) {
    double sum = 0.0;
    for( int j = 0; j < NDIMS; j++ ) {
        (*_eta)[i][j] = r[j];
        sum += r[j];
    }
    (*_eta)[i][NDIMS] = 1.0 - sum;
}


void MarkerSet::resize( const int newsize )
{
    if( newsize > _reserved_space ) {
        // enlarge arrays
        std::cout << "  Increasing marker arrays size to " << newsize << " markers.\n";
        _reserved_space = newsize;

        double *new_eta = new double[NODES_PER_ELEM * newsize];
        std::copy( (*_eta)[0], (*_eta)[_nmarkers], new_eta );
        _eta->reset( new_eta, newsize );

        _elem->resize( newsize );
        _mattype->resize( newsize );
        _id->resize( newsize );
        _time->resize( newsize );
        _z->resize( newsize );
        _distance->resize( newsize );
        _slope->resize( newsize );
        _tmp->resize( newsize );
    }
    // else if( nmarkers_new < _reserved_space ) {
    //     // TBD: shrink arrays
    //     _reserved_space = newsize;
    // }
    else {
        // New size is too close to old size, don't do anything.
    }
}


void MarkerSet::write_chkpt_file(BinaryOutput &bin) const
{
    int_vec itmp(2);
    itmp[0] = _nmarkers;
    itmp[1] = _last_id;
    bin.write_array(itmp, (_name + " size").c_str(), itmp.size());

    bin.write_array(*_eta, (_name + ".eta").c_str(), _nmarkers);
    bin.write_array(*_elem, (_name + ".elem").c_str(), _nmarkers);
    bin.write_array(*_mattype, (_name + ".mattype").c_str(), _nmarkers);
    bin.write_array(*_id, (_name + ".id").c_str(), _nmarkers);
    bin.write_array(*_time, (_name + ".time").c_str(), _nmarkers);
    bin.write_array(*_z, (_name + ".z").c_str(), _nmarkers);
    bin.write_array(*_distance, (_name + ".distance").c_str(), _nmarkers);
    bin.write_array(*_slope, (_name + ".slope").c_str(), _nmarkers);

}


void MarkerSet::read_chkpt_file(Variables &var, BinaryInput &bin)
{
    int_vec itmp(2);
    bin.read_array(itmp, (_name + " size").c_str());
    _nmarkers = itmp[0];
    _last_id = itmp[1];

    allocate_markerdata(_nmarkers);

    if (_nmarkers != 0) {
        bin.read_array(*_eta, (_name + ".eta").c_str());
        bin.read_array(*_elem, (_name + ".elem").c_str());
        bin.read_array(*_mattype, (_name + ".mattype").c_str());
        bin.read_array(*_id, (_name + ".id").c_str());
        bin.read_array(*_time, (_name + ".time").c_str());
        bin.read_array(*_z, (_name + ".z").c_str());
        bin.read_array(*_distance, (_name + ".distance").c_str());
        bin.read_array(*_slope, (_name + ".slope").c_str());
    }

    if (_name == "markerset")
        for( int i = 0; i < _nmarkers; i++ ) {
            int e = (*_elem)[i];
            int mt = (*_mattype)[i];
            ++(*var.elemmarkers)[e][mt];
            (*var.markers_in_elem)[e].push_back(i);
        }
    else if (_name == "hydrous-markerset")
        for( int i = 0; i < _nmarkers; i++ ) {
            int e = (*_elem)[i];
            ++(*var.hydrous_elemmarkers)[e][0];
            (*var.hydrous_markers_in_elem)[e].push_back(i);
        }
}


void MarkerSet::write_save_file(const Variables &var, BinaryOutput &bin) const
{
#ifdef USE_NPROF
    nvtxRangePushA("write markersets");
#endif
    int_vec itmp(1);
    itmp[0] = _nmarkers;
    bin.write_array(itmp, (_name + " size").c_str(), itmp.size());

    array_t *mcoord = calculate_marker_coord(var); // coordinate of markers

    bin.write_array(*mcoord, (_name + ".coord").c_str(), _nmarkers);
    bin.write_array(*_eta, (_name + ".eta").c_str(), _nmarkers);
    bin.write_array(*_elem, (_name + ".elem").c_str(), _nmarkers);
    bin.write_array(*_mattype, (_name + ".mattype").c_str(), _nmarkers);
    bin.write_array(*_id, (_name + ".id").c_str(), _nmarkers);
    bin.write_array(*_time, (_name + ".time").c_str(), _nmarkers);
    bin.write_array(*_z, (_name + ".z").c_str(), _nmarkers);
    bin.write_array(*_distance, (_name + ".distance").c_str(), _nmarkers);
    bin.write_array(*_slope, (_name + ".slope").c_str(), _nmarkers);

    delete mcoord;

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::get_ZPT(const Param& param, const Variables& var, int m, double &Z, double &P, double &T) const {
        // Get depth and temperature at the marker
        Z = T = 0;
        const double* eta = (*_eta)[m];
        const int *conn = (*var.connectivity)[(*_elem)[m]];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            Z += (*var.coord)[conn[i]][NDIMS-1] * eta[i];
            T += (*var.temperature)[conn[i]] * eta[i];
        }

        // Get pressure, which is constant in the element
        // P = - trace((*var.stress)[e]) / NDIMS;
        P = ref_pressure(param, Z);    
}


array_t* MarkerSet::calculate_marker_coord(const Variables &var) const {
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // const MarkerSet &ms = *var.markersets[0];
    const int nmarkers = get_nmarkers();
    array_t* points = new array_t(nmarkers);

    #pragma omp parallel for default(none) shared(var, points) firstprivate(nmarkers)
    for (int n=0; n<nmarkers; n++) {
        const int e = get_elem(n);
        const double* eta = get_eta(n);
        const int* conn = (*var.connectivity)[e];

        for(int d=0; d<NDIMS; d++) {
            double sum = 0;
            for(int k=0; k<NODES_PER_ELEM; k++)
                sum += (*var.coord)[ conn[k] ][d] * eta[k];
            (*points)[n][d] = sum;
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif
    return points;
}

namespace {

    template <class T>
    void find_markers_in_element(const Param &param, MarkerSet& ms, T& elemmarkers, int_vec2D& markers_in_elem,
                                 KDTree& kdtree, const Barycentric_transformation &bary,
                                 const array_t& old_coord, const conn_t &old_connectivity)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        const int k = std::min((std::size_t) 20, old_connectivity.size());  // how many nearest neighbors to search?

        int last_marker = ms.get_nmarkers();

#ifdef ACC
        double3_vec queries(last_marker);

        #pragma acc parallel loop async
        for (int i = 0; i < last_marker; i++) {
            bool found = false;

            // 1. Get physical coordinates, x, of an old marker.
            int eold = ms.get_elem(i);
            double x[NDIMS] = {0};
            for (int j = 0; j < NDIMS; j++)
                for (int k = 0; k < NODES_PER_ELEM; k++)
                    x[j] += ms.get_eta(i)[k]*
                        old_coord[ old_connectivity[eold][k] ][j];

            queries[i] = {x[0], x[1], x[2]};
        }


        neighbor_vec neighbors(last_marker*k);

        #pragma acc wait

        kdtree.search_grid(queries, neighbors, k, 3);
#endif

        // Loop over all the old markers and identify a containing element in the new mesh.
        int_vec removed_markers;

        #pragma omp parallel default(none) shared(param, ms, elemmarkers, markers_in_elem, kdtree, \
            bary, old_coord, old_connectivity, last_marker,removed_markers)
        {
            int_vec removed_local;
            int_map emarker_local;

            #pragma omp for nowait
            for (int i = 0; i < last_marker; i++) {
                bool found = false;

                // 1. Get physical coordinates, x, of an old marker.
                int eold = ms.get_elem(i);
                double x[NDIMS] = {0};
#ifdef ACC
                x[0] = queries[i].x;
                x[1] = queries[i].y;
#ifdef THREED
                x[2] = queries[i].z;
#endif

                // 2. look for nearby elements.
                neighbor* nn_idx = neighbors.data() + i*k;
#else
                for (int j = 0; j < NDIMS; j++)
                    for (int k = 0; k < NODES_PER_ELEM; k++)
                        x[j] += ms.get_eta(i)[k]*
                            old_coord[ old_connectivity[eold][k] ][j];

                // 2. look for nearby elements.
                size_t_vec nn_idx(k);
                double_vec out_dists_sqr(k);
                KNNResultSet resultSet(k);
                resultSet.init(nn_idx.data(), out_dists_sqr.data());

                kdtree.findNeighbors(resultSet, x);
#endif
                for( int j = 0; j < k; j++ ) {
#ifdef ACC
                    int e = nn_idx[j].idx;
#else
                    int e = nn_idx[j];
#endif
                    double r[NDIMS];

                    bary.transform(x, e, r);

                    // change this if-condition to (i == N) to debug the N-th marker
                    if (bary.is_inside(r)) {
                        ms.set_eta(i, r);
                        ms.set_elem(i, e);
                        emarker_local[e*param.mat.nmat + ms.get_mattype(i)]++;

                        #pragma omp critical
                        markers_in_elem[e].push_back(i);
                
                        found = true;
                        ms.set_tmp(i, 1.0);
                        break;
                    }
                }

                if( found ) continue;

                /* not found */
                // Since no containing element has been found, delete this marker.
                removed_local.emplace_back(i);
                ms.set_tmp(i, 0.0); // mark this marker for removal
            }

            #pragma omp critical
            removed_markers.insert(removed_markers.end(), removed_local.begin(), removed_local.end());

            #pragma omp critical
            for (const auto& pair : emarker_local) {
                int e = pair.first / param.mat.nmat;
                int mt = pair.first % param.mat.nmat;
                elemmarkers[e][mt] += pair.second;
            }
        }

        ms.remove_markers(removed_markers, markers_in_elem);

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_0(const Param& param, const Variables &var,
                                          int_pair_vec &unplenished_elems)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        for (const auto& pair : unplenished_elems) {
            int e = pair.first;
            int num_marker_in_elem = pair.second;

            while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
                const int mt = 0;
                var.markersets[0]->append_random_marker_in_elem(e, mt);
                if (DEBUG) {
                    std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
                }

                ++(*var.elemmarkers)[e][mt];
                (*var.markers_in_elem)[e].push_back(var.markersets[0]->get_nmarkers()-1);
                ++num_marker_in_elem;
            }
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_from_cpdf(const Param& param, const Variables &var,
                                                  int_pair_vec &unplenished_elems)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        for (const auto& pair : unplenished_elems) {
            int e = pair.first;
            int num_marker_in_elem = pair.second;

            // cummulative probability density function of mattype
            double_vec cpdf(param.mat.nmat, 0);

            if (num_marker_in_elem > 0) {
                // cpdf of this element
                cpdf[0] = (*(var.elemmarkers))[e][0] / double(num_marker_in_elem);
                for( int i = 1; i < param.mat.nmat - 1; i++ )
                    cpdf[i] = cpdf[i-1] + (*(var.elemmarkers))[e][i] / double(num_marker_in_elem);
            }
            else {
                // No markers in this element.
                // Construct cpdf from neighboring elements
                int num_markers_in_nbr_elems = 0;
                // Looping over all neighboring elements (excluding self)
                for( int kk = 0; kk < NODES_PER_ELEM; kk++) {
                    int n = (*var.connectivity)[e][kk]; // node of this element
                    for( auto ee = (*var.support)[n].begin(); ee < (*var.support)[n].end(); ++ee) {
                        if (*ee == e) continue;
                        // Note: some (NODES_PER_ELEM) elements will be iterated
                        // more than once (NDIMS times). These elements are direct neighbors,
                        // i.e. they share facets (3D) or edges (2D) with element e.
                        // So they are counted multiple times to reprensent a greater weight.
                        for( int i = 0; i < param.mat.nmat; i++ ) {
                            cpdf[i] += (*(var.elemmarkers))[*ee][i];
                            num_markers_in_nbr_elems += (*(var.elemmarkers))[*ee][i];
                        }
                    }
                }
                for( int i = 1; i < param.mat.nmat - 1; i++ )
                    cpdf[i] += cpdf[i-1];
                for( int i = 0; i < param.mat.nmat - 1; i++ )
                    cpdf[i] = cpdf[i] / double(num_markers_in_nbr_elems);
            }

            cpdf[param.mat.nmat - 1] = 1; // fix to 1 to avoid round-off error
            if (DEBUG > 1) {
                std::cout << num_marker_in_elem << " markers in element " << e << '\n'
                        << "  cpdf: ";
                print(std::cout, cpdf);
                std::cout << '\n';
            }

            while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
                // Determine new marker's matttype based on cpdf
                auto upper = std::upper_bound(cpdf.begin(), cpdf.end(), rand()/(double)RAND_MAX);
                const int mt = upper - cpdf.begin();
                var.markersets[0]->append_random_marker_in_elem(e, mt);
                if (DEBUG) {
                    std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
                }

                ++(*var.elemmarkers)[e][mt];
                (*var.markers_in_elem)[e].push_back(var.markersets[0]->get_nmarkers()-1);
                ++num_marker_in_elem;
            }
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_from_nn(const Param& param, const Variables &var,
                                                int_pair_vec &unplenished_elems)
    {
        if (unplenished_elems.empty()) return;

#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        array_t *points = var.markersets[0]->calculate_marker_coord(var); // coordinate of markers
        PointCloud cloud(*points);

#ifdef USE_NPROF
        nvtxRangePushA("create kdtree for markers");
#endif
        NANOKDTree kdtree(NDIMS, cloud);
        kdtree.buildIndex();
#ifdef USE_NPROF
        nvtxRangePop();
#endif
        MarkerSet &ms = *var.markersets[0];
        AMD_vec marker_data_all;
        int ne = unplenished_elems.size();

        #pragma omp parallel default(none) shared(param, var, unplenished_elems, kdtree, ms, marker_data_all,ne)
        {
            AMD_vec marker_data_local;

            #pragma omp for nowait
            for (int i=0; i<ne; ++i) {
                int e = unplenished_elems[i].first;
                int num_marker_in_elem = unplenished_elems[i].second;

                while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
                    double eta[NODES_PER_ELEM];
                    ms.random_eta_seed(eta, e+num_marker_in_elem+var.steps);

                    double x[NDIMS] = {0};
                    const int *conn = (*var.connectivity)[e];
                    for (int d=0; d<NDIMS; d++) {
                        for (int ii=0; ii<NODES_PER_ELEM; ii++) {
                            x[d] += (*var.coord)[ conn[ii] ][d] * eta[ii];
                        }
                    }

                    const int k = 1;  // how many nearest neighbors to search?
                    size_t_vec nn_idx(k);
                    double_vec dd(k);
                    KNNResultSet resultSet(k);
                    resultSet.init(nn_idx.data(), dd.data());

                    // Look for nearest marker.
                    kdtree.findNeighbors(resultSet, x);

                    int m = nn_idx[0]; // nearest marker
                    const int mt = ms.get_mattype(m);
    //            const double ti = ms.get_time(m);

                    // ms.append_marker(eta, e, mt, 0., 0., 0., 0.);
                    // (*var.markers_in_elem)[e].push_back(ms.get_nmarkers()-1);

                    double_vec eta_tmp(eta, eta + NODES_PER_ELEM);
                    marker_data_local.push_back({eta_tmp, e, mt, 0., 0., 0., 0.});

                    ++(*var.elemmarkers)[e][mt];
                    ++num_marker_in_elem;
                }
            }

            #pragma omp critical
            marker_data_all.insert(marker_data_all.end(), marker_data_local.begin(), marker_data_local.end());
        }

        std::stable_sort(marker_data_all.begin(), marker_data_all.end(),
            [](const AppendMarkerData& a, const AppendMarkerData& b) {
                return a.elem < b.elem;
            });

        // Append new markers to the end of the marker set.
        ms.append_markers(marker_data_all);

        int nnew = marker_data_all.size();
        int nmarkers = ms.get_nmarkers();
        for (int i=0; i<nnew; ++i)
            (*var.markers_in_elem)[marker_data_all[i].elem].push_back(nmarkers-nnew+i);

        delete points;

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

} // anonymous namespace


void MarkerSet::check_marker_elem_consistency(const Variables &var) const
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    #pragma acc serial
    int ncount = 0, is_error = 0;
#ifndef ACC
    #pragma omp parallel for reduction(+:ncount,is_error) default(none) shared(var,std::cerr)
#endif
    #pragma acc parallel loop reduction(+:ncount,is_error)
    for (int e=0; e<var.nelem; ++e) {
        int nmarker_mat = std::accumulate((*var.elemmarkers)[e].begin(), (*var.elemmarkers)[e].end(), 0);
        int elenmarkers = (*var.markers_in_elem)[e].size();

        if (elenmarkers != nmarker_mat) {
            printf("Error: markers count mismatch in element %d: %d vs %d\n", e, elenmarkers, nmarker_mat);
            ++is_error;
        }
        for (int i=0; i<elenmarkers; ++i) {
            int m = (*var.markers_in_elem)[e][i];
            if (m < 0 || m >= _nmarkers) {
                printf("Error: melt production marker %d in element %d is out of range [0,%d)\n", m, e, _nmarkers);
                ++is_error;
            }

            int me = (*_elem)[m];
            if (me != e) {
                printf("Error: marker %d in element %d is not in element %d\n", m, me, e);
                ++is_error;
            }
        }
        ncount += elenmarkers;
    }

    if (_nmarkers != ncount) {
        std::cerr << "Error: markers count mismatch: " << _nmarkers << " vs " << ncount << '\n';
        std::exit(1);
    } else if (is_error > 0) {
        std::exit(1);
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


// surface processes correcttion of marker
void MarkerSet::correct_surface_marker(const Param &param, const Variables& var, const double_vec& dhacc, int_vec2D &elemmarkers, int_vec2D &markers_in_elem)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // correct surface marker.
    Barycentric_transformation bary(*var.top_elems, *var.coord, *var.connectivity, *var.volume);

    array_t coord0s(var.ntop_elems*NODES_PER_ELEM, 0.);

    int_vec delete_marker;
    int nchange = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var,coord0s,dhacc)
#endif
    #pragma acc parallel loop async
    for (int i=0;i<var.ntop_elems;i++) {
        int* tnodes = (*var.connectivity)[(*var.top_elems)[i]];

        double *c00 = coord0s[i*NODES_PER_ELEM];
        double *c01 = coord0s[i*NODES_PER_ELEM+1];
        double *c02 = coord0s[i*NODES_PER_ELEM+2];
#ifdef THREED
        double *c03 = coord0s[i*NODES_PER_ELEM+3];
#endif

        // restore the reference node locations before deposition/erosion 
        c00[0] = (*var.coord)[tnodes[0]][0];
        c01[0] = (*var.coord)[tnodes[1]][0];
        c02[0] = (*var.coord)[tnodes[2]][0];

        c00[NDIMS-1] = (*var.coord)[tnodes[0]][NDIMS-1] - dhacc[tnodes[0]];
        c01[NDIMS-1] = (*var.coord)[tnodes[1]][NDIMS-1] - dhacc[tnodes[1]];
        c02[NDIMS-1] = (*var.coord)[tnodes[2]][NDIMS-1] - dhacc[tnodes[2]];
#ifdef THREED
        c00[1] = (*var.coord)[tnodes[0]][1];
        c01[1] = (*var.coord)[tnodes[1]][1];
        c02[1] = (*var.coord)[tnodes[2]][1];

        c03[0] = (*var.coord)[tnodes[3]][0];
        c03[1] = (*var.coord)[tnodes[3]][1];
        c03[2] = (*var.coord)[tnodes[3]][2] - dhacc[tnodes[3]];
#endif
    }

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,coord0s,bary,markers_in_elem) reduction(+:nchange)
#endif
    #pragma acc parallel loop reduction(+:nchange) async
    for (int i=0; i<var.ntop_elems;i++) {
        int e = (*var.top_elems)[i];
        int_vec &markers = markers_in_elem[e];
        int nmarker = markers.size();
        for (int j=0; j<nmarker;j++) {
            int m = markers[j];
            (*_tmp)[m] = 0.;
            double m_coord[NDIMS], new_eta[NDIMS];
            for (int k=0; k<NDIMS; k++) {
                m_coord[k] = 0.;
                    for (int l=0; l<NODES_PER_ELEM; l++)
                        m_coord[k] += (*_eta)[m][l] * coord0s[i*NODES_PER_ELEM+l][k];
            }
            // check if the marker is still in original element
            bary.transform(m_coord,i,new_eta);
            if (bary.is_inside(new_eta)) {
                set_eta(m, new_eta);
            } else {
                ++(*_tmp)[m];
                ++nchange;
            }
        }
    }

    #pragma acc wait

    if (nchange > 0) {
        delete_marker.reserve(nchange);

        std::vector<MarkerUpdate> updates;

        #pragma omp parallel default(none) \
                shared(var,coord0s,markers_in_elem, updates, elemmarkers)
        {
            std::vector<MarkerUpdate> local_updates;

            #pragma omp for nowait
            for (int i=0; i<var.ntop_elems;i++) {
                int e = (*var.top_elems)[i];
                int_vec &markers = markers_in_elem[e];
                for (int j=0; j<(int)markers.size();++j) {
                    int m = markers[j];
                    if ((*_tmp)[m] < 1.) continue;

                    double m_coord[NDIMS], new_eta[NDIMS];

                    for (int k=0; k<NDIMS; k++) {
                        m_coord[k] = 0.;
                        for (int l=0; l<NODES_PER_ELEM; l++)
                            m_coord[k] += (*_eta)[m][l] * coord0s[i*NODES_PER_ELEM+l][k];
                    }

                    int inc, new_elem;
                    int mat = (*_mattype)[m];
                    // find new element of the marker
                    remap_marker(var,m_coord,e,new_elem,new_eta,inc);

                    MarkerUpdate u;
                    u.m = m;
                    u.src_elem = e;
                    u.dst_elem = new_elem;
                    u.inc = inc;
                    local_updates.push_back(u);

                    if (inc) {
                        set_eta(m, new_eta);
                        set_elem(m, new_elem);
                        #pragma omp atomic update
                        ++elemmarkers[new_elem][mat];
                    }
                    --elemmarkers[e][mat];
                }
            }

            #pragma omp critical
            updates.insert(updates.end(), local_updates.begin(), local_updates.end());
        }

        // vector operations cannot be parallelized
        for (const auto& u : updates) {
            int_vec& old_markers = markers_in_elem[u.src_elem];
            old_markers.erase(std::find(old_markers.begin(), old_markers.end(), u.m));

            if (u.inc) {
                int_vec& new_markers = markers_in_elem[u.dst_elem];
                auto it = std::lower_bound(new_markers.begin(), new_markers.end(), u.m);
                new_markers.insert(it, u.m);
            } else {
                delete_marker.push_back(u.m);
            }
        }

        // delete recorded marker
        if (!delete_marker.empty())
            remove_markers(delete_marker, markers_in_elem);

        #pragma acc wait

        int_pair_vec unplenished_elems;
        for (int i=0; i<var.ntop_elems; i++) {
            int e = (*var.top_elems)[i];
            int  nmarkers = markers_in_elem[e].size();

            if (nmarkers < param.markers.min_num_markers_in_element)
                unplenished_elems.emplace_back(e, nmarkers);
        }

        if (unplenished_elems.size() > 0) {
            printf("replenish markers in %d elements.\n", (int)unplenished_elems.size());
            for (int i=0; i<(int)unplenished_elems.size(); i++) {
                int e = unplenished_elems[i].first;
                int nmarkers = unplenished_elems[i].second;
                printf("  Element %d has %d markers.\n", e, nmarkers);
            }
        }

        switch (param.markers.replenishment_option) {
        case 0:
            replenish_markers_with_mattype_0(param, var, unplenished_elems);
            break;
        case 1:
            replenish_markers_with_mattype_from_cpdf(param, var, unplenished_elems);
            break;
        case 2:
            replenish_markers_with_mattype_from_nn(param, var, unplenished_elems);
            break;
        default:
            std::cerr << "Error: unknown markers.replenishment_option: " << param.markers.replenishment_option << '\n';
            std::exit(1);
        }

        for (int i=0; i<var.ntop_elems; i++) {
            int e = (*var.top_elems)[i];
            int nemarker0 = std::accumulate((*var.elemmarkers)[e].begin(), (*var.elemmarkers)[e].end(), 0);
            int  nmarkers = (*var.markers_in_elem)[e].size();
            if (nemarker0 != nmarkers) {
                std::cerr << "Error: number of markers in element " << e
                        << " does not match number of elemmarkers: "
                        << nmarkers << " vs. " << nemarker0 << '\n';
                std::exit(1);
            }
        }
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void remap_markers(const Param& param, Variables &var, const array_t &old_coord,
                   const conn_t &old_connectivity)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // Re-create elemmarkers
    delete var.elemmarkers;
    delete var.markers_in_elem;
    if (param.control.has_hydration_processes) {
        delete var.hydrous_elemmarkers;
        delete var.hydrous_markers_in_elem;
    }
    create_elemmarkers( param, var );

    // Locating markers in new elements
    {
        double_vec new_volume( var.nelem );
        compute_volume( *var.coord, *var.connectivity, new_volume );

        Barycentric_transformation bary( *var.coord, *var.connectivity, new_volume );

        // nearest-neighbor search structure
        array_t *centroid = elem_center(*var.coord, *var.connectivity); // centroid of elements

#ifdef USE_NPROF
        nvtxRangePushA("create kdtree for new elements");
#endif

#ifdef ACC
        double3_vec points(var.nelem);
        elem_center3(*var.coord, *var.connectivity,points); // centroid of elements
        CudaKNN kdtree(param, points);
#else
        PointCloud cloud(*centroid);
        KDTree kdtree(NDIMS, cloud);
        kdtree.buildIndex();
#endif

#ifdef USE_NPROF
        nvtxRangePop();
#endif

        find_markers_in_element(param, *var.markersets[0], *var.elemmarkers, *var.markers_in_elem,
                                kdtree, bary, old_coord, old_connectivity);
        if (param.control.has_hydration_processes)
            find_markers_in_element(param, *var.markersets[var.hydrous_marker_index], *var.hydrous_elemmarkers, *var.hydrous_markers_in_elem,
                                    kdtree, bary, old_coord, old_connectivity);

        #pragma acc wait

        delete centroid;
    }

    // If any new element has too few markers, generate markers in them.
#ifdef USE_NPROF
    nvtxRangePushA("find unplenished elements");
#endif
    // unplenish markers
    int_pair_vec unplenished_elems;

    #pragma omp parallel default(none) shared(param, var, unplenished_elems)
    {
        // local pairs
        int_pair_vec local_pairs;

        #pragma omp for nowait
        for( int e = 0; e < var.nelem; e++ ) {
            int num_marker_in_elem = std::accumulate((*var.elemmarkers)[e].begin(), (*var.elemmarkers)[e].end(), 0);

            if (num_marker_in_elem < param.markers.min_num_markers_in_element)
                local_pairs.emplace_back(e, num_marker_in_elem);
        }

        #pragma omp critical
        unplenished_elems.insert(unplenished_elems.end(), 
                                local_pairs.begin(),
                                local_pairs.end());
        #pragma omp barrier

        #pragma omp single
        {
        std::sort(unplenished_elems.begin(), unplenished_elems.end(),
            [](const int_pair &a, const int_pair &b) {
                return a.first < b.first;
            });
        }
    }
#ifdef USE_NPROF
    nvtxRangePop();
#endif

    switch (param.markers.replenishment_option) {
    case 0:
        replenish_markers_with_mattype_0(param, var, unplenished_elems);
        break;
    case 1:
        replenish_markers_with_mattype_from_cpdf(param, var, unplenished_elems);
        break;
    case 2:
        replenish_markers_with_mattype_from_nn(param, var, unplenished_elems);
        break;
    default:
        std::cerr << "Error: unknown markers.replenishment_option: " << param.markers.replenishment_option << '\n';
        std::exit(1);
    }

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


namespace {

    Barycentric_transformation* get_bary_from_cache(std::unordered_map<int, Barycentric_transformation*> &cache,
                                                    int el, const array_t &coordinate, const int *conn,
                                                    double_vec &volume)
    {
        Barycentric_transformation *bary;
        auto search = cache.find(el);
        if (search == cache.end()) {
            const double *coord[NODES_PER_ELEM];
            for(int j=0; j<NODES_PER_ELEM; j++) {
                coord[j] = coordinate[ conn[j] ];
            }
            bary = new Barycentric_transformation(coord, volume[el]);
            cache[el] =  bary;
        }
        else {
            bary = search->second;
        }
        return bary;
    }
}


void advect_hydrous_markers(const Param& param, const Variables& var, double dt_subtotal,
                            MarkerSet& hydms, Array2D<int,1>& hydem)
{
    std::unordered_map<int, Barycentric_transformation*> cache;
    Barycentric_transformation *bary;

    #pragma acc wait // here is not ACC parallelized yet

    int last_marker = hydms.get_nmarkers();
    int m = 0;
    while (m < last_marker) {
        // Find physical coordinate of the marker
        int el = hydms.get_elem(m);
        double x[NDIMS] = {0};
        const int *conn = (*var.connectivity)[el];
        for(int i=0; i<NDIMS; i++) {
            for(int j=0; j<NODES_PER_ELEM; j++)
                x[i] += hydms.get_eta(m)[j] * (*var.coord)[ conn[j] ][i];
        }

        // Advect the marker
        x[NDIMS-1] += dt_subtotal * param.control.hydration_migration_speed;

        // Transform back to barycentric coordinate
        double r[NDIMS];
        bary = get_bary_from_cache(cache, el, *var.coord, conn, *var.volume);
        bary->transform(x, 0, r); // always (local) 0-th element for bary

        if (bary->is_inside(r)) {
            hydms.set_eta(m, r);
            ++m;
            goto next;
        }
        else {
            // Marker has moved out of el. Find the new containing element.
            for(int j=0; j<NODES_PER_ELEM; j++) {
                const int_vec& supp = (*var.support)[ conn[j] ];
                for (std::size_t k=0; k<supp.size(); k++) {
                    int ee = supp[k];
                    const int *conn2 = (*var.connectivity)[ee];
                    bary = get_bary_from_cache(cache, ee, *var.coord, conn2, *var.volume);
                    bary->transform(x, 0, r);
                    if (bary->is_inside(r)) {
                        hydms.set_elem(m, ee);
                        hydms.set_eta(m, r);
                        --hydem[el][0];
                        ++hydem[ee][0];
                        ++m;
                        goto next;
                    }
                }
            }
            // Since no containing element has been found, delete this marker.
            // Note m is not inc'd.
            if (DEBUG) {
                std::cout << m << "-th hydrous marker not in any element\n";
            }
            --last_marker;
            hydms.remove_marker(m);
            --hydem[el][0];
        }
    next:;
    }

    // clean up all Barycentric_transformation instances
    for (auto i=cache.begin(); i!=cache.end(); ++i) {
        delete i->second;
    }
}


