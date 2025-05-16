#include <cstring>
#include <iostream> // for std::cerr
#include <time.h> // for time()
#include <assert.h>
#include <unordered_map>
#include <set>
#include <numeric>
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
void MarkerSet::set_surface_marker(const Param& param,const Variables& var, const double smallest_size, const int mattype, double_vec& edvacc, int_vec2D& elemmarkers)
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

// surface processes correcttion of marker
void MarkerSet::correct_surface_marker(const Variables& var, array_t& coord0s, Barycentric_transformation &bary, int_vec2D &elemmarkers)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif

    int_vec delete_marker;
    delete_marker.reserve(100);

    #pragma omp parallel default(none) shared(var,coord0s,bary,delete_marker,elemmarkers)
    {
        int_vec delete_local;
        int_pair_vec switch_local;

        #pragma omp for nowait
        for (int i=0; i<_nmarkers;i++) {
            double m_coord[NDIMS], new_eta[NDIMS];
            int e = (*_elem)[i];
            auto it = find(var.top_elems->begin(), var.top_elems->end(), e);
            // If element was not found
            if (it == var.top_elems->end())
                continue;
            int e_ind = it - var.top_elems->begin();

            for (int k=0; k<NDIMS; k++) {
                m_coord[k] = 0.;
                for (int j=0; j<NODES_PER_ELEM; j++)
                    m_coord[k] += (*_eta)[i][j] * coord0s[e_ind*NODES_PER_ELEM+j][k];
            }
            // check if the marker is still in original element
            bary.transform(m_coord,e_ind,new_eta);
            if (bary.is_inside(new_eta))
                set_eta(i, new_eta);
            else {
                int inc, new_elem;
                // std::cout << "Marker " << i << " in element " << e_ind << " is trying to remap in element ";
                // find new element of the marker
                remap_marker(var,m_coord,e,new_elem,new_eta,inc);
                if (inc) {
                    set_eta(i, new_eta);
                    set_elem(i, new_elem);
                    switch_local.emplace_back(e, new_elem);
                    // std::cout << "... Success!.\n";
                }
                else {
                    delete_local.emplace_back(i);
                    // std::cout << "... Fail!. (Erosion might have occurred)\n";
                }
            }
        }

        #pragma omp critical
        delete_marker.insert(delete_marker.end(), delete_local.begin(), delete_local.end());

        for (auto it = delete_local.begin(); it != delete_local.end(); ++it)
            #pragma omp atomic
            elemmarkers[(*_elem)[*it]][(*_mattype)[*it]]--;

        for (auto it = switch_local.begin(); it != switch_local.end(); ++it) {
            #pragma omp atomic
            elemmarkers[it->first][(*_mattype)[it->first]]--;
            #pragma omp atomic
            elemmarkers[it->second][(*_mattype)[it->second]]++;
        }

    }

    // delete recorded marker
    remove_markers(delete_marker);

#ifdef USE_NPROF
    nvtxRangePop();
#endif
}


void MarkerSet::remap_marker(const Variables &var, const double *m_coord, \
                    const int e, int& new_elem, double *new_eta, int& inc)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    const double *coord[NODES_PER_ELEM];
    double volume, eta[NODES_PER_ELEM];
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

            compute_volume(coord,volume);
            Barycentric_transformation bary(coord,volume);
            bary.transform(m_coord,0,eta);

            if (bary.is_inside(eta)) {
                for (int j=0; j<NDIMS; j++)
                    new_eta[j] = eta[j];

                new_elem = *ee;
                inc = 1;
#ifdef USE_NPROF
                nvtxRangePop();
#endif
                return;
            }
            for (int j=0; j<NODES_PER_ELEM; j++)   
                coord[j] = NULL;
        }
    }
    inc = 0;
#ifdef USE_NPROF
    nvtxRangePop();
#endif
}

void MarkerSet::append_random_marker_in_elem( int el, int mt)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    double eta[NODES_PER_ELEM];
    random_eta(eta);
    append_marker(eta, el, mt, 0., 0., 0., 0.);
#ifdef USE_NPROF
    nvtxRangePop();
#endif
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
    KDTree kdtree(NDIMS, cloud);
    kdtree.buildIndex();

    const int k = std::min(20, var.nelem);  // how many nearest neighbors to search?

    double_vec new_volume( var.nelem );
    compute_volume( *var.coord, *var.connectivity, new_volume );
    Barycentric_transformation bary(*var.coord, *var.connectivity, new_volume);

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

void MarkerSet::remove_markers(int_vec& markers)
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

    #pragma omp parallel default(none) shared(replaced_markers, markers, n, a_out, b_out)
    {
        #pragma omp for
        for (int i=0; i<n; i++)
            replaced_markers[i] = _nmarkers - n + i;

        #pragma omp sections
        {
            #pragma omp section
            {
                std::set_difference(markers.begin(), markers.end(),
                            replaced_markers.begin(), replaced_markers.end(),
                            std::back_inserter(a_out));
            }
            #pragma omp section
            {
                std::set_difference(replaced_markers.begin(), replaced_markers.end(),
                            markers.begin(), markers.end(),
                            std::back_inserter(b_out));
            }
        }

        #pragma omp for
        for (int i = 0; i < a_out.size(); i++)
            remove_marker_data(a_out[i],b_out[i]);
    }

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
        }
    else if (_name == "hydrous-markerset")
        for( int i = 0; i < _nmarkers; i++ ) {
            int e = (*_elem)[i];
            ++(*var.hydrous_elemmarkers)[e][0];
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
    void find_markers_in_element(MarkerSet& ms, T& elemmarkers,
                                 KDTree& kdtree, const Barycentric_transformation &bary,
                                 const array_t& old_coord, const conn_t &old_connectivity)
    {
#ifdef USE_NPROF
        nvtxRangePushA(__FUNCTION__);
#endif
        const int k = std::min((std::size_t) 20, old_connectivity.size());  // how many nearest neighbors to search?

        // Loop over all the old markers and identify a containing element in the new mesh.
        int last_marker = ms.get_nmarkers();
        int_vec removed_markers;

        #pragma omp parallel default(none) shared(ms, elemmarkers, kdtree, bary, old_coord, old_connectivity, last_marker,removed_markers)
        {
            int_vec removed_local;

            #pragma omp for nowait
            for (int i = 0; i < last_marker; i++) {
                bool found = false;

                // 1. Get physical coordinates, x, of an old marker.
                int eold = ms.get_elem(i);
                double x[NDIMS] = {0};
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

                for( int j = 0; j < k; j++ ) {
                    int e = nn_idx[j];
                    double r[NDIMS];

                    bary.transform(x, e, r);

                    // change this if-condition to (i == N) to debug the N-th marker
                    if (bary.is_inside(r)) {
                        ms.set_eta(i, r);
                        ms.set_elem(i, e);
                        ++elemmarkers[e][ms.get_mattype(i)];
                
                        found = true;
                        break;
                    }
                }

                if( found ) continue;

                /* not found */
                // Since no containing element has been found, delete this marker.
                removed_local.emplace_back(i);
            }

            #pragma omp critical
            removed_markers.insert(removed_markers.end(), removed_local.begin(), removed_local.end());
        }

        ms.remove_markers(removed_markers);

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_0(const Param& param, Variables &var,
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
                ++num_marker_in_elem;
            }
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_from_cpdf(const Param& param, Variables &var,
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
                ++num_marker_in_elem;
            }
        }
#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

    void replenish_markers_with_mattype_from_nn(const Param& param, Variables &var,
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
        KDTree kdtree(NDIMS, cloud);
        kdtree.buildIndex();
#ifdef USE_NPROF
        nvtxRangePop();
#endif
        MarkerSet &ms = *var.markersets[0];
        for (const auto& pair : unplenished_elems) {
            int e = pair.first;
            int num_marker_in_elem = pair.second;

            while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
                double eta[NODES_PER_ELEM];
                ms.random_eta(eta);

                double x[NDIMS] = {0};
                const int *conn = (*var.connectivity)[e];
                for (int d=0; d<NDIMS; d++) {
                    for (int i=0; i<NODES_PER_ELEM; i++) {
                        x[d] += (*var.coord)[ conn[i] ][d] * eta[i];
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
                const double fmelt = (*var.melt_fraction)[e];
    //            const double ti = ms.get_time(m);

                ms.append_marker(eta, e, mt, 0., 0., 0., 0., fmelt);
                if (DEBUG) {
                    std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
                }

                ++(*var.elemmarkers)[e][mt];
                ++num_marker_in_elem;
            }
        }

        delete points;

#ifdef USE_NPROF
        nvtxRangePop();
#endif
    }

} // anonymous namespace


void remap_markers(const Param& param, Variables &var, const array_t &old_coord,
                   const conn_t &old_connectivity)
{
#ifdef USE_NPROF
    nvtxRangePushA(__FUNCTION__);
#endif
    // Re-create elemmarkers
    delete var.elemmarkers;
    if (param.control.has_hydration_processes)
        delete var.hydrous_elemmarkers;
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
        PointCloud cloud(*centroid);
        KDTree kdtree(NDIMS, cloud);
        kdtree.buildIndex();
#ifdef USE_NPROF
        nvtxRangePop();
#endif

        find_markers_in_element(*var.markersets[0], *var.elemmarkers,
                                kdtree, bary, old_coord, old_connectivity);
        if (param.control.has_hydration_processes)
            find_markers_in_element(*var.markersets[var.hydrous_marker_index], *var.hydrous_elemmarkers,
                                    kdtree, bary, old_coord, old_connectivity);

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


