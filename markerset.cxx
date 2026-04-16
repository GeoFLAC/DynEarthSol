#include <cstring>
#include <iostream> // for std::cerr
#include <time.h> // for time()
#include <assert.h>
#include <unordered_map>
#include <set>
#include <numeric>
#include <random>

#include "constants.hpp"
#include "parameters.hpp"
#include "barycentric-fn.hpp"
#include "binaryio.hpp"
#include "markerset.hpp"
#include "mesh.hpp"
#include "geometry.hpp"
#include "utils.hpp"
#include "knn.hpp"
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

#ifdef HDF5
template
MarkerSet::MarkerSet(const Param& param, Variables& var, HDF5Input& bin_save, HDF5Input& bin_chkpt, const std::string& name);
#else
template
MarkerSet::MarkerSet(const Param& param, Variables& var, BinaryInput& bin_save, BinaryInput& bin_chkpt, const std::string& name);
#endif

template <class T>
MarkerSet::MarkerSet(const Param& param, Variables& var, T& bin_save, T& bin_chkpt, const std::string& name) :
    _name(name)
{
    // init from checkpoint file
    read_chkpt_file(var, bin_save, bin_chkpt);

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
    _genesis = new int_vec( max_markers );
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

void MarkerSet::random_eta_seed_surface( double *eta, int seed )
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double sum = 0;
    for(int n = 0; n < NDIMS; n++) {
        eta[n] = dist(gen);
        sum += eta[n];
    }
    for(int n = 0; n < NDIMS; n++)
        eta[n] /= sum;
}


void MarkerSet::random_eta_seed( ShapefnAccessor eta, int seed )
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
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    int nmarker = md.size();
    // Ensure sufficient array size
    while ( _nmarkers + nmarker > _reserved_space ) {
        // Resize the marker-related arrays if necessary.
        const int newsize = _reserved_space * over_alloc_ratio;
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
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void MarkerSet::append_marker_at_i(AppendMarkerData &md, int idx, int last_id)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    int m = idx;
    (*_eta)[m].copy_from(md.eta);
    (*_elem)[m] = md.elem;
    (*_mattype)[m] = md.mattype;
    (*_id)[m] = last_id;
    (*_time)[m] = md.time;
    (*_z)[m] = md.depth;
    (*_distance)[m] = md.distance;
    (*_slope)[m] = md.slope;
    (*_genesis)[m] = md.genesis;

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

template <typename T>
void MarkerSet::append_marker(T eta, int el, int mt , double ti, double z, double distance, double slope, int genesis )
{
    // Ensure sufficient array size
    if( _nmarkers == _reserved_space ) {
        // Resize the marker-related arrays if necessary.
        const int newsize = _nmarkers * over_alloc_ratio;
        resize( newsize );
    }

    int m = _nmarkers;
    (*_eta)[m].copy_from(eta);
    (*_elem)[m] = el;
    (*_mattype)[m] = mt;
    (*_id)[m] = _last_id;
    (*_time)[m] = ti;
    (*_z)[m] = z;
    (*_distance)[m] = distance;
    (*_slope)[m] = slope;
    (*_genesis)[m] = genesis;

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

template
void MarkerSet::append_marker<double*>(double* eta, int el, int mt , double ti, double z, double distance, double slope, int genesis );
template
void MarkerSet::append_marker<ConstShapefnAccessor>(ConstShapefnAccessor eta, int el, int mt , double ti, double z, double distance, double slope, int genesis );

// set marker  generated by deposit
void MarkerSet::set_surface_marker(const Param& param,const Variables& var, const double smallest_size, \
            const int mattype, double_vec& edvacc, int_vec2D& elemmarkers, int_vec2D& markers_in_elem)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

    double marker_dh_applied_ratio = 0.8;
    // double_vec &src_locs = *var.surfinfo.src_locs;

    // if (var.steps%1000==0) {
    //     double max_edvacc = 0.0;
    //     #pragma omp parallel for reduction(max:max_edvacc)
    //     for (int i=0; i<var.surfinfo.etop; i++) {
    //         if (edvacc[i] > max_edvacc) {
    //             max_edvacc = edvacc[i];
    //         }
    //     }

    //     printf("Surface accu mat #%d: ", mattype);
    //     for (int i=0; i<var.surfinfo.etop; i++) {
    //         int e = (*var.surfinfo.top_facet_elems)[i];
    //         if (edvacc[i] > max_edvacc - 1.) {
    //             double ratio = (param.markers.markers_per_element * edvacc[i]) / (*var.volume)[e];
    //             printf("elem[%6d] = %.1e, ratio = %5.1f%% ", e, edvacc[i], ratio*100.);
    //         }
    //     }
    //     printf("\n");
    // }

    #pragma omp parallel for default(none) shared(param, var, \
        edvacc, elemmarkers, markers_in_elem,marker_dh_applied_ratio) firstprivate(mattype)
    for (int i=0; i<var.surfinfo.etop; i++) {
        int e = (*var.surfinfo.top_facet_elems)[i];
        (*var.etmp_int)[i] = -1;

        if (edvacc[i] < 0.) {
            edvacc[i] = 0.;
            continue;
        }

        // if the ratio of increased Vol. and element Vol. up to certain level,
        // set a marker in the area between edvacc.
        double nmarkers;
        if (param.mesh.remeshing_option == 13 && param.mesh.meshing_elem_shape > 0) {
            nmarkers = param.markers.markers_per_element;
        } else {
            nmarkers = (*var.markers_in_elem)[e].size();
        }

        if ((nmarkers * edvacc[i]) < (*var.volume)[e]) continue;

        double eta[NDIMS], mcoord[NDIMS] = {0.};

        // random horizontal coordinate (2D:x or 3D:x-y) 
        random_eta_seed_surface(eta, e+var.steps);

        ConstArrayIndirectAccessor coord_surf = var.coord->view_const((*var.connectivity_surface)[i]);

        for (int j=0; j<NDIMS; j++) {
            for (int k=0; k<NDIMS; k++)
                mcoord[j] += coord_surf[k][j] * eta[k];
        }

        double base = compute_area_facet(coord_surf);

        // height of marker
        double dv_apply = (*var.volume)[e] / nmarkers;
        edvacc[i] -= dv_apply;

        double edh = dv_apply / base;
        mcoord[NDIMS-1] -= edh * marker_dh_applied_ratio;

        ConstArrayIndirectAccessor coord1 = var.coord->view_const((*var.connectivity)[e]);

        Barycentric_transformation bary(coord1, (*var.volume)[e]);

        double eta0[NDIMS];
        bary.transform(mcoord, 0, eta0);

        int elem_dest = e;

        if (!bary.is_inside(eta0)) {
            int inc = 0;
            remap_marker(var, mcoord, e, elem_dest, eta0, inc);

            if (!inc) {
                // msg += "... Success!\n";
                // printf("%s", msg.c_str());
            // } else {
                char buffer[200];
                sprintf(buffer, "  A generated marker (mat=%d) in element %7d is trying to remap in elements ",
                        mattype, e);
                std::string msg(buffer);
                printf("%s", msg.c_str());
                printf("... Surface marker generated fail!\n Coordinate: ");
                for (int j=0; j<NDIMS; j++) printf(" %f", mcoord[j]);
                printf("\neta: ");
                for (int j=0; j<NDIMS; j++) printf(" %d %f", j, eta0[j]); 
                printf("\n");

                std::exit(168);
            }
        }

        double eta1[NODES_PER_ELEM];

        eta1[NDIMS] = 1.;
        for (int j=0; j<NDIMS; j++) {
            eta1[j] = eta0[j];
            eta1[NDIMS] -= eta0[j];
        }
        #pragma omp atomic update
        ++elemmarkers[elem_dest][mattype];

        // the surface slope for marker
#ifdef THREED
        double slope; //, aspect_deg;
        ConstArrayAccessor p0 = coord_surf[0];
        ConstArrayAccessor p1 = coord_surf[1];
        ConstArrayAccessor p2 = coord_surf[2];
        double v1x = p1[0] - p0[0], v1y = p1[1] - p0[1], v1z = p1[2] - p0[2];
        double v2x = p2[0] - p0[0], v2y = p2[1] - p0[1], v2z = p2[2] - p0[2];

        // norm n = v1 x v2
        double nx = v1y * v2z - v1z * v2y;
        double ny = v1z * v2x - v1x * v2z;
        double nz = v1x * v2y - v1y * v2x;

        double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (norm == 0.0) {
            slope = 0.0;
            // aspect_deg = -1.0;
        } else {
            // slope: normal vector and z-axis angle
            slope = std::acos(nz / norm);
            // asmith: direction of maximum slope (relative to north, clockwise)
            // aspect_deg = std::atan2(ny, nx) * 180.0 / M_PI;
            // aspect_deg = 90.0 - aspect_deg;
            // if (aspect_deg < 0) aspect_deg += 360.0;
        }
#else
        double slope = ( coord_surf[0][NDIMS-1] - coord_surf[1][NDIMS-1] ) / ( coord_surf[0][0] - coord_surf[1][0] );
#endif
        double water_depth = var.surfinfo.base_level - mcoord[NDIMS-1];

        // double distance = 0.;
        // if (mcoord[0] >= src_locs[0] && mcoord[0] <= src_locs[1]) {
        //     double dx0 = mcoord[0] - src_locs[0];
        //     double dx1 = src_locs[1] - mcoord[0];
        //     if (slope < 0)
        //         distance = dx0;
        //     else if (slope > 0)
        //         distance = dx1;
        //     else
        //         distance =  std::min(dx0, dx1);
        // }

        ElemCacheAccessor marker_data = (*var.tmp_result)[i];

        (*var.etmp_int)[i] = elem_dest;
        for (int j=0; j<NODES_PER_ELEM; j++)
            marker_data[j] = eta1[j];
        marker_data[NODES_PER_ELEM] = water_depth;
        marker_data[NODES_PER_ELEM + 1] = 0.; // distance;
        marker_data[NODES_PER_ELEM + 2] = slope;
    }

    AMD_vec marker_data_all;

    for (int i=0; i<var.surfinfo.etop; i++) {
        if ((*var.etmp_int)[i] < 0) continue;

        int elem = (*var.etmp_int)[i];
        ConstElemCacheAccessor data = (*var.tmp_result)[i];
        double depth = data[NODES_PER_ELEM];
        double distance = data[NODES_PER_ELEM + 1];
        double slope = data[NODES_PER_ELEM + 2];
        int genesis = 2; // deposition

        marker_data_all.emplace_back(data, elem, mattype, 
                var.time / YEAR2SEC, depth, distance, slope, genesis);
    }

    append_markers(marker_data_all);

    int nnew = marker_data_all.size();
    for (int i=0; i<nnew; ++i)
        markers_in_elem[marker_data_all[i].elem].push_back(_nmarkers-nnew+i);

    // if (nnew > 0)
    //     printf("Set surface markers: %d (mat: %d)\n", nnew, mattype);

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif

}


void MarkerSet::remap_marker(const Variables &var, const double *m_coord, \
                    const int e, int& new_elem, double *new_eta, int& inc)
{
    const double *coord[NODES_PER_ELEM];
    double eta[NODES_PER_ELEM];
    int_vec searched(var.nelem,0);
    searched[e]=1;
    ConstConnAccessor conn = (*var.connectivity)[e];

//    std::cout << "Try to remap in ";
    for (int i = 0; i < NODES_PER_ELEM; i++) {
        int n = conn[i];        
        for(auto ee = (*var.support)[n].begin(); ee < (*var.support)[n].end(); ++ee) {
            if (searched[*ee]) continue;
            searched[*ee]=1;

            ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[*ee]);

            double volume = compute_volume(coord);
            Barycentric_transformation bary(coord, volume);
            bary.transform(m_coord,0,eta);

            if (bary.is_inside(eta)) {
                for (int j=0; j<NDIMS; j++)
                    new_eta[j] = eta[j];

                new_elem = *ee;
                inc = 1;
                return;
            }
        }
    }
    inc = 0;
}

void MarkerSet::append_random_marker_in_elem( int el, int mt, int genesis)
{
    double eta[NODES_PER_ELEM];
    random_eta(eta);
    append_marker(eta, el, mt, 0., 0., 0., 0., genesis);
}

void MarkerSet::append_random_marker_in_elem( int el, int mt, double time, int genesis)
{
    double eta[NODES_PER_ELEM];
    random_eta(eta);
    append_marker(eta, el, mt, time, 0., 0., 0., genesis);
}

void MarkerSet::random_markers( const Param& param, Variables &var, int genesis )
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
            append_marker(eta, e, mt, 0., 0., 0., 0., genesis);
            ++(*var.elemmarkers)[e][mt];
            (*var.markers_in_elem)[e].push_back(_nmarkers-1);
        }
}


void MarkerSet::regularly_spaced_markers( const Param& param, Variables &var, int genesis )
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
    array_t centroid(var.nelem);
    elem_center(*var.coord, *var.connectivity, centroid);

    PointCloud cloud(centroid);
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
                append_marker(eta, e, mt, 0., 0., 0., 0., genesis);
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
            ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[elem]);
            for(int i=0; i<NDIMS; i++) {
                for(int j=0; j<NODES_PER_ELEM; j++)
                    p[i] += coord[j][i] * eta[j];
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


int MarkerSet::custom_initial_mattype(const Param& param, const Variables &var,
                                       int elem, const double eta[NODES_PER_ELEM],
                                       const double *x )
{
    /* User defined function */
    int mt = 0;

    return mt;
}

void MarkerSet::remove_markers(const Param& param, const Variables &var, int_vec& markers, int_vec2D& markers_in_elem)
{
    if (markers.empty()) return;

#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
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
        #pragma acc parallel loop gang vector async
#endif
        for (int i = 0; i < a_out.size(); i++) {
            int_vec& emarkers = markers_in_elem[(*_elem)[b_out[i]]];
            auto it = std::find(emarkers.begin(), emarkers.end(), b_out[i]);
            emarkers[it - emarkers.begin()] = a_out[i];
            remove_marker_data(a_out[i],b_out[i]);
        }
    }

    // printf("    Removed %d markers from markerset (", n);

    // int *ndelete = var.etmp_int->data();
    // for (int i=0; i<param.mat.nmat; ++i) {
    //     ndelete[i]  = 0;
    // }
    // for (int i=0; i<markers.size();++i) {
    //     ndelete[(*_mattype)[markers[i]]]++;
    // }
    // for (int i=0; i<param.mat.nmat; ++i) {
    //     if (ndelete[i] > 0) {
    //         printf("mat #%d: %d; ", i, ndelete[i]);
    //     }
    // }
    // printf(")\n");

    #pragma acc wait

    _nmarkers -= n;

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void MarkerSet::remove_marker(int i)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // Replace marker i by the last marker.
    --_nmarkers;
    remove_marker_data(i, _nmarkers);
#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void MarkerSet::remove_marker_data(int is, int ie)
{
    // Replace marker i by the target marker ie.
    (*_eta)[is] = (*_eta)[ie];
    (*_id)[is] = (*_id)[ie];
    (*_elem)[is] = (*_elem)[ie];
    (*_mattype)[is] = (*_mattype)[ie];
    (*_time)[is] = (*_time)[ie];
    (*_z)[is] = (*_z)[ie];
    (*_distance)[is] = (*_distance)[ie];
    (*_slope)[is] = (*_slope)[ie];
    (*_genesis)[is] = (*_genesis)[ie];
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

        _eta->resize( newsize );
        _elem->resize( newsize );
        _mattype->resize( newsize );
        _id->resize( newsize );
        _time->resize( newsize );
        _z->resize( newsize );
        _distance->resize( newsize );
        _slope->resize( newsize );
        _genesis->resize( newsize );
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

#ifdef HDF5
template
void MarkerSet::write_chkpt_file(HDF5Output &bin) const;
#else
template
void MarkerSet::write_chkpt_file(BinaryOutput &bin) const;
#endif

template <class T>
void MarkerSet::write_chkpt_file(T &bin) const
{
#ifdef HDF5
    bin.write_scalar(_last_id, _name + ".last_id");
    bin.write_scalar(_reserved_space, _name + ".reserved_space");
#else
    int_vec itmp(3);
    itmp[0] = _nmarkers;
    itmp[1] = _last_id;
    itmp[2] = _reserved_space;
    bin.write_array(itmp, (_name + " size").c_str(), itmp.size());
#endif
    bin.write_array(*_genesis, (_name + ".genesis").c_str(), _nmarkers);
}

#ifdef HDF5
template
void MarkerSet::read_chkpt_file(Variables &var, HDF5Input &bin_save, HDF5Input &bin_chkpt);
#else
template
void MarkerSet::read_chkpt_file(Variables &var, BinaryInput &bin_save, BinaryInput &bin_chkpt);
#endif

template <class T>
void MarkerSet::read_chkpt_file(Variables &var, T &bin_save, T &bin_chkpt)
{
#ifdef HDF5
    bin_save.read_scaler(_nmarkers, _name + ".nmarkers");
    bin_chkpt.read_scaler(_last_id, _name + ".last_id");
    bin_chkpt.read_scaler(_reserved_space, _name + ".reserved_space");
#else
    int_vec itmp(3);
    bin_chkpt.read_array(itmp, (_name + " size").c_str());
    _nmarkers = itmp[0];
    _last_id = itmp[1];
    _reserved_space = itmp[2];
#endif
    allocate_markerdata(_reserved_space);

    if (_nmarkers != 0) {
        bin_save.read_array(*_eta, (_name + ".eta").c_str(), _nmarkers);
        bin_save.read_array(*_elem, (_name + ".elem").c_str(), _nmarkers);
        bin_save.read_array(*_mattype, (_name + ".mattype").c_str(), _nmarkers);
        bin_save.read_array(*_id, (_name + ".id").c_str(), _nmarkers);
        bin_save.read_array(*_time, (_name + ".time").c_str(), _nmarkers);
        bin_save.read_array(*_z, (_name + ".z").c_str(), _nmarkers);
        bin_save.read_array(*_distance, (_name + ".distance").c_str(), _nmarkers);
        bin_save.read_array(*_slope, (_name + ".slope").c_str(), _nmarkers);

        bin_chkpt.read_array(*_genesis, (_name + ".genesis").c_str(), _nmarkers);
    }
}

#ifdef HDF5
template
void MarkerSet::write_save_file(const Variables &var, HDF5Output &bin) const;
#else
template
void MarkerSet::write_save_file(const Variables &var, BinaryOutput &bin) const;
#endif

template <class T>
void MarkerSet::write_save_file(const Variables &var, T &bin) const
{
#ifdef NPROF_DETAIL
    nvtxRangePush("write markersets");
#endif

#ifndef HDF5
    int_vec itmp(1);
    itmp[0] = _nmarkers;
    bin.write_array(itmp, (_name + " size").c_str(), itmp.size());

    array_t mcoord(_nmarkers);
    calculate_marker_coord(var, mcoord); // coordinate of markers
    bin.write_array(mcoord, (_name + ".coord").c_str(), _nmarkers);
#endif
    bin.write_array(*_eta, (_name + ".eta").c_str(), _nmarkers);
    bin.write_array(*_elem, (_name + ".elem").c_str(), _nmarkers);
    bin.write_array(*_mattype, (_name + ".mattype").c_str(), _nmarkers);
    bin.write_array(*_id, (_name + ".id").c_str(), _nmarkers);
    bin.write_array(*_time, (_name + ".time").c_str(), _nmarkers);
    bin.write_array(*_z, (_name + ".z").c_str(), _nmarkers);
    bin.write_array(*_distance, (_name + ".distance").c_str(), _nmarkers);
    bin.write_array(*_slope, (_name + ".slope").c_str(), _nmarkers);

    // temporary show results
    bin.write_array(*_genesis, (_name + ".genesis").c_str(), _nmarkers);

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void MarkerSet::get_ZPT(const Param& param, const Variables& var, int m, double &Z, double &P, double &T) const {
        // Get depth and temperature at the marker
        Z = T = 0;
        ConstShapefnAccessor eta = (*_eta)[m];
        ConstConnAccessor conn = (*var.connectivity)[(*_elem)[m]];
        for (int i=0; i<NODES_PER_ELEM; ++i) {
            Z += (*var.coord)[conn[i]][NDIMS-1] * eta[i];
            T += (*var.temperature)[conn[i]] * eta[i];
        }

        // Get pressure, which is constant in the element
        // P = - trace((*var.stress)[e]) / NDIMS;
        P = ref_pressure(param, Z);    
}


void MarkerSet::calculate_marker_coord(const Variables &var, array_t &points) const {
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif

#ifndef ACC
    #pragma omp parallel for default(none) shared(var, points)
#endif
    #pragma acc parallel loop gang vector async
    for (int n=0; n<_nmarkers; n++) {
        ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[(*_elem)[n]]);

        for(int d=0; d<NDIMS; d++) {
            double sum = 0;
            for(int k=0; k<NODES_PER_ELEM; k++)
                sum += coord[k][d] * (*_eta)[n][k];
            points[n][d] = sum;
        }
    }

    #pragma acc wait

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}

namespace {

    template <class T>
    void find_markers_in_element(const Param &param, const Variables& var, MarkerSet& ms,
                                 T& elemmarkers, int_vec2D& markers_in_elem,
                                 KNN& kdtree, const Barycentric_transformation &bary,
                                 const array_t& old_coord, const conn_t &old_connectivity)
    {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif
        const int k_neig = std::min((std::size_t) 20, old_connectivity.size());  // how many nearest neighbors to search?

        int nmarkers = ms.get_nmarkers();

        int markers_per_block_max = std::max(kdtree.max_batch_size(k_neig), 1);
        int markers_per_block = std::min(markers_per_block_max, nmarkers);
        int nblocks = (nmarkers + markers_per_block - 1) / markers_per_block;

        {
            char knn_be[64]; kdtree.backend_str(knn_be, sizeof(knn_be));
            char mem_str[32] = "";
            size_t qm = KNN::query_mem_bytes(markers_per_block, k_neig);
            if (qm) snprintf(mem_str, sizeof(mem_str), ", %6.1f MB", qm / 1048576.0);

            printf("    Marker remapping: knn: %10s pts (%s) for %11s q w/ %2d k%s",
                   format_with_commas((unsigned long)kdtree.get_npoints()).c_str(), knn_be,
                   format_with_commas((unsigned long)nmarkers).c_str(), k_neig, mem_str);

            if (nblocks > 1) {
                printf("/blk (%d x %sq)", nblocks,
                        format_with_commas((unsigned long)markers_per_block).c_str());
                fflush(stdout);
            } else {
                printf("\n");
            }
        }

        array_t queries(markers_per_block);

        for (int b=0; b<nblocks; ++b) {
            int start = b * markers_per_block;
            int end = std::min(start + markers_per_block, nmarkers);

            if (nblocks > 1) {printf(" %d", b+1); fflush(stdout);}

            queries.resize(end - start);

#ifndef ACC
            #pragma omp parallel for default(none) shared(ms, queries, old_coord, \
                old_connectivity, start, end)
#endif
            #pragma acc parallel loop gang vector
            for (int i = start; i < end; i++) {
                int idx_q = i - start;
                // 1. Get physical coordinates, x, of an old marker.
                ConstConnAccessor conn = old_connectivity[ms.get_elem(i)];
                ConstShapefnAccessor eta = ms.get_eta(i);
                for (int j = 0; j < NDIMS; j++) {
                    queries[idx_q][j] = 0.;
                    for (int l = 0; l < NODES_PER_ELEM; l++) {
                        queries[idx_q][j] += eta[l] * old_coord[ conn[l] ][j];
                    }
                }
            }

            neighbor* neighbors = kdtree.search(queries, (end - start), k_neig, false);

            // Loop over all the old markers and identify a containing element in the new mesh.
#ifndef ACC
            #pragma omp parallel for default(none) shared(param, ms, bary, old_coord, old_connectivity, \
                nmarkers, queries, neighbors, start, end) firstprivate(k_neig)
#endif
            #pragma acc parallel loop gang vector deviceptr(neighbors)
            for (int i = start; i < end; i++) {
                int idx_q = i - start;
                bool found = false;

                // 2. look for nearby elements.
                neighbor* nn_idx = neighbors + idx_q * k_neig;

                for( int j = 0; j < k_neig; j++ ) {
                    int e = nn_idx[j].idx;

                    if (e < 0) break; // no more neighbors

                    double r[NDIMS];

                    bary.transform(queries[idx_q], e, r);

                    // change this if-condition to (i == N) to debug the N-th marker
                    if (bary.is_inside(r)) {
                        ms.set_eta(i, r);
                        ms.set_elem(i, e);            
                        found = true;
                        break;
                    }
                }

                if( found ) continue;

                /* not found */
                // Since no containing element has been found, delete this marker.
                ms.set_elem(i, -1.); // mark this marker for removal
            }
        }

        if (nblocks > 1) printf("\n");

        int_vec removed_markers;

        #pragma acc wait

        for (int i = 0; i < nmarkers; i++) {
            int e = ms.get_elem(i);
            if (e >= 0)
                markers_in_elem[e].push_back(i);
            else
                removed_markers.push_back(i);
        }

#ifndef ACC
        #pragma omp parallel for default(none) shared(var, ms, markers_in_elem, elemmarkers)
#endif
        #pragma acc parallel loop gang vector async
        for (int e = 0; e < var.nelem; e++) {
            int_vec &markers = markers_in_elem[e];
            for (int i = 0; i < markers.size(); i++) {
                int m = markers[i];
                if (ms.get_elem(m) >= 0) {
                    // This marker is not removed, so we need to update the element markers.
                    elemmarkers[e][ms.get_mattype(m)]++;
                }
            }
        }

        #pragma acc wait

        ms.remove_markers(param, var, removed_markers, markers_in_elem);

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_0(const Param& param, const Variables &var,
                                          int_pair_vec &unplenished_elems, int genesis)
    {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif
        for (const auto& pair : unplenished_elems) {
            int e = pair.first;
            int num_marker_in_elem = pair.second;

            while( num_marker_in_elem < param.markers.min_num_markers_in_element ) {
                const int mt = 0;
                var.markersets[0]->append_random_marker_in_elem(e, mt, genesis);
                if (DEBUG) {
                    std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
                }

                ++(*var.elemmarkers)[e][mt];
                (*var.markers_in_elem)[e].push_back(var.markersets[0]->get_nmarkers()-1);
                ++num_marker_in_elem;
            }
        }
#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }


    void replenish_markers_with_mattype_from_cpdf(const Param& param, const Variables &var,
                                                  int_pair_vec &unplenished_elems, int genesis)
    {
#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
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
                var.markersets[0]->append_random_marker_in_elem(e, mt, genesis);
                if (DEBUG) {
                    std::cout << "Add marker with mattype " << mt << " in element " << e << '\n';
                }

                ++(*var.elemmarkers)[e][mt];
                (*var.markers_in_elem)[e].push_back(var.markersets[0]->get_nmarkers()-1);
                ++num_marker_in_elem;
            }
        }
#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }

    double interpolate_nd(const double *p1, double v1, const double *p2, double v2, ConstArrayAccessor pq)
    {
        double dir[NDIMS], rel[NDIMS];
        double len2 = 0.0;

        for (int i = 0; i < NDIMS; ++i) {
            dir[i] = p2[i] - p1[i];
            rel[i] = pq[i] - p1[i];
            len2 += dir[i] * dir[i];
        }

        if (len2 == 0.0) return (v1 + v2) / 2.;

        double t = 0.;
        for (int i = 0; i < NDIMS; ++i) t += dir[i]*rel[i];

        return v1 + (v2 - v1) * t / len2;
    }

    // Interpolate deposition time and set genesis=4 for a new sediment marker placed in element e.
    // Averages existing sediment markers in two groups (reworked: genesis 2|4 → ind=0, others → ind=1)
    // and picks the first non-empty group for interpolation.  If no suitable reference exists,
    // ti and ge are left unchanged.
    void compute_sed_marker_time_genesis(
        const Param& param, const Variables& var, const MarkerSet& ms,
        int e, int e_local, const ElemMarkerInfo* emi_ptr,
        ConstArrayAccessor x, double& ti, int& ge)
    {
        if (e_local == -1 || emi_ptr[e_local].nmarker == 0) return;

        ElemMarkerInfo emi_base[2];
        for (int m : (*var.markers_in_elem)[e]) {
            if (ms.get_mattype(m) != param.mat.mattype_sed) continue;
            int ind = (ms.get_genesis(m) == 2 || ms.get_genesis(m) == 4) ? 0 : 1;
            emi_base[ind].value += ms.get_time(m);
            ConstShapefnAccessor eta = ms.get_eta(m);
            ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[e]);
            for (int d = 0; d < NDIMS; ++d)
                for (int kk = 0; kk < NODES_PER_ELEM; ++kk)
                    emi_base[ind].coord[d] += coord[kk][d] * eta[kk];
            emi_base[ind].nmarker++;
        }

        for (int k = 0; k < 2; ++k) {
            if (emi_base[k].nmarker == 0) continue;
            if (emi_base[k].nmarker > 1) {
                emi_base[k].value /= emi_base[k].nmarker;
                for (int d = 0; d < NDIMS; ++d)
                    emi_base[k].coord[d] /= emi_base[k].nmarker;
            }
            ti = std::max(interpolate_nd(emi_base[k].coord, emi_base[k].value,
                                         emi_ptr[e_local].coord, emi_ptr[e_local].value, x), 0.0);
            ge = 4; // sediment with interpolation
            break;
        }
    }

    void replenish_markers_with_mattype_from_nn(const Param& param, const Variables &var,
            int_pair_vec &unplenished_elems, int genesis, bool is_surface = false, const EMI_vec& emi = EMI_vec())
    {
        if (unplenished_elems.empty()) return;

#ifdef NPROF_DETAIL
        nvtxRangePush(__FUNCTION__);
#endif

        int k_neig = 1, nnew = 0;
        int ne = unplenished_elems.size();
        int_vec nneed_mk(ne), mk_start(ne, 0);
        MarkerSet &ms = *var.markersets[0];

#ifndef ACC
        #pragma omp parallel for default(none) shared(param, unplenished_elems, ne, nneed_mk)
#endif
        #pragma acc parallel loop gang vector
        for (int i=0; i<ne; ++i)
            nneed_mk[i] = param.markers.min_num_markers_in_element - unplenished_elems[i].second;

        for (int i=1; i<ne; ++i)
            mk_start[i] = mk_start[i-1] + nneed_mk[i-1];
        nnew = mk_start[ne-1] + nneed_mk[ne-1];

        array_t queries(nnew);
        shapefn etas(nnew);

        // normal random cannot be parallelized in acc
        #pragma omp parallel for default(none) \
                shared(ms, unplenished_elems, nneed_mk, etas, ne, var, mk_start)
        for (int i=0; i<ne; ++i) {
            int e = unplenished_elems[i].first;
            int num_marker_in_elem = unplenished_elems[i].second;

            for (int j=0; j<nneed_mk[i]; ++j) {
                ShapefnAccessor eta = etas[mk_start[i] + j];
                ms.random_eta_seed(eta, e + num_marker_in_elem + var.steps);

                num_marker_in_elem++;
            }
        }

#ifndef ACC
        #pragma omp parallel for default(none) \
                shared(var, unplenished_elems, nneed_mk, etas, queries, ne, mk_start)
#endif
        #pragma acc parallel loop gang vector async
        for (int i=0; i<ne; ++i) {
            int e = unplenished_elems[i].first;
            for (int j=0; j<nneed_mk[i]; ++j) {
                int m_idx = mk_start[i] + j;
                ArrayAccessor x = queries[m_idx];
                ConstArrayIndirectAccessor coord = var.coord->view_const((*var.connectivity)[e]);
                for (int d=0; d<NDIMS; d++) {
                    x[d] = 0.;
                    for (int ii=0; ii<NODES_PER_ELEM; ii++) {
                        x[d] += coord[ii][d] * etas[m_idx][ii];
                    }
                }
            }
        }

        AMD_vec marker_data_all(nnew);
        array_t points(ms.get_nmarkers());

        ms.calculate_marker_coord(var, points); // coordinate of markers

#ifdef ACC
        array_t point_tmp(1);
        PointCloud cloud(point_tmp);
#else
        PointCloud cloud(points);
#endif

#ifdef NPROF_DETAIL
        nvtxRangePush("create kdtree for markers");
#endif
        NANOKDTree nano_kdtree(NDIMS, cloud);

        int MAX_EXPECTED_N = -1;
        // Estimate the maximum number of points by the end of the iteration cycle.
        // We add a 20% + 10,000 baseline padding to prevent reallocations.
        if (is_surface) {
            MAX_EXPECTED_N = ms.get_nmarkers() * 1.2 + 10000;
            if (MAX_EXPECTED_N < ms.get_nmarkers() + nnew) {
                MAX_EXPECTED_N = ms.get_nmarkers() + nnew + 10000;
            }
        }
        
        KNN kdtree(param, points, nano_kdtree, false, MAX_EXPECTED_N);

        if (!is_surface && nnew > 0) {
            char knn_be[64]; kdtree.backend_str(knn_be, sizeof(knn_be));
            char mem_str[32] = "";
            size_t qm = KNN::query_mem_bytes(nnew, k_neig);
            if (qm) snprintf(mem_str, sizeof(mem_str), ", %6.1f MB", qm / 1048576.0);
            printf("    Marker replenish: knn: %10s pts (%s) for %11s q w/ %2d k%s\n",
                    format_with_commas((unsigned long)kdtree.get_npoints()).c_str(), knn_be,
                    format_with_commas((unsigned long)nnew).c_str(), k_neig, mem_str);
        }

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif

        neighbor* neighbors = kdtree.search(queries, nnew, k_neig, false);

        #pragma acc wait
        const ElemMarkerInfo* emi_ptr = emi.data();
        const int* top_elems_ptr = var.top_elems->data();
        const int ntop_elems = var.ntop_elems;
        bool is_no_nn = false;

#ifndef ACC
        #pragma omp parallel for default(none) shared(param, var, unplenished_elems, \
                ms, marker_data_all, ne, is_surface, emi_ptr, top_elems_ptr, ntop_elems, \
                genesis, nneed_mk, mk_start, neighbors, queries, etas) reduction(||:is_no_nn)
#endif
        #pragma acc parallel loop gang vector deviceptr(neighbors) reduction(||:is_no_nn)
        for (int i=0; i<ne; ++i) {
            int e = unplenished_elems[i].first;
            int start = mk_start[i];

            for (int j=0; j<nneed_mk[i]; ++j) {
                int m_idx = start + j;
                int m = neighbors[m_idx].idx;
                if (m < 0) {
                    is_no_nn = true;
                    break;
                }

                ConstArrayAccessor x = queries[m_idx];

                int mt = ms.get_mattype(m);
                double ti = ms.get_time(m);
                int ge = genesis;

                if (is_surface && mt == param.mat.mattype_sed) {
                    int e_local = binary_search_index(top_elems_ptr, ntop_elems, e);
                    compute_sed_marker_time_genesis(param, var, ms, e, e_local, emi_ptr, x, ti, ge);
                }

                marker_data_all[m_idx] = AppendMarkerData(etas[m_idx], e, mt, ti, 0., 0., 0., ge);
                ++(*var.elemmarkers)[e][mt];
            }
        }

        if (is_no_nn) {
            std::cerr << "Error: no nearest neighbor found for some new markers.\n";
            std::exit(168);
        }

        // Append new markers to the end of the marker set.
        ms.append_markers(marker_data_all);

        int nmarkers = ms.get_nmarkers();
        for (int i=0; i<nnew; ++i)
            (*var.markers_in_elem)[marker_data_all[i].elem].push_back(nmarkers-nnew+i);

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif
    }

} // anonymous namespace


void MarkerSet::check_marker_elem_consistency(const Variables &var) const
{
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif
    #pragma acc serial
    int ncount = 0, is_error = 0;
#ifndef ACC
    #pragma omp parallel for reduction(+:ncount,is_error) default(none) shared(var,std::cerr)
#endif
    #pragma acc parallel loop gang vector reduction(+:ncount,is_error)
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

#ifdef NPROF
    nvtxRangePop();
#endif
}


// surface processes correcttion of marker
void MarkerSet::correct_surface_marker(const Param &param, const Variables& var, const double_vec& dhacc, int_vec2D &elemmarkers, int_vec2D &markers_in_elem)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
#endif
    // correct surface marker.
    Barycentric_transformation bary(*var.top_elems, *var.coord, *var.connectivity, *var.volume);

    array_t coord0s(var.ntop_elems*NODES_PER_ELEM, 0.);

    int_vec delete_marker;
    int nchange = 0;
#ifndef ACC
    #pragma omp parallel for default(none) shared(var,coord0s,dhacc)
#endif
    #pragma acc parallel loop gang vector async
    for (int i=0;i<var.ntop_elems;i++) {
        ConstConnAccessor tnodes = (*var.connectivity)[(*var.top_elems)[i]];
        ConstArrayIndirectAccessor coord = var.coord->view_const(tnodes);

        ArrayAccessor c00 = coord0s[i*NODES_PER_ELEM];
        ArrayAccessor c01 = coord0s[i*NODES_PER_ELEM+1];
        ArrayAccessor c02 = coord0s[i*NODES_PER_ELEM+2];
#ifdef THREED
        ArrayAccessor c03 = coord0s[i*NODES_PER_ELEM+3];
#endif

        // restore the reference node locations before deposition/erosion 
        c00[0] = coord[0][0];
        c01[0] = coord[1][0];
        c02[0] = coord[2][0];

        c00[NDIMS-1] = coord[0][NDIMS-1] - dhacc[tnodes[0]];
        c01[NDIMS-1] = coord[1][NDIMS-1] - dhacc[tnodes[1]];
        c02[NDIMS-1] = coord[2][NDIMS-1] - dhacc[tnodes[2]];
#ifdef THREED
        c00[1] = coord[0][1];
        c01[1] = coord[1][1];
        c02[1] = coord[2][1];

        c03[0] = coord[3][0];
        c03[1] = coord[3][1];
        c03[2] = coord[3][2] - dhacc[tnodes[3]];
#endif
    }

#ifndef ACC
    #pragma omp parallel for default(none) shared(var,coord0s,bary,markers_in_elem) reduction(+:nchange)
#endif
    #pragma acc parallel loop gang vector reduction(+:nchange) async
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

        EMI_vec markers_in_elem_info(var.ntop_elems);
        MU_vec updates;

        #pragma omp parallel default(none) shared(param, var, coord0s, markers_in_elem, \
                updates, elemmarkers, markers_in_elem_info)
        {
            MU_vec local_updates;

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

                    local_updates.emplace_back(m, e, new_elem, inc);

                    if (inc) {
                        set_eta(m, new_eta);
                        set_elem(m, new_elem);
                        #pragma omp atomic update
                        ++elemmarkers[new_elem][mat];
                    } else {
                        // record time information on deleted sediment marker
                        if (mat == param.mat.mattype_sed) {
                            for (int d=0; d<NDIMS; ++d) {
                                markers_in_elem_info[i].coord[d] += m_coord[d];
                            }
                            markers_in_elem_info[i].value += (*_time)[m];
                            markers_in_elem_info[i].nmarker++;
                        }
                    }
                    --elemmarkers[e][mat];
                }
                int n = markers_in_elem_info[i].nmarker;
                if (n > 1) {
                    markers_in_elem_info[i].value /= n;
                    for (int d=0; d<NDIMS; ++d) {
                        markers_in_elem_info[i].coord[d] /= n;
                    }
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
        remove_markers(param, var, delete_marker, markers_in_elem);

        int_pair_vec unplenished_elems;
        for (int i=0; i<var.ntop_elems; i++) {
            int e = (*var.top_elems)[i];
            int  nmarkers = markers_in_elem[e].size();

            if (nmarkers < param.markers.min_num_markers_in_element)
                unplenished_elems.emplace_back(e, nmarkers);
        }

        // if (unplenished_elems.size() > 0) {
        //     printf("replenish markers in %d elements.\n", (int)unplenished_elems.size());
        //     for (int i=0; i<(int)unplenished_elems.size(); i++) {
        //         int e = unplenished_elems[i].first;
        //         int nmarkers = unplenished_elems[i].second;
        //         printf("  Element %d has %d markers.\n", e, nmarkers);
        //     }
        // }

        switch (param.markers.replenishment_option) {
        case 0:
            replenish_markers_with_mattype_0(param, var, unplenished_elems, 3);
            break;
        case 1:
            replenish_markers_with_mattype_from_cpdf(param, var, unplenished_elems, 3);
            break;
        case 2:
            replenish_markers_with_mattype_from_nn(param, var, unplenished_elems, 3, true, markers_in_elem_info);
            break;
        default:
            std::cerr << "Error: unknown markers.replenishment_option: " << param.markers.replenishment_option << '\n';
            std::exit(1);
        }

        // for (int i=0; i<var.ntop_elems; i++) {
        //     int e = (*var.top_elems)[i];
        //     int nemarker0 = std::accumulate((*var.elemmarkers)[e].begin(), (*var.elemmarkers)[e].end(), 0);
        //     int  nmarkers = (*var.markers_in_elem)[e].size();
        //     if (nemarker0 != nmarkers) {
        //         std::cerr << "Error: number of markers in element " << e
        //                 << " does not match number of elemmarkers: "
        //                 << nmarkers << " vs. " << nemarker0 << '\n';
        //         std::exit(1);
        //     }
        // }
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


void remap_markers(const Param& param, Variables &var, const array_t &old_coord,
                   const conn_t &old_connectivity)
{
#ifdef NPROF_DETAIL
    nvtxRangePush(__FUNCTION__);
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
#ifdef NPROF_DETAIL
        nvtxRangePush("create kdtree for new elements");
#endif

        array_t points(var.nelem);
        elem_center(*var.coord, *var.connectivity, points); // centroid of elements

#ifdef ACC
        array_t point_tmp(1);
        PointCloud cloud(point_tmp);
#else
        PointCloud cloud(points);
#endif

        NANOKDTree nano_kdtree(NDIMS, cloud);
        KNN kdtree(param, points, nano_kdtree, false);

#ifdef NPROF_DETAIL
        nvtxRangePop();
#endif

        find_markers_in_element(param, var, *var.markersets[0], *var.elemmarkers, *var.markers_in_elem,
                                kdtree, bary, old_coord, old_connectivity);
        if (param.control.has_hydration_processes)
            find_markers_in_element(param, var, *var.markersets[var.hydrous_marker_index], *var.hydrous_elemmarkers, *var.hydrous_markers_in_elem,
                                    kdtree, bary, old_coord, old_connectivity);

        #pragma acc wait
    }

    // If any new element has too few markers, generate markers in them.
#ifdef NPROF_DETAIL
    nvtxRangePush("find unplenished elements");
#endif


    int nunplenished = 0;

#ifndef ACC
    #pragma omp parallel default(none) shared(param, var) reduction(+:nunplenished)
#endif
    #pragma acc parallel loop gang vector reduction(+:nunplenished)
    for (int e = 0; e < var.nelem; e++) {
        int num_marker_in_elem = (*var.markers_in_elem)[e].size();
        if (num_marker_in_elem < param.markers.min_num_markers_in_element) {
            (*var.etmp_int)[e] = num_marker_in_elem;
            ++nunplenished;            
        } else {
            (*var.etmp_int)[e] = -1;
        }
    }

    // unplenish markers
    int_pair_vec unplenished_elems;
    unplenished_elems.reserve(nunplenished);

    #pragma omp parallel default(none) shared(var, unplenished_elems)
    {
        int_pair_vec local_pairs;

        #pragma omp for nowait
        for (int e = 0; e < var.nelem; e++) {
            if ((*var.etmp_int)[e] >= 0) {
                local_pairs.emplace_back(e, (*var.etmp_int)[e]);
            }
        }
        #pragma omp critical
        unplenished_elems.insert(unplenished_elems.end(), 
                                local_pairs.begin(),
                                local_pairs.end());
    }

    std::sort(unplenished_elems.begin(), unplenished_elems.end(),
        [](const int_pair &a, const int_pair &b) {
            return a.first < b.first;
        });

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif

    switch (param.markers.replenishment_option) {
    case 0:
        replenish_markers_with_mattype_0(param, var, unplenished_elems, 1);
        break;
    case 1:
        replenish_markers_with_mattype_from_cpdf(param, var, unplenished_elems, 1);
        break;
    case 2:
        replenish_markers_with_mattype_from_nn(param, var, unplenished_elems, 1);
        break;
    default:
        std::cerr << "Error: unknown markers.replenishment_option: " << param.markers.replenishment_option << '\n';
        std::exit(1);
    }

#ifdef NPROF_DETAIL
    nvtxRangePop();
#endif
}


namespace {

    Barycentric_transformation* get_bary_from_cache(std::unordered_map<int, Barycentric_transformation*> &cache,
                                                    int el, const array_t &coordinate, conn_t::ConstAccessor conn,
                                                    double_vec &volume)
    {
        Barycentric_transformation *bary;
        auto search = cache.find(el);
        if (search == cache.end()) {
            ConstArrayIndirectAccessor coord = coordinate.view_const(conn);
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
        ConstConnAccessor conn = (*var.connectivity)[el];
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
                    ConstConnAccessor conn2 = (*var.connectivity)[ee];
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


