#ifndef DYNEARTHSOL3D_MARKERSET_HPP
#define DYNEARTHSOL3D_MARKERSET_HPP

#include <string>
#include "barycentric-fn.hpp"

class MarkerSet
{

public:
    MarkerSet( const std::string& name );
    MarkerSet( const Param& param, Variables& var, const std::string& name );
    template <class T>
    MarkerSet( const Param& param, Variables& var, T& bin_save, T& bin_chkpt, const std::string& name );
    ~MarkerSet()
    {
        delete _slope;
        delete _distance;
        delete _z;
        delete _time;
        delete _id;
        delete _eta; 
        delete _elem;
        delete _mattype;
        delete _genesis;
        delete _tmp;
    }

    static void random_eta( double* ); // class method
    static void random_eta_seed(double*, int);
    static void random_eta_seed_surface(double*, int);
    void check_marker_elem_consistency(const Variables &var) const;
    void correct_surface_marker(const Param& param, const Variables& var, const double_vec& dhacc, int_vec2D &elemmarkers, int_vec2D &markers_in_elem);
    void set_surface_marker(const Param& param ,const Variables& var, const double smallest_size, \
                        const int mattype_sed, double_vec& edvacc, int_vec2D& elemmarkers, int_vec2D& markers_in_elem);
    void remap_marker(const Variables &var, const double *m_coord, const int e, int &new_elem, double *new_eta, int &inc);
    void append_random_marker_in_elem( int el, int mt, int genesis);
    void append_random_marker_in_elem( int el, int mt, double time, int genesis);
    void append_marker(const double *eta, int el, int mt, double time, double depth, double distance, double slope, int genesis);
    void append_marker_at_i(AppendMarkerData &md, int idx, int last_id);
    void append_markers(AMD_vec &md);
    void remove_marker(int i);
    void remove_marker_data(int is, int ie);
    void remove_markers(const Param& param, const Variables& var, int_vec& markers, int_vec2D& markers_in_elem);
    void resize(const int);
    template <class T>
    void write_chkpt_file(T &bin) const;
    template <class T>
    void read_chkpt_file(Variables &var, T &bin_save, T &bin_chkpt);
    template <class T>
    void write_save_file(const Variables &var, T &bin) const;
    array_t* calculate_marker_coord(const Variables &var) const;

    inline std::string get_name() const { return _name; }

    inline int get_nmarkers() const { return _nmarkers; }
    inline void set_nmarkers(int n) { _nmarkers = n; }

    inline bool if_melt(const int mat) const { return (std::find((*_mattype).begin(), (*_mattype).end(), mat) != (*_mattype).end()); }

    inline int get_id(int m) const { return (*_id)[m]; }
    inline void set_id(const int m, const int i) { (*_id)[m] = i; }

    inline int get_elem(int m) const { return (*_elem)[m]; }
    inline void set_elem(const int m, const int e) { (*_elem)[m] = e; }

    inline int get_mattype(int m) const { return (*_mattype)[m]; }
    inline void set_mattype(const int m, const int mt) { (*_mattype)[m] = mt; }

    inline double get_time(int m) const { return (*_time)[m]; }
    inline void set_time(const int m, const double ti) { (*_time)[m] = ti; }

    inline double get_z(int m) const { return (*_z)[m]; }
    inline void set_z(const int m, const double z) { (*_z)[m] = z; }

    inline double get_distance(int m) const { return (*_distance)[m]; }
    inline void set_distance(const int m, const double d) { (*_distance)[m] = d; }

    inline double get_slope(int m) const { return (*_slope)[m]; }
    inline void set_slope(const int m, const double s) { (*_slope)[m] = s; }

    inline const double *get_eta(int m) const { return (*_eta)[m]; }
    inline void set_eta( const int i, const double r[NDIMS] );

    inline double get_genesis(int m) const { return (*_genesis)[m]; }

    inline double get_tmp(int m) const { return (*_tmp)[m]; }
    inline void set_tmp(const int m, const double tmp) { (*_tmp)[m] = tmp; }

    #pragma acc routine seq
    void get_ZPT(const Param& param, const Variables& var, int m, double &Z, double &P, double &T) const;

private:
    const std::string _name;

    // Didn't create a data type for an individual marker to follow the "structure of arrays" concept.
    // Number of markers (may change during remeshing)

    int _nmarkers;
    int _reserved_space;
    int _last_id;

    // Barycentric (local) coordinate within the reference element
    shapefn *_eta;
    // Containing element
    int_vec *_elem;
    // Material type
    int_vec *_mattype;
    // Unique id
    int_vec *_id;
    // Cearte time
    double_vec *_time;
    // Sedimentation: Depth
    double_vec *_z;
    // Sedimentation: Distance to coastline
    double_vec *_distance;
    // Sedimentation: Slope of surface
    double_vec *_slope;
    // origin of markers:
    // 0: IC,
    // 1: remeshing replenishment,
    // 2: depositon,
    // 3: erosional replenishment with nn,
    // 4: erosional replenishment with interpolation
    int_vec *_genesis;
    // temporary storage
    double_vec *_tmp;

    void random_markers( const Param&, Variables&, int = 0 );
    void regularly_spaced_markers( const Param&, Variables&, int = 0 );
    void allocate_markerdata( const int );

    int initial_mattype( const Param&, const Variables&,
                         int elem, const double eta[NODES_PER_ELEM],
                         const double *x=NULL );
    int layered_initial_mattype( const Param& param, const Variables &var,
                                 int elem, const double eta[NODES_PER_ELEM],
                                 const double *x);
    int custom_initial_mattype( const Param& param, const Variables &var,
                                int elem, const double eta[NODES_PER_ELEM],
                                const double *x );

};

void remap_markers(const Param&, Variables &, 
                   const array_t &, const conn_t &);
void advect_hydrous_markers(const Param &, const Variables &, double,
                            MarkerSet &, Array2D<int,1> &);

#endif
