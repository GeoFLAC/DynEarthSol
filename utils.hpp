#ifndef DYNEARTHSOL3D_UTILS_HPP
#define DYNEARTHSOL3D_UTILS_HPP

#include "parameters.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <cfloat>
#include <math.h>
#include <iomanip>
#if defined(_WIN32)
#include <windows.h>
#else
#include <time.h>
#endif

// Process exit codes: first digit the category, second the cause, so a bare status
// says where to look -- 1x the user to fix, 2x the environment, 3x-6x ours. Keep in
// step with DEVELOPING.md. Nothing above 99; the shell reserves 126+.
enum ExitCode {
    EXIT_OK                   =  0,

    EXIT_CONFIG               = 10,  // generic config error
    EXIT_CONFIG_VALUE         = 11,  // value out of range / unknown option
    EXIT_CONFIG_DATA          = 12,  // malformed input data file (.poly, .exo)

    EXIT_IO_OPEN              = 20,  // cannot open file
    EXIT_IO_RW                = 21,  // read/write failed (incl. HDF5)
    EXIT_IO_RESTART           = 22,  // restart/checkpoint mismatch

    EXIT_UNSUPPORTED_DIM      = 30,  // not implemented for this NDIMS
    EXIT_UNSUPPORTED_LIB      = 31,  // optional library not compiled in

    EXIT_MESH_TETGEN          = 40,  // Triangle / TetGen
    EXIT_MESH_MMG             = 41,  // MMG
    EXIT_MESH_QUALITY         = 42,  // mesh quality / topology

    EXIT_RUNTIME_NAN          = 50,  // NaN or non-finite state
    EXIT_RUNTIME_LOOKUP       = 51,  // marker/geometry lookup failed
    EXIT_RUNTIME_RESOURCE     = 52,  // resource exhausted

    EXIT_INTERNAL_ASSERT      = 60,  // assertion / invariant violated
    EXIT_INTERNAL_UNREACHABLE = 61   // unreachable branch
};

// Terminate with `code`, naming the category so the number never has to be looked up.
// One-argument form where the site already printed its diagnostic. [[noreturn]] is what
// keeps -Wreturn-type and -Wsometimes-uninitialized alive across the call sites.
//
// Defined in utils.cxx, not inline: this header also carries binary_search_index() with an
// `acc routine seq`, and where device code calls that routine nvc++ walks the header's
// statics and gives die() a device version too, which its std::cerr cannot support
// (NVC++-W-1053). A declaration alone leaves nothing to promote.
[[noreturn]] void die(ExitCode code);
[[noreturn]] void die(ExitCode code, const char* msg);

static void print(std::ostream& os, const double& x)
{
  os << x;
}


static void print(std::ostream& os, const int& x)
{
  os << x;
}


static void print(std::ostream& os, const std::size_t& x)
{
  os << x;
}


template <typename T1, typename T2>
void print(std::ostream& os, const std::pair<T1,T2>& x)
{
    os << x.first << ':' << x.second;
}


template <typename T>
void print(std::ostream& os, const T& A, std::size_t size)
{
  os << "[";
  for (std::size_t i = 0; i != size; ++i) {
    print(os, A[i]);
    if (i+1 != size)
      os << ", ";
  }
  os << "]";
}


template <typename Array>
void print(std::ostream& os, const Array& A)
{
  typename Array::const_iterator i;
  os << "[";
  for (i = A.begin(); i != A.end(); ++i) {
    print(os, *i);
    os << ", ";
  }
  os << "]";
}

#ifdef NPROF

#pragma acc routine seq
static inline double log_lookup(const double_vec& log_table, double x) {
    // Error: log_lookup called with x <= 0;
    if (x <= 0.0)
       return NAN;

    int exponent = 0;
    while (x < LOG_XMIN) {
        x *= 10.0;
        exponent -= 1;
    }
    while (x > LOG_XMAX) {
        x *= 0.1;
        exponent += 1;
    }

    int idx = static_cast<int>((x - LOG_XMIN) / LOG_XDELTA);
    double dx = x - (LOG_XMIN + idx * LOG_XDELTA);
    double slope = (log_table[idx + 1] - log_table[idx]) / LOG_XDELTA;
    return log_table[idx] + slope * dx + exponent * 2.302585092994046; // LN_10
}

#pragma acc routine seq
static inline double tan_lookup(const double_vec& tan_table, const double x) {
    // Error: tan_lookup called with x out of range;
    if (x < TAN_XMIN || x > TAN_XMAX) {
       return NAN;
    }

    int idx = static_cast<int>((x - TAN_XMIN) / TAN_XDELTA);
    double dx = x - (TAN_XMIN + idx * TAN_XDELTA);
    double slope = (tan_table[idx + 1] - tan_table[idx]) / TAN_XDELTA;
    return tan_table[idx] + slope * dx;
}

#pragma acc routine seq
static inline double sin_lookup(const double_vec& sin_table, const double x) {
    // Error: sin_lookup called with x out of range;
    if (x < SIN_XMIN || x > SIN_XMAX) {
       return NAN;
    }

    int idx = static_cast<int>((x - SIN_XMIN) / SIN_XDELTA);
    double dx = x - (SIN_XMIN + idx * SIN_XDELTA);
    double slope = (sin_table[idx + 1] - sin_table[idx]) / SIN_XDELTA;
    return sin_table[idx] + slope * dx;
}

#endif

#pragma acc routine seq
static inline double tan_safe(const double_vec& tan_table, const double x) {
#ifdef NPROF
    return tan_lookup(tan_table, x);
#else
    return std::tan(x);
#endif
}

#pragma acc routine seq
static inline double log_safe(const double_vec& log_table, const double x) {
#ifdef NPROF
    return log_lookup(log_table, x);
#else
    return std::log(x);
#endif
}

#pragma acc routine seq
static inline double sin_safe(const double_vec& sin_table, const double x) {
#ifdef NPROF
    return sin_lookup(sin_table, x);
#else
    return std::sin(x);
#endif
}

#pragma acc routine seq
static double pow_safe(const double_vec& log_table, const double x, const double y) {
#ifdef NPROF
    // Error: pow_safe called with x < 0;
    if (x < 0.) {
       return NAN;
    } else if (x == 0.) {
        return 0.0; // Avoid log(0) which is undefined
    } else
        return exp(y * log_lookup(log_table, x));
#else
    return std::pow(x, y);
#endif
}

#pragma acc routine seq
static double pow_1_5(const double x) {
    return x * sqrt(x);
}

#pragma acc routine seq
static double pow_2(const double x) {
    return x * x;
}

#pragma acc routine seq
template <typename T>
static double trace(T s)
{
#ifdef THREED
    return s[0] + s[1] + s[2];
#else
    return s[0] + s[1];
#endif
}

template <typename T>
static double second_invariant2(T t)
{
#ifdef THREED
    double a = (t[0] + t[1] + t[2]) / 3;
    return ( 0.5 * ((t[0]-a)*(t[0]-a) + (t[1]-a)*(t[1]-a) + (t[2]-a)*(t[2]-a))
             + t[3]*t[3] + t[4]*t[4] + t[5]*t[5] );
#else
    return 0.25*(t[0]-t[1])*(t[0]-t[1]) + t[2]*t[2];
#endif
}

template <typename T>
static double second_invariant(T t)
{
    /* second invariant of the deviatoric part of tensor t
     * defined as: td = deviatoric(t); sqrt( td(i,j) * td(i,j) / 2)
     */
    return std::sqrt(second_invariant2(t));
}


#pragma acc routine seq
static inline int binary_search_index(const int* arr, int size, int target) {
    int low = 0, high = size - 1;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (arr[mid] == target) return mid;
        if (arr[mid] < target) low = mid + 1;
        else high = mid - 1;
    }
    return -1;
}


static int findNearestNeighbourIndex( double x_new, const double_vec& x )
{
    /* find nearest neighbour index for interpolation
     * x vector only can be ascending
     */
    double dist = DBL_MAX;
    int idx = -1;
    for (size_t i = 0; i < x.size(); ++i ) {
        double newDist = x_new - x[i];
        if ( newDist >= 0 && newDist <= dist ) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}


static double interp1(const double_vec& x, const double_vec& y, double x_new)
{
    int idx = findNearestNeighbourIndex( x_new, x);
    double slope = 0;

    if (idx < 0)
        idx = 0;
    else if ( idx < static_cast<int>(x.size()-1) )
        slope = (y[idx+1] - y[idx]) / (x[idx+1] - x[idx]);

    return slope * (x_new-x[idx]) + y[idx];
}

static int64_t get_nanoseconds() {
    #pragma acc wait

    #if defined(_WIN32)
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return (int64_t)((double)counter.QuadPart / frequency.QuadPart * 1e9);
    #else
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1e9 + ts.tv_nsec;
    #endif
}

static void print_time_ns(const int64_t duration) {
    int hours = duration / (int64_t)3600000000000;
    int minutes = (duration % (int64_t)3600000000000) / (int64_t)60000000000;
    double seconds = (duration % (int64_t)60000000000) / 1e9;
    std::cout << std::setw(3) << std::setfill('0') << hours << ":"
    << std::setw(2) << std::setfill('0') << minutes << ":"
    << std::setw(9) << std::fixed << std::setprecision(6) << std::setfill('0') << seconds;
}


#pragma acc routine seq
static int out_nan_error(const char* msg, const int idx0, const int idx1 = -1) {
    if (idx1 >= 0)
        printf("Error: %s[%d][%d] becomes NaN\n", msg, idx0, idx1);
    else
        printf("Error: %s[%d] becomes NaN\n", msg, idx0);
    return 1;
}

static void check_nan(const Variables& var, const char* func_name = nullptr) {
#ifdef NPROF
    nvtxRangePush(__FUNCTION__);
#endif

    // Host scalar: the two loops below reduce into it and the host reads it.
    int is_nan = 0;

#ifndef ACC
    #pragma omp parallel default(none) shared(var,is_nan)
#endif
    {
#ifndef ACC
        #pragma omp for reduction(+:is_nan)
#endif
        #pragma acc parallel loop gang vector reduction(+:is_nan)
        for (int e=0; e<var.nelem;e++) {
            if (std::isnan((*var.volume)[e]))
                is_nan += out_nan_error("volume", e);
            
            if (std::isnan((*var.dpressure)[e]))
                is_nan += out_nan_error("dpressure", e);

            if (std::isnan((*var.viscosity)[e]))
               is_nan +=  out_nan_error("viscosity", e);
            
            for (int i=0; i<NODES_PER_ELEM;i++)
                if(std::isnan((*var.connectivity)[e][i]))
                    is_nan += out_nan_error("connectivity", e, i);

            for (int i=0; i<NSTR; i++)
                if (std::isnan((*var.stress)[e][i]))
                    is_nan += out_nan_error("stress", e, i);

        }

#ifndef ACC
        #pragma omp for reduction(+:is_nan)
#endif
        #pragma acc parallel loop gang vector reduction(+:is_nan)
        for (int n=0; n<var.nnode; n++) {
            if (std::isnan((*var.temperature)[n]))
                is_nan += out_nan_error("temperature", n);

            if (std::isnan((*var.tmass)[n]))
                is_nan += out_nan_error("tmass", n);

            for (int i=0; i<NDIMS; i++) {
                if (std::isnan((*var.force)[n][i]))
                    is_nan += out_nan_error("force", n, i);

                if (std::isnan((*var.vel)[n][i]))
                    is_nan += out_nan_error("vel", n, i);

                if (std::isnan((*var.coord)[n][i]))
                    is_nan += out_nan_error("coord", n, i);                
            }
        }
    }

    if (is_nan > 0) {
        if (func_name) {
            std::cerr << "Error: " << is_nan << " NaN values found in the variables in " << func_name << "." << std::endl;
        } else {
            std::cerr << "Error: " << is_nan << " NaN values found in the variables." << std::endl;
        }
        die(EXIT_RUNTIME_NAN);
    }
#ifdef NPROF
    nvtxRangePop();
#endif
}

static std::string format_with_commas(unsigned long value) {
    std::string s = std::to_string(value);

    int insert_position = s.length() - 3;
    while (insert_position > 0) {
        s.insert(insert_position, ",");
        insert_position -= 3;
    }
    return s;
}

#endif // DYNEARTHSOL3D_UTILS_HPP