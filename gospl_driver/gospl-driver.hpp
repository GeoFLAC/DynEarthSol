#ifndef GOSPL_DRIVER_HPP
#define GOSPL_DRIVER_HPP

#ifdef HAS_GOSPL_CPP_INTERFACE

#include <vector>
#include <string>

// Forward declaration for gospl_extensions types
typedef int ModelHandle;

/**
 * Structure for elevation statistics
 */
struct ElevationStats {
    double min_elev;
    double max_elev;
    double mean_elev;
    double rms_change;
    int significant_changes;
    
    ElevationStats() : min_elev(0), max_elev(0), mean_elev(0), rms_change(0), significant_changes(0) {}
};

/**
 * GoSPL Driver Class
 * 
 * This class provides an interface for integrating goSPL (Landscape Evolution Model)
 * with DynEarthSol through the gospl_extensions C++ interface.
 * 
 * Features:
 * - Model initialization and cleanup
 * - Time-stepping control
 * - Velocity field application
 * - Elevation tracking and analysis
 * - Interpolation capabilities
 */
class GoSPLDriver {
public:
    ModelHandle model_handle;
    bool initialized;
    bool python_initialized;
    std::vector<std::vector<double>> elevation_history;
    std::vector<double> time_history;
    
    // GoSPL mesh domain bounds [x_min, x_max, y_min, y_max]
    double mesh_bounds[4];
    bool mesh_bounds_valid;
    
    // Coupling frequency control: couple every N DES steps
    int coupling_frequency;      // Run GoSPL every N steps (default: 1)
    int step_counter;            // Current step count since last coupling
    double accumulated_dt;       // Accumulated time since last coupling (in years)

    // --- Gradual application state ---
    std::vector<double> pending_erosion_rate;    // dh/dt rate to apply per DES step (per node, in m/s)
    std::vector<double> prev_erosion_rate;       // Previous cycle's rate (for extrapolation)
    bool has_prev_rate;                          // Whether prev_erosion_rate is valid (false on first cycle)
    int remaining_steps;                         // Steps remaining to apply pending rate

    // --- Adaptive coupling frequency ---
    int base_coupling_frequency;                 // User-specified baseline N
    int adaptive_coupling_frequency;             // Currently active N (may differ from base)
    double rate_change_tolerance;                // Threshold for adapting N (default: 0.3 = 30% change)
    double rate_change_metric;                   // Last computed rate-change metric (for diagnostics)

    // --- Velocity coupling ---
    bool velocity_coupling;                      // Send DES vertical velocities to GoSPL as upsub (default: false)
    
    /**
     * Constructor
     */
    GoSPLDriver();
    
    /**
     * Destructor - automatically calls cleanup
     */
    ~GoSPLDriver();
    
    /**
     * Initialize Python and gospl_extensions only (without loading a model)
     * This is needed before calling generate_mesh()
     * 
     * @return true on success, false on failure
     */
    bool init_python();
    
    /**
     * Initialize the GoSPL model with a configuration file
     * 
     * @param config_path Path to goSPL configuration YAML file
     * @return true on success, false on failure
     */
    bool initialize(const std::string& config_path);
    
    /**
     * Clean up resources and finalize GoSPL
     */
    void cleanup();
    
    /**
     * Check if the driver is initialized
     * 
     * @return true if initialized, false otherwise
     */
    bool is_initialized() const { return initialized; }
    
    /**
     * Run GoSPL processes for a specific time step
     *
     * @param dt Time step size
     * @param verbose Print progress information (default: false)
     * @param skip_tectonics Skip internal getTectonics() call (default: false).
     *        Set true when using apply_velocity_data() to prevent GoSPL from
     *        overwriting upsub with its own tectonic file.
     * @return Elapsed time on success, -1.0 on error
     */
    double run_processes_for_dt(double dt, bool verbose = false, bool skip_tectonics = false);
    
    /**
     * Run GoSPL processes for a specified number of steps
     * 
     * @param num_steps Number of steps to run
     * @param dt Time step size
     * @param verbose Print progress information (default: false)
     * @return Number of steps completed on success, -1 on error
     */
    int run_processes_for_steps(int num_steps, double dt, bool verbose = false);
    
    /**
     * Run GoSPL processes until target time is reached
     * 
     * @param target_time Target simulation time
     * @param dt Time step size
     * @param verbose Print progress information (default: false)
     * @return Number of steps completed on success, -1 on error
     */
    int run_processes_until_time(double target_time, double dt, bool verbose = false);
    
    /**
     * Apply elevation data to update GoSPL mesh with external (DES) topography
     * 
     * This must be called BEFORE run_processes_for_dt() to ensure GoSPL
     * starts from the current DynEarthSol surface topography.
     * 
     * @param coords Pointer to coordinates array (num_points * 3)
     * @param elevations Pointer to elevation values array (num_points)
     * @param num_points Number of data points
     * @param k Number of nearest neighbors for interpolation (default: 3)
     * @param power Inverse distance weighting power (default: 1.0)
     * @return 0 on success, -1 on error
     */
    int apply_elevation_data(const double* coords, const double* elevations,
                            int num_points, int k = 3, double power = 1.0);
    
    /**
     * Apply vertical velocity data to GoSPL's internal uplift/subsidence field (upsub).
     *
     * Call this BEFORE run_processes_for_dt() and with skip_tectonics=true so that
     * GoSPL uses the DES-computed vertical velocity rather than its own tectonic file.
     *
     * @param coords Pointer to coordinates array (num_points * 3), in metres
     * @param velocities Pointer to velocities array (num_points * 3), in m/yr.
     *        Only the vertical component (index [i*3+2] for 3-D) is used by GoSPL.
     * @param num_points Number of data points
     * @param timer Current simulation time in years (used for bookkeeping)
     * @param k Number of nearest neighbours for interpolation (default: 3)
     * @param power Inverse distance weighting power (default: 1.0)
     * @return 0 on success, -1 on error
     */
    int apply_velocity_data(const double* coords, const double* velocities,
                            int num_points, double timer,
                            int k = 3, double power = 1.0);

    /**
     * Interpolate elevation field to external points
     * 
     * @param coords Pointer to coordinates array (num_points * 3)
     * @param num_points Number of points to interpolate to
     * @param elevations Output elevations array (must be pre-allocated)
     * @param k Number of nearest neighbors for interpolation (default: 3)
     * @param power Inverse distance weighting power (default: 1.0)
     * @return 0 on success, -1 on error
     */
    int interpolate_elevation_to_points(const double* coords, int num_points,
                                       double* elevations, int k = 3, double power = 1.0);
    
    /**
     * Get current simulation time
     * 
     * @return Current time on success, -1.0 on error
     */
    double get_current_time();
    
    /**
     * Get model time step
     * 
     * @return Time step on success, -1.0 on error
     */
    double get_time_step();
    
    /**
     * Calculate elevation statistics from an array of elevation values
     * 
     * @param elevations Vector of elevation values
     * @return ElevationStats structure with computed statistics
     */
    ElevationStats calculate_elevation_stats(const std::vector<double>& elevations);
    
    /**
     * Analyze elevation changes between two time points
     * 
     * @param z_before Elevation values before
     * @param z_after Elevation values after
     * @param step_info Optional step information for logging (default: "")
     * @return ElevationStats structure with change analysis
     */
    ElevationStats analyze_elevation_changes(const std::vector<double>& z_before, 
                                           const std::vector<double>& z_after,
                                           const std::string& step_info = "");
    
    /**
     * Set GoSPL verbose mode
     * 
     * This controls the detailed progress output from GoSPL processes.
     * 
     * @param verbose true to enable verbose output, false to suppress
     * @return 0 on success, -1 on error
     */
    int set_verbose(bool verbose);
    
    /**
     * Generate a GoSPL mesh from surface node coordinates
     * 
     * Creates a regular mesh using Delaunay triangulation that covers
     * the same domain as the DynEarthSol surface nodes. The mesh has
     * approximately sqrt(n_nodes) x sqrt(n_nodes) vertices.
     * 
     * @param x_coords X coordinates of surface nodes
     * @param y_coords Y coordinates of surface nodes
     * @param output_file Path to output NPZ file for GoSPL mesh
     * @param resolution Desired mesh spacing (in same units as coordinates). If <= 0, uses sqrt(n_nodes) approach.
     * @param initial_topo_amplitude Amplitude of random initial topography (in meters). 0 = flat.
     * @param mesh_perturbation Fraction of grid spacing to randomly perturb nodes (0-1). 0 = regular grid, 0.3 = moderate perturbation.
     * @return 0 on success, -1 on error
     */
    int generate_mesh(const std::vector<double>& x_coords,
                      const std::vector<double>& y_coords,
                      const std::string& output_file,
                      double resolution = -1.0,
                      double initial_topo_amplitude = 100.0,
                      double mesh_perturbation = 0.0);
    
    /**
     * Demonstrate elevation interpolation capabilities
     * This creates a test grid and shows interpolation results with different parameters
     */
    void demonstrate_elevation_interpolation();
    
    /**
     * Run a controlled simulation with elevation tracking
     * 
     * @param duration Total simulation duration (default: 5.0)
     * @param dt Time step size (default: 1.0)
     */
    void run_controlled_simulation_with_elevation_tracking(double duration = 5.0, double dt = 1.0);

private:
    /**
     * Create a time-dependent velocity field at specified coordinates
     * 
     * @param t Current time
     * @param coords Input coordinates array (num_points * 3)
     * @param velocities Output velocities array (num_points * 3)
     * @param num_points Number of points
     */
    void create_velocity_field_at_coords(double t, double* coords, double* velocities, int num_points);
    
    /**
     * Print final analysis of the simulation
     * 
     * @param num_steps Number of simulation steps completed
     */
    void print_final_analysis(int num_steps);
};

#endif // HAS_GOSPL_CPP_INTERFACE

#endif // GOSPL_DRIVER_HPP
