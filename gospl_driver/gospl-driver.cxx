#include "gospl-driver.hpp"

#ifdef HAS_GOSPL_CPP_INTERFACE

#include "/home/eunseo/opt/gospl_extensions/cpp_interface/gospl_extensions.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// Constructor
GoSPLDriver::GoSPLDriver() : model_handle(-1), initialized(false) {
}

// Destructor
GoSPLDriver::~GoSPLDriver() {
    cleanup();
}

bool GoSPLDriver::initialize(const std::string& config_path) {
    std::cout << "GoSPL Driver: Initializing with " << config_path << std::endl;
    
    // Initialize gospl extensions
    if (initialize_gospl_extensions() != 0) {
        std::cerr << "Failed to initialize gospl extensions" << std::endl;
        return false;
    }
    
    // Create enhanced model
    model_handle = create_enhanced_model(config_path.c_str());
    
    if (model_handle < 0) {
        std::cerr << "Failed to create enhanced model" << std::endl;
        return false;
    }
    
    double current_time = get_current_time();
    double dt = get_time_step();
    
    std::cout << "GoSPL model initialized at t=" << current_time 
              << ", dt=" << dt << std::endl;
    
    initialized = true;
    return true;
}

void GoSPLDriver::cleanup() {
    if (model_handle >= 0) {
        destroy_model(model_handle);
        model_handle = -1;
    }
    
    if (initialized) {
        finalize_gospl_extensions();
        initialized = false;
    }
}

double GoSPLDriver::run_processes_for_dt(double dt, bool verbose) {
    if (!initialized) return -1.0;
    return ::run_processes_for_dt(model_handle, dt, verbose ? 1 : 0);
}

int GoSPLDriver::run_processes_for_steps(int num_steps, double dt, bool verbose) {
    if (!initialized) return -1;
    return ::run_processes_for_steps(model_handle, num_steps, dt, verbose ? 1 : 0);
}

int GoSPLDriver::run_processes_until_time(double target_time, double dt, bool verbose) {
    if (!initialized) return -1;
    return ::run_processes_until_time(model_handle, target_time, dt, verbose ? 1 : 0);
}

int GoSPLDriver::apply_velocity_data(const double* coords, const double* velocities,
                                    int num_points, double timer, int k, double power) {
    if (!initialized) return -1;
    return ::apply_velocity_data(model_handle, coords, velocities, num_points, timer, k, power);
}

int GoSPLDriver::interpolate_elevation_to_points(const double* coords, int num_points,
                                                double* elevations, int k, double power) {
    if (!initialized) return -1;
    return ::interpolate_elevation_to_points(model_handle, coords, num_points, elevations, k, power);
}

double GoSPLDriver::get_current_time() {
    if (!initialized) return -1.0;
    return ::get_current_time(model_handle);
}

double GoSPLDriver::get_time_step() {
    if (!initialized) return -1.0;
    return ::get_time_step(model_handle);
}

ElevationStats GoSPLDriver::calculate_elevation_stats(const std::vector<double>& elevations) {
    ElevationStats stats;
    
    if (elevations.empty()) return stats;
    
    stats.min_elev = *std::min_element(elevations.begin(), elevations.end());
    stats.max_elev = *std::max_element(elevations.begin(), elevations.end());
    
    double sum = 0.0;
    for (double elev : elevations) {
        sum += elev;
    }
    stats.mean_elev = sum / elevations.size();
    
    return stats;
}

ElevationStats GoSPLDriver::analyze_elevation_changes(const std::vector<double>& z_before, 
                                                     const std::vector<double>& z_after,
                                                     const std::string& step_info) {
    ElevationStats stats_before = calculate_elevation_stats(z_before);
    ElevationStats stats_after = calculate_elevation_stats(z_after);
    
    // Calculate changes
    std::vector<double> changes(z_before.size());
    double sum_sq_change = 0.0;
    
    for (size_t i = 0; i < z_before.size(); i++) {
        changes[i] = z_after[i] - z_before[i];
        sum_sq_change += changes[i] * changes[i];
    }
    
    double rms_change = std::sqrt(sum_sq_change / changes.size());
    
    // Calculate standard deviation for threshold
    double mean_change = 0.0;
    for (double change : changes) {
        mean_change += change;
    }
    mean_change /= changes.size();
    
    double variance = 0.0;
    for (double change : changes) {
        variance += (change - mean_change) * (change - mean_change);
    }
    double std_dev = std::sqrt(variance / changes.size());
    double threshold = std_dev * 2.0;
    
    int significant_changes = 0;
    for (double change : changes) {
        if (std::abs(change) > threshold) {
            significant_changes++;
        }
    }
    
    // Print analysis
    std::cout << "  Elevation Analysis" << step_info << ":" << std::endl;
    std::cout << "    Before - Min: " << std::fixed << std::setprecision(6) 
              << stats_before.min_elev << ", Max: " << stats_before.max_elev 
              << ", Mean: " << stats_before.mean_elev << std::endl;
    std::cout << "    After  - Min: " << stats_after.min_elev 
              << ", Max: " << stats_after.max_elev 
              << ", Mean: " << stats_after.mean_elev << std::endl;
    
    auto min_change = *std::min_element(changes.begin(), changes.end());
    auto max_change = *std::max_element(changes.begin(), changes.end());
    
    std::cout << "    Change - Min: " << min_change << ", Max: " << max_change 
              << ", Mean: " << mean_change << std::endl;
    std::cout << "    RMS change: " << rms_change << std::endl;
    
    if (significant_changes > 0) {
        std::cout << "    Points with significant change (>" << threshold 
                  << "): " << significant_changes << "/" << changes.size() << std::endl;
    }
    
    ElevationStats result;
    result.min_elev = stats_after.min_elev;
    result.max_elev = stats_after.max_elev;
    result.mean_elev = stats_after.mean_elev;
    result.rms_change = rms_change;
    result.significant_changes = significant_changes;
    
    return result;
}

void GoSPLDriver::demonstrate_elevation_interpolation() {
    if (!initialized) return;
    
    std::cout << "\nDemonstrating elevation interpolation:" << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    // Create test points across the domain (11x11 grid = 121 points)
    const int grid_size = 11;
    const int num_test_points = grid_size * grid_size;
    std::vector<double> test_coords(num_test_points * 3);
    std::vector<double> test_elevations(num_test_points);
    
    // Generate test grid
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            int idx = i * grid_size + j;
            test_coords[idx * 3 + 0] = i * 10.0 / (grid_size - 1);  // x: 0 to 10
            test_coords[idx * 3 + 1] = j * 10.0 / (grid_size - 1);  // y: 0 to 10
            test_coords[idx * 3 + 2] = 0.0;                          // z: 0
        }
    }
    
    std::cout << "Interpolating elevation at " << num_test_points << " test points" << std::endl;
    
    // Test different interpolation parameters
    std::vector<int> k_values = {1, 3, 5};
    for (int k : k_values) {
        int result = interpolate_elevation_to_points(test_coords.data(), 
                                                   num_test_points, test_elevations.data(), 
                                                   k, 1.0);
        if (result == 0) {
            ElevationStats stats = calculate_elevation_stats(test_elevations);
            std::cout << "  k=" << k << ": Min=" << std::fixed << std::setprecision(6) 
                      << stats.min_elev << ", Max=" << stats.max_elev 
                      << ", Mean=" << stats.mean_elev << std::endl;
        }
    }
    
    std::vector<double> power_values = {0.5, 1.0, 2.0};
    for (double power : power_values) {
        int result = interpolate_elevation_to_points(test_coords.data(), 
                                                   num_test_points, test_elevations.data(), 
                                                   3, power);
        if (result == 0) {
            ElevationStats stats = calculate_elevation_stats(test_elevations);
            std::cout << "  power=" << power << ": Min=" << stats.min_elev 
                      << ", Max=" << stats.max_elev 
                      << ", Mean=" << stats.mean_elev << std::endl;
        }
    }
}

void GoSPLDriver::run_controlled_simulation_with_elevation_tracking(double duration, double dt) {
    if (!initialized) return;
    
    std::cout << "\nRunning controlled simulation with elevation tracking" << std::endl;
    std::cout << "Duration: " << duration << " time units, dt: " << dt << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    double start_time = get_current_time();
    double target_time = start_time + duration;
    int step = 0;
    
    // Create velocity sampling points (8x8 grid = 64 points, avoiding domain edges)
    const int grid_size = 8;
    const int num_points = grid_size * grid_size;
    std::vector<double> coords(num_points * 3);
    std::vector<double> velocities(num_points * 3);
    std::vector<double> elevations(num_points);
    
    // Initialize coordinates
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            int idx = i * grid_size + j;
            coords[idx * 3 + 0] = 1.0 + i * 8.0 / (grid_size - 1);  // x: 1 to 9
            coords[idx * 3 + 1] = 1.0 + j * 8.0 / (grid_size - 1);  // y: 1 to 9
            coords[idx * 3 + 2] = 0.0;                               // z: will be updated
        }
    }
    
    // Get initial elevations
    int elev_result = interpolate_elevation_to_points(coords.data(), 
                                                    num_points, elevations.data(), 5, 1.0);
    if (elev_result == 0) {
        elevation_history.push_back(elevations);
        time_history.push_back(get_current_time());
        
        ElevationStats initial_stats = calculate_elevation_stats(elevations);
        std::cout << "Initial elevation - Min: " << std::fixed << std::setprecision(6)
                  << initial_stats.min_elev << ", Max: " << initial_stats.max_elev
                  << ", Mean: " << initial_stats.mean_elev << std::endl;
    }
    
    // Run simulation with elevation tracking
    double current_time = start_time;
    while (current_time < target_time) {
        std::cout << "\nStep " << (step + 1) << " - Running processes for dt=" << dt << std::endl;
        
        // Store elevations before step
        std::vector<double> z_before = elevations;
        
        // Create time-dependent velocity field
        create_velocity_field_at_coords(current_time, coords.data(), velocities.data(), num_points);
        
        // Apply velocity data
        int vel_result = apply_velocity_data(coords.data(), velocities.data(), num_points, 
                                           current_time + dt, 3, 1.0);
        if (vel_result != 0) {
            std::cerr << "Failed to apply velocity data at step " << step << std::endl;
        }
        
        // Run processes for one time step
        double elapsed = run_processes_for_dt(dt, false);
        if (elapsed < 0) {
            std::cerr << "Failed to run processes at step " << step << std::endl;
            break;
        }
        
        // Get updated elevations
        std::vector<double> current_elevations(num_points);
        int elev_result = interpolate_elevation_to_points(coords.data(), num_points, 
                                                        current_elevations.data(), 5, 1.0);
        if (elev_result == 0) {
            // Update z-coordinates with new elevations
            for (int i = 0; i < num_points; i++) {
                coords[i * 3 + 2] = current_elevations[i];
            }
            
            // Analyze elevation changes
            std::string step_info = " (Step " + std::to_string(step + 1) + ")";
            analyze_elevation_changes(z_before, current_elevations, step_info);
            
            // Store elevation history
            elevation_history.push_back(current_elevations);
            time_history.push_back(get_current_time());
            
            elevations = current_elevations; // Update for next iteration
        }
        
        current_time = get_current_time();
        step++;
    }
    
    // Final analysis
    print_final_analysis(step);
}

void GoSPLDriver::create_velocity_field_at_coords(double t, double* coords, double* velocities, int num_points) {
    const double center_x = 5.0;
    const double center_y = 5.0;
    const double amplitude = 0.1;
    
    for (int i = 0; i < num_points; i++) {
        double x = coords[i * 3 + 0];
        double y = coords[i * 3 + 1];
        
        // Create rotational velocity field that changes with time
        double dx = x - center_x;
        double dy = y - center_y;
        
        // Time-dependent rotation with varying amplitude
        double omega = 0.1 * (1.0 + 0.5 * std::sin(t * 0.1));
        double vx = -dy * omega * amplitude;
        double vy = dx * omega * amplitude;
        double vz = 0.01 * std::sin(x + t * 0.05) * amplitude;
        
        velocities[i * 3 + 0] = vx;
        velocities[i * 3 + 1] = vy;
        velocities[i * 3 + 2] = vz;
    }
}

void GoSPLDriver::print_final_analysis(int num_steps) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "FINAL ELEVATION ANALYSIS" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    if (elevation_history.size() < 2) return;
    
    const auto& initial_elevations = elevation_history[0];
    const auto& final_elevations = elevation_history.back();
    
    double final_time = get_current_time();
    double start_time = time_history[0];
    
    std::cout << "Total simulation time: " << std::fixed << std::setprecision(2) 
              << (final_time - start_time) << " time units" << std::endl;
    std::cout << "Number of steps: " << num_steps << std::endl;
    
    analyze_elevation_changes(initial_elevations, final_elevations, " (Total)");
    
    // Analyze elevation evolution over time
    std::cout << "\nElevation Evolution Summary:" << std::endl;
    for (size_t i = 0; i < elevation_history.size(); i++) {
        ElevationStats stats = calculate_elevation_stats(elevation_history[i]);
        
        if (i == 0) {
            std::cout << "  t=" << std::fixed << std::setprecision(2) << time_history[i] 
                      << ": Mean elevation = " << std::setprecision(6) << stats.mean_elev 
                      << " (initial)" << std::endl;
        } else {
            ElevationStats prev_stats = calculate_elevation_stats(elevation_history[i-1]);
            double change = stats.mean_elev - prev_stats.mean_elev;
            std::cout << "  t=" << std::fixed << std::setprecision(2) << time_history[i] 
                      << ": Mean elevation = " << std::setprecision(6) << stats.mean_elev 
                      << " (Δ=" << std::showpos << change << std::noshowpos << ")" << std::endl;
        }
    }
    
    std::cout << "\n✓ Simulation completed! Ran for " 
              << std::fixed << std::setprecision(2) << (final_time - start_time) 
              << " time units in " << num_steps << " steps" << std::endl;
}

#endif // HAS_GOSPL_CPP_INTERFACE
