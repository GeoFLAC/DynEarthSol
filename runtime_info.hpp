#ifndef DYNEARTHSOL_RUNTIME_INFO_HPP
#define DYNEARTHSOL_RUNTIME_INFO_HPP

#include "parameters.hpp"

// Selects and initializes the OpenACC device; a non-ACC build does nothing. Call it
// exactly once and before any compute.
void init_offload_device();
void report_host_runtime_status();
// Describes the device this run actually got -- not a property of the build: an
// openacc=1 binary falls back to the host when no NVIDIA device is visible. Queries
// the OpenACC runtime only, never selects, so it can be skipped or repeated freely.
void report_device_runtime_status();
void report_mesh_info(const Variables& var, const char* tag);

#endif
