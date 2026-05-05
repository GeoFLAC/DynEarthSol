#ifndef DYNEARTHSOL_RUNTIME_INFO_HPP
#define DYNEARTHSOL_RUNTIME_INFO_HPP

#include "parameters.hpp"

void report_cpu_runtime_status();
void report_openacc_runtime_status();
void report_mesh_info(const Variables& var, const char* tag);

#endif
