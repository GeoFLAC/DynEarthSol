#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/../.." && pwd)

cd "${script_dir}"
PYTHONDONTWRITEBYTECODE=1 python3 "${repo_root}/benchmarks/simple_shear_rsf/check_simple_shear_benchmark.py" \
  --exe "${repo_root}/dynearthsol2d"
