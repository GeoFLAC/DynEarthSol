#!/usr/bin/env python3
"""Check non-water hydrostats and exact restart using the production executable.

Usage: check_reference_density.py EXECUTABLE OUTPUT_DIRECTORY
The output directory is new; retain its inputs, logs and physical output.
"""

from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from Dynearthsol import Dynearthsol


from des_case import digest, read_config, restart_check, write_json

def main():
    exe, output = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    summary = dict(executable=str(exe), executable_sha256=digest(exe), cases={})
    for density in (800, 1200):
        cfg = read_config(ROOT / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
        cfg['mat']['fluid_rho0'] = f'[{density}]'
        cfg['control']['gravity'] = '10'
        case = output / f'rho-{density}'
        restart_check(exe, cfg, case, 3, 1, require_stationary=True, timeout=120)
        model = Dynearthsol(str(case / 'fresh'))
        coordinate = model.read_field(0, 'coordinate')
        expected = -density * 10 * coordinate[:, -1]
        errors = []
        for frame in model.frames:
            pressure = model.read_field(frame, 'pore pressure')
            np.testing.assert_allclose(pressure, expected, atol=1e-6, rtol=1e-9)
            errors.append(float(np.max(np.abs(pressure - expected))))
        summary['cases'][str(density)] = dict(ndims=int(model.ndims),
            max_absolute_pressure_error_Pa=max(errors),
            stationary_and_exact_restart=True)
    write_json(output / 'summary.json', summary)
    print('PASS: rho=800/1200 hydrostats, stationary mechanics/pressure, exact restart')


if __name__ == '__main__':
    main()
