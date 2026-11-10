#!/usr/bin/env python3
"""Check invalid active hydraulic material properties at the executable boundary."""

import configparser
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    exe = Path(sys.argv[1]).resolve()
    root = Path(__file__).resolve().parents[3]
    cases = [('porosity', -0.1), ('porosity', 1.1), ('biot_coeff', 0.2),
             ('biot_coeff', 1.1), ('bulk_modulus', 0), ('shear_modulus', -1),
             ('bulk_modulus_s', 0), ('fluid_bulk_modulus', 0),
             ('fluid_visc', 0), ('fluid_rho0', -1), ('hydraulic_perm', -1)]
    with tempfile.TemporaryDirectory(prefix='des-hydraulic-input-') as tmp:
        for key, value in cases:
            cfg = configparser.ConfigParser()
            cfg.optionxform = str
            cfg.read(root / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
            cfg['mat'][key] = f'[{value}]'
            path = Path(tmp) / 'invalid.cfg'
            with path.open('w') as stream:
                cfg.write(stream)
            result = subprocess.run([str(exe), str(path)], cwd=tmp,
                                    env=dict(os.environ, OMP_NUM_THREADS='1'),
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, timeout=10)
            assert result.returncode == 11, (key, value, result.returncode, result.stdout)
            assert 'Invalid hydraulic properties for material 0' in result.stdout, result.stdout
        extra = [
            ('ic', {'water_table_z': 'nan'}),
            ('ic', {'water_table_z': 'inf'}),
            ('ic', {'water_table_z': '-inf'}),
            ('mat', {'permeability_tensor': '[1,1,1,2,0,0]'}),
            ('mat', {'permeability_model': '3'}),
            ('mat', {'permeability_exponent': '[0]'}),
            ('mat', {'permeability_exponent': '[-1]'}),
            ('mat', {'permeability_exponent': '[nan]'}),
            ('mat', {'permeability_exponent': '[inf]'}),
            ('mat', {'pore_healing_time': '[-1]'}),
            ('mat', {'pore_healing_time': '[nan]'}),
            ('mat', {'pore_healing_time': '[inf]'}),
            ('mat', {'pore_healing_temperature': '[0]'}),
            ('mat', {'pore_healing_temperature': '[nan]'}),
            ('mat', {'pore_healing_energy': '[-1]'}),
            ('mat', {'pore_healing_energy': '[inf]'}),
            ('mat', {'pore_healing_time': '[1]', 'porosity_model': '0'}),
            ('mat', {'pore_healing_time': '[1]', 'porosity_model': '1'}),
            ('mat', {'matrix_exchange': '[1]', 'matrix_storage': '[0]'}),
            ('mat', {'matrix_exchange': '[-1]'}),
            ('bc', {'hbc_z1': '4'}),
            ('bc', {'num_hbc_schedule': '2', 'hbc_schedule_time_in_yr': '[1,0]',
                    'hbc_schedule_ratio': '[1,1]'}),
            ('control', {'has_thermal_pressurization': 'true', 'has_thermal_diffusion': 'false'}),
            ('mat', {'retention_model': '[4]'}),
            ('mat', {'retention_model': '[1]', 'retention_alpha': '[0]'}),
            ('mat', {'retention_model': '[1]', 'retention_n': '[1]'}),
            ('mat', {'retention_model': '[2]', 'retention_n': '[0]'}),
            ('mat', {'retention_model': '[3]', 'residual_saturation': '[1]'}),
            ('mat', {'retention_model': '[1]', 'porosity': '[0]'}),
            # Valid laws in combinations whose content balance is not yet accounted.
            ('mat', {'retention_model': '[1]', 'retention_n': '[1.6]', 'porosity_model': '1'}),
            ('mat', {'retention_model': '[3]', 'matrix_storage': '[1e-9]', 'matrix_exchange': '[1e-6]'}),
        ]
        for section, updates in extra:
            cfg = configparser.ConfigParser()
            cfg.optionxform = str
            cfg.read(root / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
            cfg[section].update(updates)
            with path.open('w') as stream:
                cfg.write(stream)
            result = subprocess.run([str(exe), str(path)], cwd=tmp,
                                    env=dict(os.environ, OMP_NUM_THREADS='1'),
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=10)
            assert result.returncode == 11, (section, updates, result.returncode, result.stdout)
    print(f'PASS: {len(cases)+len(extra)} invalid hydraulic inputs rejected before initialization')


if __name__ == '__main__':
    main()
