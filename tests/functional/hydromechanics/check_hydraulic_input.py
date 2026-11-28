#!/usr/bin/env python3
"""Check invalid active hydraulic material properties at the executable boundary."""

from pathlib import Path
import sys
import tempfile


from des_case import expect_exit, read_config

def main():
    exe = Path(sys.argv[1]).resolve()
    root = Path(__file__).resolve().parents[3]
    cases = [('porosity', -0.1), ('porosity', 1.1), ('biot_coeff', 0.2),
             ('biot_coeff', 1.1), ('bulk_modulus', 0), ('shear_modulus', -1),
             ('bulk_modulus_s', 0), ('fluid_bulk_modulus', 0),
             ('fluid_visc', 0), ('fluid_rho0', -1), ('hydraulic_perm', -1)]
    with tempfile.TemporaryDirectory(prefix='des-hydraulic-input-') as tmp:
        for key, value in cases:
            cfg = read_config(root / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
            cfg['mat'][key] = f'[{value}]'
            expect_exit(exe, tmp, 'invalid', cfg, 11, timeout=10,
                        message='Invalid hydraulic properties for material 0')
        extra = [
            ({'hbc_z1': '4'}, 'Hydraulic boundary types must be 0, 1, 2 or 3'),
            ({'num_hbc_period_z1': '2', 'hbc_period_z1_time_in_yr': '[1,0]',
              'hbc_period_z1_ratio': '[1,1]'}, 'bc.hbc_period_z1_time_in_yr must be non-decreasing'),
        ]
        for updates, message in extra:
            cfg = read_config(root / 'tests/functional/pt_hydrostatic_equilibrium.cfg')
            cfg['bc'].update(updates)
            expect_exit(exe, tmp, 'invalid', cfg, 11, timeout=10, message=message)
    print(f'PASS: {len(cases)+len(extra)} invalid hydraulic inputs rejected before initialization')


if __name__ == '__main__':
    main()
