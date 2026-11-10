#!/usr/bin/env python3
"""Plot the drained plastic-opening hold and its analytical recovery reference.

Usage: recovery.py CHECK_RECOVERY_OUTPUT FIGURE_DIRECTORY
The analytical hold starts from the accepted loaded state; no parameters fit.
"""
import argparse
import configparser
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import prepare_figures as pf
import diagram_style as ds
import setup_drawings
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from Dynearthsol import Dynearthsol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    data, output = args.data.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = json.loads((data/'summary.json').read_text())
    assert report['passed']
    folder = data/'plastic-extension'
    cfg = configparser.ConfigParser()
    cfg.read(folder/'hold.cfg')
    value = lambda k: json.loads(cfg['mat'][k])[0]
    model = Dynearthsol(str(folder/'hold'))
    time = np.array(model.time)-model.time[0]
    eta, phi, k = (np.array([model.read_field(f, key).mean() for f in model.frames])
                   for key in ('pore dilation', 'porosity', 'permeability'))
    theta = np.exp(-time/value('pore_healing_time'))
    eta_ref = eta[0]*theta
    phi_ref = value('biot_coeff')-(value('biot_coeff')-phi[0])*np.exp(eta[0]*(1-theta))
    exponent = value('permeability_exponent') if cfg.has_option('mat','permeability_exponent') else 3
    k_ref = value('hydraulic_perm')*(phi_ref/value('porosity'))**exponent
    table = np.column_stack((time, eta, eta_ref, phi, phi_ref, k, k_ref))
    np.savetxt(output/'hold.csv', table, delimiter=',', comments='',
        header='hold_time_s,remaining_dilation,reference_dilation,porosity,reference_porosity,permeability_m2,reference_permeability_m2')
    with plt.rc_context(ds.RC):
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.25), constrained_layout=True)
        sample = np.arange(0, len(time), 5)
        for ax, actual, reference, ylabel in (
            (axes[0], eta*1e4, eta_ref*1e4, r'Recoverable dilation, $10^{-4}$'),
            (axes[1], k/value('hydraulic_perm'), k_ref/value('hydraulic_perm'), r'$k/k_0$')):
            ax.plot(time, reference, color=ds.INK, lw=1.4, label='Analytical hold')
            ax.plot(time[sample], actual[sample], 'o', color=ds.WATER, mfc='white', ms=4,
                    label=f'DES, {model.ndims}D PT')
            ax.set(xlabel='Time since loading stopped (s)', ylabel=ylabel, xlim=(0, time[-1]))
            ax.grid(alpha=.16)
            ax.spines[['top', 'right']].set_visible(False)
        axes[0].legend(frameon=False)
        pf.save_figure(fig, output, 'result')
    setup_drawings.pore_recovery(output)
    errors = dict(max_dilation_error=float(abs(eta-eta_ref).max()),
                  max_porosity_error=float(abs(phi-phi_ref).max()),
                  max_relative_permeability_error=float((abs(k-k_ref)/value('hydraulic_perm')).max()))
    assert errors['max_dilation_error'] < 1e-14 and errors['max_porosity_error'] < 3e-14
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    (output/'manifest.json').write_text(json.dumps(dict(data=str(data), errors=errors,
        source_sha256={str(p): digest(p) for p in (Path(__file__), Path(setup_drawings.__file__), Path(ds.__file__))},
        files_sha256={p.name: digest(p) for p in output.iterdir() if p.is_file()}), indent=2)+'\n')


if __name__ == '__main__':
    main()
