#!/usr/bin/env python3
"""Plot paired OpenMP/MPI medians and sample ranges from benchmark_scaling.py."""
import argparse
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.results.read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(args.output / '.matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    resolutions = sorted({r['resolution'] for r in data['measurements']}, reverse=True)
    fig, axes = plt.subplots(len(resolutions), 2, figsize=(11, 3.5 * len(resolutions)),
                             squeeze=False, constrained_layout=True)
    for row, resolution in zip(axes, resolutions):
        records = [r for r in data['measurements'] if r['resolution'] == resolution]
        baseline = next(r['median_seconds'] for r in records
                        if r['mode'] == 'openmp' and r['workers'] == 1)
        for mode, label, marker in [('openmp', 'OpenMP threads', 'o'), ('mpi', 'MPI ranks, 1 thread each', 's')]:
            samples = sorted((r for r in records if r['mode'] == mode), key=lambda r: r['workers'])
            workers = np.array([r['workers'] for r in samples])
            median = np.array([r['median_seconds'] for r in samples])
            low = np.array([min(s['simulation_s'] for s in r['samples']) for r in samples])
            high = np.array([max(s['simulation_s'] for s in r['samples']) for r in samples])
            row[0].errorbar(workers, median, yerr=[median-low, high-median],
                            marker=marker, capsize=3, label=label)
            row[1].errorbar(workers, baseline/median,
                            yerr=[baseline/median-baseline/high, baseline/low-baseline/median],
                            marker=marker, capsize=3, label=label)
        for ax in row:
            ax.set(xlabel='CPU budget (threads or processes)', xticks=workers)
            ax.grid(alpha=.25)
            ax.legend(fontsize=8)
        row[0].set(ylabel='Simulation time (s)', title=f'{resolution} m resolution: 1,000 fixed-mesh steps')
        row[1].set(ylabel='Speedup over one-thread parent', title='Equal allowed CPU IDs; verified affinity')
        row[1].axhline(1, color='gray', linewidth=.7)
    fig.suptitle(f"CPU OpenMP versus spatial MPI: medians of {data['repeats']} runs; bars show sample ranges")
    for suffix in ('png', 'pdf'):
        fig.savefig(args.output / f'cpu-scaling.{suffix}', dpi=180)
    plt.close(fig)
    print(args.output / 'cpu-scaling.png')


if __name__ == '__main__':
    main()
