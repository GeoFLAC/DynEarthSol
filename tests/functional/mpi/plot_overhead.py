#!/usr/bin/env python3
"""Plot controlled MPI overhead timings, payloads and per-rank profiles."""
import argparse
import json
import os
from pathlib import Path
import statistics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.results.read_text())
    if not data.get('completed'):
        parser.error('The benchmark must finish before plotting')
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(args.output / '.matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    cases = list(dict.fromkeys(row['case'] for row in data['measurements']))
    workers = sorted({row['workers'] for row in data['measurements']})
    modes = [('before', 'Previous MPI'), ('blocking', 'Updated, blocking'),
             ('overlap', 'Updated, overlap')]
    records = {(r['case'], r['mode'], r['workers']): r for r in data['measurements']}
    summary = {'timing': [], 'profiles': []}
    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 3.1 * len(cases)),
                             squeeze=False, constrained_layout=True)
    for case, row in zip(cases, axes):
        for mode, label in modes:
            samples = [records[case, mode, n] for n in workers]
            medians = np.array([r['median_seconds'] for r in samples])
            low = np.array([min(s['simulation_s'] for s in r['samples']) for r in samples])
            high = np.array([max(s['simulation_s'] for s in r['samples']) for r in samples])
            row[0].errorbar(workers, medians, yerr=[medians-low, high-medians],
                            marker='o', capsize=4, label=label)
            root_mb = [statistics.median(s['communication']['root_collection_bytes']
                                        for s in r['samples']) / 1e6 for r in samples]
            row[1].plot(workers, root_mb, 'o-', label=label)
            for n, record in zip(workers, samples):
                before = records[case, 'before', n]['median_seconds']
                summary['timing'].append(dict(case=case, mode=mode, ranks=n,
                    median_seconds=record['median_seconds'],
                    speedup=before/record['median_seconds'],
                    min_seconds=min(s['simulation_s'] for s in record['samples']),
                    max_seconds=max(s['simulation_s'] for s in record['samples'])))
        row[0].set(title=case, ylabel='Simulation time (s)')
        row[1].set(title='Logical root collection payload', ylabel='MB sent to/from root')
        for ax in row:
            ax.set(xlabel='MPI ranks (one CPU thread each)', xticks=workers)
            ax.grid(alpha=.25)
            ax.legend(fontsize=8)
    fig.suptitle('Same CPU IDs; three interleaved samples, ranges shown; physical outputs exact')
    for extension in ('png', 'pdf'):
        fig.savefig(args.output / f'mpi-overhead.{extension}', dpi=180)
    plt.close(fig)

    operations = ['Stress', 'Temperature', 'Force', 'Mass']
    profile_cases = [case for case in cases if case != 'small-fixed']
    fig, axes = plt.subplots(len(profile_cases), len(workers), squeeze=False,
                             figsize=(5.5*len(workers), 3.6*len(profile_cases)),
                             constrained_layout=True)
    for case, row in zip(profile_cases, axes):
        for n, ax in zip(workers, row):
            samples = records[case, 'overlap', n]['samples']
            for index, operation in enumerate(operations):
                if not any(p['operation'] == operation for s in samples for p in s['rank_profile']):
                    continue  # Fixed meshes do not recompute mass during timesteps.
                times = []
                comm = []
                for rank in range(n):
                    profiles = [next(p for p in s['rank_profile']
                                     if int(p['rank']) == rank and p['operation'] == operation)
                                for s in samples]
                    times.append(statistics.median(p['outside_communication_s'] for p in profiles))
                    comm.append(statistics.median(p['communication_s'] for p in profiles))
                imbalance = max(times) / statistics.mean(times) if sum(times) else 1
                x = np.arange(n)+(index-1.5)*.19
                ax.bar(x, times, width=.19, label=f'{operation} ({imbalance:.2f} max/mean)')
                summary['profiles'].append(dict(case=case, ranks=n, operation=operation,
                    median_outside_communication_seconds=times,
                    median_communication_seconds=comm, max_over_mean=imbalance))
            ax.set(title=f'{case}: {n} ranks', xlabel='Rank', xticks=range(n),
                   ylabel='Seconds outside instrumented communication')
            ax.grid(axis='y', alpha=.25)
            ax.legend(fontsize=8)
    fig.suptitle('Per-rank medians, overlap enabled: includes packing, scheduling and cache effects')
    for extension in ('png', 'pdf'):
        fig.savefig(args.output / f'rank-profile.{extension}', dpi=180)
    plt.close(fig)
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(args.output / 'mpi-overhead.png')
    print(args.output / 'rank-profile.png')


if __name__ == '__main__':
    main()
