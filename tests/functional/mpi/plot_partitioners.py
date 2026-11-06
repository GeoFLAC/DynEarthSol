#!/usr/bin/env python3
"""Plot paired RCB/METIS timings, partition cost and logical communication."""
import argparse
import json
import os
from pathlib import Path
import re
import statistics


def median(rows, function):
    return statistics.median(function(row) for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.results.read_text())
    if not data.get('completed'):
        parser.error('Benchmark must complete before plotting')
    for row in data['measurements']:
        if 'phase_seconds' not in row:
            log = (Path(row['directory'])/'run.log').read_text()
            row['phase_seconds'] = {name: 3600*int(h)+60*int(m)+float(s)
                for name, h, m, s in re.findall(
                    r'^  (Execute|Initiate|Compute|Remesh|Output)\s*:\s*(\d+):(\d+):([\d.]+)',
                    log, re.MULTILINE)}
            if not all(name in row['phase_seconds'] for name in ['Initiate', 'Compute', 'Output']):
                raise ValueError('Missing solver phase timings')
    args.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(args.output/'.matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    ranks = sorted({r['ranks'] for r in data['measurements']})
    x = np.arange(len(ranks))
    rows_by_key = {}
    for row in data['measurements']:
        rows_by_key.setdefault((row['case'], row['ranks'], row['method']), []).append(row)
    summary = dict(completed=True, repeats=data['repeats'], comparisons=[])

    def bars(ax, name, function, ylabel, title):
        for offset, method in [(-.18, 'rcb'), (.18, 'metis')]:
            rows = [rows_by_key[name, n, method] for n in ranks]
            values = np.array([median(group, function) for group in rows])
            lower = np.array([min(function(r) for r in group) for group in rows])
            upper = np.array([max(function(r) for r in group) for group in rows])
            ax.bar(x+offset, values, width=.36, label=method.upper(),
                   yerr=[values-lower, upper-values], capsize=3)
        ax.set(title=title, ylabel=ylabel, xlabel='MPI ranks (one CPU thread each)',
               xticks=x, xticklabels=ranks)
        ax.grid(axis='y', alpha=.25)
        ax.legend(fontsize=8)

    timing, axes = plt.subplots(len(data['cases']), 3, squeeze=False,
                               figsize=(15, 2.8*len(data['cases'])), constrained_layout=True)
    traffic, traffic_axes = plt.subplots(len(data['cases']), 2, squeeze=False,
                               figsize=(11, 2.8*len(data['cases'])), constrained_layout=True)
    for case, row_axes, comm_axes in zip(data['cases'], axes, traffic_axes):
        name = case['name']
        label = f"{name}: initially {case['initial_elements']:,} elements, {case['steps']} steps"
        bars(row_axes[0], name, lambda r: r['simulation_s'], 'Seconds', label)
        bars(row_axes[1], name, lambda r: r['phase_seconds']['Compute'],
             'Seconds', 'Solver-reported timestep computation')
        bars(row_axes[2], name, lambda r: 1000*r['partition']['seconds'],
             'Milliseconds', 'All partition rebuilds, included in total runtime')
        bars(comm_axes[0], name, lambda r: r['communication']['halo_sent_bytes']/1e6,
             'MB sent across all ranks', label)
        bars(comm_axes[1], name, lambda r: r['communication']['seconds'],
             'Seconds on rank zero', 'Instrumented communication, packing and waiting')
        for n in ranks:
            rcb = rows_by_key[name, n, 'rcb']
            metis = rows_by_key[name, n, 'metis']
            pairs = [(next(r for r in rcb if r['repeat'] == i),
                      next(r for r in metis if r['repeat'] == i)) for i in range(data['repeats'])]
            ratios = [r['simulation_s']/m['simulation_s'] for r, m in pairs]
            item = dict(case=name, ndims=case['ndims'], ranks=n, methods={},
                        metis_faster_pairs=sum(ratio > 1 for ratio in ratios),
                        paired_speedups=ratios, median_paired_speedup=statistics.median(ratios))
            for method, samples in [('rcb', rcb), ('metis', metis)]:
                totals = [r['simulation_s'] for r in samples]
                item['methods'][method] = dict(
                    median_seconds=statistics.median(totals), min_seconds=min(totals),
                    max_seconds=max(totals),
                    median_wall_with_launcher_seconds=median(samples, lambda r: r['wall_with_launcher_s']),
                    median_initialization_seconds=median(samples, lambda r: r['phase_seconds']['Initiate']),
                    median_compute_seconds=median(samples, lambda r: r['phase_seconds']['Compute']),
                    median_output_seconds=median(samples, lambda r: r['phase_seconds']['Output']),
                    median_partition_seconds=median(samples, lambda r: r['partition']['seconds']),
                    median_root_communication_seconds=median(samples, lambda r: r['communication']['seconds']),
                    halo_bytes=median(samples, lambda r: r['communication']['halo_sent_bytes']),
                    root_bytes=median(samples, lambda r: r['communication']['root_collection_bytes']),
                    remeshes=samples[0]['remeshes'],
                    final_owned_elements=[r['owned_elements'] for r in samples[0]['partition']['generations'][-1]],
                    final_ghost_elements=sum(r['ghost_elements'] for r in samples[0]['partition']['generations'][-1]),
                    final_neighbors=sum(r['neighbors'] for r in samples[0]['partition']['generations'][-1]))
            a = item['methods']['rcb']['median_seconds']
            b = item['methods']['metis']['median_seconds']
            item['metis_runtime_reduction_percent'] = 100*(1-b/a)
            item['halo_reduction_percent'] = 100*(1-item['methods']['metis']['halo_bytes']/item['methods']['rcb']['halo_bytes'])
            summary['comparisons'].append(item)
    timing.suptitle(f"RCB versus METIS: {data['repeats']} paired rounds; medians and full sample ranges")
    traffic.suptitle('Logical halo payload and root communication timers; physical outputs exact')
    for name, fig in [('partitioner-timing', timing), ('partitioner-communication', traffic)]:
        for extension in ['png', 'pdf']:
            fig.savefig(args.output/f'{name}.{extension}', dpi=165)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 5.5), constrained_layout=True)
    for offset, n, marker in zip(np.linspace(-.14, .14, len(ranks)), ranks,
                                 ['o', 's', '^', 'D', 'v', 'P', 'X']*len(ranks)):
        records = [next(r for r in summary['comparisons']
                        if r['case'] == case['name'] and r['ranks'] == n) for case in data['cases']]
        values = np.array([r['methods']['rcb']['median_seconds']/r['methods']['metis']['median_seconds']
                           for r in records])
        low = np.array([min(r['paired_speedups']) for r in records])
        high = np.array([max(r['paired_speedups']) for r in records])
        ax.errorbar(values, np.arange(len(records))+offset, xerr=[values-low, high-values],
                    fmt=marker, capsize=4, label=f'{n} ranks')
    ax.axvline(1, color='black', linewidth=1)
    ax.set(yticks=np.arange(len(data['cases'])),
           yticklabels=[case['name'] for case in data['cases']],
           xlabel='RCB time / METIS time: below 1 favors RCB; above 1 favors METIS',
           title='Total simulation time: ratio of medians, bars show the paired sample range')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=.25)
    ax.legend()
    for extension in ['png', 'pdf']:
        fig.savefig(args.output/f'partitioner-speedup.{extension}', dpi=180)
    plt.close(fig)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (args.output/'timings.json').write_text(json.dumps(data, indent=2)+'\n')
    print(args.output/'partitioner-timing.png')
    print(args.output/'partitioner-communication.png')
    print(args.output/'partitioner-speedup.png')


if __name__ == '__main__':
    main()
