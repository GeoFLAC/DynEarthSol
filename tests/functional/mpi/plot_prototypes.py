#!/usr/bin/env python3
"""Plot measured partition maps/halo counts and native remesher compatibility."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from plot_partitions import load
from remesh_prototype import read_mesh


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native', type=Path, required=True)
    p.add_argument('--replay', type=Path, required=True)
    p.add_argument('--partitions', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(a.output/'.matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    def save(fig, stem):
        for suffix in ['png', 'pdf']: fig.savefig(a.output/f'{stem}.{suffix}', dpi=170, bbox_inches='tight')
        plt.close(fig)
    native = json.loads(a.native.read_text()); replay = json.loads(a.replay.read_text())
    rows = [f'{dim}d-{metric}-{constraint}' for dim in (2, 3) for metric in ['scalar', 'tensor'] for constraint in ['free', 'frozen']]
    columns = ['mmg1', 'parmmg2', 'parmmg3', 'parmmg4']
    matrix = np.zeros((len(rows), len(columns)))
    for i, row in enumerate(rows):
        for j, column in enumerate(columns):
            matrix[i, j] = next(r['passed'] for r in native['cases'] if r['case'] == row+'-'+column)
    captures = list(dict.fromkeys(r['capture'] for r in replay['cases']))
    lost = np.zeros((len(captures), 4)); passed = lost.copy()
    for i, capture in enumerate(captures):
        for j, column in enumerate(columns):
            backend, rank = ('mmg', 1) if j == 0 else ('parmmg', j+1)
            r = next(r for r in replay['cases'] if r['capture'] == capture and r['backend'] == backend and r['ranks'] == rank)
            lost[i, j] = r.get('lost_required_elements', 0); passed[i, j] = r['passed']
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    axes[0].imshow(matrix, cmap=ListedColormap(['#d95b59', '#61ad83']), vmin=0, vmax=1, aspect='auto')
    for (i, j), value in np.ndenumerate(matrix): axes[0].text(j, i, 'PASS' if value else 'FAIL', ha='center', va='center', fontsize=9)
    axes[0].set(yticks=range(len(rows)), yticklabels=rows, xticks=range(4), xticklabels=['Serial MMG', 'Par 2', 'Par 3', 'Par 4'], title='Native library checks: two remeshing passes')
    axes[1].imshow(passed, cmap=ListedColormap(['#d95b59', '#61ad83']), vmin=0, vmax=1, aspect='auto')
    for (i, j), value in np.ndenumerate(lost): axes[1].text(j, i, f'{int(value)} lost', ha='center', va='center')
    labels = [f"{next(r['ndims'] for r in replay['cases'] if r['capture']==c)}D capture {i%2+1}" for i, c in enumerate(captures)]
    axes[1].set(yticks=range(len(labels)), yticklabels=labels, xticks=range(4), xticklabels=['Serial MMG', 'Par 2', 'Par 3', 'Par 4'], title='Actual DES inputs: missing required elements')
    fig.suptitle('Compatibility checks, not a speed benchmark; red also includes label/metric/vertex failures')
    save(fig, 'remesher-compatibility')

    # Plot native output labels directly; do not use the optional repair adapter.
    directory = next(a.native.parent/f'2d-scalar-free-parmmg{n}'/'cycle0' for n in (2,3,4)
                     if (a.native.parent/f'2d-scalar-free-parmmg{n}'/'cycle0'/'output.mesh').exists())
    serial = a.native.parent/'2d-scalar-free-mmg1'/'cycle0'
    fig, axes = plt.subplots(3, 1, figsize=(12, 6), constrained_layout=True)
    palette = {1:'#2673b8', 2:'#e8862c', 16:'#535353', 32:'#279746'}
    for ax, path, title in zip(axes, [directory/'input.mesh', serial/'output.mesh', directory/'output.mesh'],
                               ['Original DES snapshot', 'Serial MMG output', 'Native ParMmg2D output']):
        points, cells, faces, refs = read_mesh(path)
        points = points/1000
        ax.triplot(mtri.Triangulation(points[:, 0], points[:, 1], cells), color='lightgray', linewidth=.25)
        ax.add_collection(LineCollection(points[faces], colors=[palette.get(int(r), '#cf18a1') for r in refs], linewidths=2))
        ax.autoscale(); ax.set(aspect='equal', ylabel='Elevation (km)', title=f'{title}: {len(set(refs))} distinct boundary labels')
    axes[-1].set_xlabel('Distance (km)')
    fig.legend(handles=[Patch(color=c, label=k) for k,c in [('Left',palette[1]),('Right',palette[2]),('Bottom',palette[16]),('Top',palette[32]),('Unexpected label','#cf18a1')]], loc='outside lower center', ncol=5)
    save(fig, 'boundary-labels')

    comparisons = json.loads(a.partitions.read_text())
    fig, axes = plt.subplots(2, 2, figsize=(13, 5), constrained_layout=True)
    for row, ranks in enumerate((3, 4)):
        for col, method in enumerate(('rcb', 'metis')):
            record = next(r for r in comparisons if r['ndims']==2 and r['ranks']==ranks and r['method']==method)
            state = load(Path(record['directory']))[-1]
            points = state['points']; tri = mtri.Triangulation(points[:,0], points[:,1], state['cells'])
            axes[row, col].tripcolor(tri, facecolors=state['owners'], cmap=ListedColormap(plt.get_cmap('tab10').colors[:ranks]), vmin=-.5, vmax=ranks-.5)
            axes[row, col].set(aspect='equal', xlabel='Distance (km)', ylabel='Elevation (km)', title=f'{method.upper()}, {ranks} ranks: after MMG')
    fig.suptitle('Actual DES ownership with RCB versus optional METIS; physical outputs identical')
    save(fig, 'partitioner-maps')
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    groups = [(d,n) for d in (2,3) for n in (3,4)]
    for method, offset, color in [('rcb',-.18,'#2673b8'),('metis',.18,'#e8862c')]:
        ghosts, nodes = [], []
        for dim, rank in groups:
            directory = Path(next(r['directory'] for r in comparisons if r['ndims']==dim and r['ranks']==rank and r['method']==method))
            # The map reader is 2D-only; rank log counts work for either dimension.
            import re
            lines = re.findall(r'\[mpi\]\[partition\] (.*)', (directory/'run.log').read_text())
            records = [{k:int(v) for k,v in re.findall(r'(\w+)=(\d+)', line)} for line in lines]
            generation = max(r['generation'] for r in records); records = [r for r in records if r['generation']==generation]
            ghosts.append(sum(r['ghost_elements'] for r in records)); nodes.append(sum(r['ghost_nodes'] for r in records))
        for ax, values, title in zip(axes, [ghosts,nodes], ['Total ghost elements', 'Total ghost nodes']):
            ax.bar(np.arange(4)+offset, values, .36, label=method.upper(), color=color)
            ax.set(xticks=range(4), xticklabels=[f'{d}D / {n} ranks' for d,n in groups], title=title)
            ax.grid(axis='y', alpha=.25); ax.legend()
    fig.suptitle('Measured halo storage after remeshing: METIS is not automatically better')
    save(fig, 'partitioner-halos')


if __name__ == '__main__': main()
