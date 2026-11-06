#!/usr/bin/env python3
"""Plot actual DES_MPI_PARTITION_TRACE CSVs and matching MPI partition log rows."""
import argparse
import csv
import json
import os
from pathlib import Path
import re

import numpy as np


def load(directory):
    log=(directory/'run.log').read_text()
    generations={}
    for row in re.findall(r'\[mpi\]\[partition\] (.*)',log):
        data={k:int(v) for k,v in re.findall(r'(\w+)=(\d+)',row)}
        generations.setdefault(data['generation'],[]).append(data)
    files=sorted(directory.glob('partition.*.csv'))
    if len(files)!=len(generations):
        raise ValueError('Trace/generation count mismatch')
    result=[]
    for path in files:
        generation=int(path.name.split('.')[-2]);counts=sorted(generations[generation],key=lambda r:r['rank'])
        with path.open() as f:rows=list(csv.DictReader(f))
        if rows and 'n3' in rows[0]:
            raise ValueError('This plotter draws 2D triangles; use a 3D viewer for tetrahedral traces')
        nodes=sorted((r for r in rows if r['entity']=='node'),key=lambda r:int(r['id']))
        elements=sorted((r for r in rows if r['entity']=='element'),key=lambda r:int(r['id']))
        for group in (nodes,elements):
            if [int(r['id']) for r in group]!=list(range(len(group))):
                raise ValueError(f'Missing or duplicate global ownership in {path}')
        nproc=len(counts)
        if [r['rank'] for r in counts]!=list(range(nproc)):raise ValueError('Rank list mismatch')
        for kind,group in [('nodes',nodes),('elements',elements)]:
            actual=np.bincount([int(r['owner']) for r in group],minlength=nproc)
            if list(actual)!=[r['owned_'+kind] for r in counts]:raise ValueError('Counts do not match actual maps')
        result.append(dict(generation=generation,step=counts[0]['step'],counts=counts,
                           points=np.array([[float(n['x']),float(n['z'])] for n in nodes])/1000,
                           cells=np.array([[int(e['n0']),int(e['n1']),int(e['n2'])] for e in elements]),
                           owners=np.array([int(e['owner']) for e in elements])))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs',type=Path,nargs='+')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR',str(a.output/'.matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.colors import ListedColormap,BoundaryNorm
    data=[load(path) for path in a.runs]
    fig,axes=plt.subplots(len(data),2,figsize=(13,2.3*len(data)),squeeze=False,constrained_layout=True)
    for axes_row,history in zip(axes,data):
        nproc=len(history[0]['counts']);colors=ListedColormap(plt.get_cmap('tab10').colors[:nproc]);norm=BoundaryNorm(np.arange(nproc+1)-.5,nproc)
        for ax,record,label in zip(axes_row,[history[0],history[-1]],['Initial ownership','Rebalanced ownership']):
            points=record['points'];tri=mtri.Triangulation(points[:,0],points[:,1],record['cells'])
            color=ax.tripcolor(tri,facecolors=record['owners'],cmap=colors,norm=norm,rasterized=True)
            ax.set(aspect='equal',xlabel='Distance (km)',ylabel='Elevation (km)',
                   title=f"{nproc} processes: {label}; step {record['step']}; {len(record['cells'])} elements")
        fig.colorbar(color,ax=list(axes_row),ticks=range(nproc),label='MPI owner',fraction=.02,pad=.02)
    fig.suptitle('Measured mesh ownership before and after MMG remeshing')
    fig.savefig(a.output/'partition-maps.png',dpi=180,bbox_inches='tight');fig.savefig(a.output/'partition-maps.pdf',bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(len(data),4,figsize=(14,3.2*len(data)),squeeze=False,constrained_layout=True)
    metrics=[]
    for axes_row,history in zip(axes,data):
        first,last=history[0],history[-1];nproc=len(first['counts']);ranks=np.arange(nproc)
        for ax,field,title in zip(axes_row[:3],['owned_elements','owned_nodes','owned_markers'],['Owned elements','Owned nodes','Owned markers']):
            ax.bar(ranks-.18,[r[field] for r in first['counts']],.36,label='Initial')
            ax.bar(ranks+.18,[r[field] for r in last['counts']],.36,label='After remeshing')
            ax.set(title=f'{nproc} processes: {title}',xlabel='MPI rank',xticks=ranks);ax.grid(axis='y',alpha=.3)
        ax=axes_row[3];e=np.array([r['owned_elements'] for r in last['counts']]);g=np.array([r['ghost_elements'] for r in last['counts']])
        ax.bar(ranks,e,label='Owned');ax.bar(ranks,g,bottom=e,label='Ghost');ax.set(title='Local element storage',xlabel='MPI rank',xticks=ranks);ax.legend(fontsize=8)
        axes_row[0].legend(fontsize=8)
        metrics.append(dict(processes=nproc,generations=len(history),remeshes=len(history)-1,
                            final_element_imbalance=float(e.max()/e.mean()),
                            owned_nodes=[r['owned_nodes'] for r in last['counts']],
                            owned_elements=list(map(int,e)),ghost_elements=list(map(int,g))))
    fig.suptitle('Measured per-rank loads and halo storage; counts checked against ownership CSVs')
    fig.savefig(a.output/'partition-loads.png',dpi=160);fig.savefig(a.output/'partition-loads.pdf');plt.close(fig)
    (a.output/'partition-metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics,indent=2))


if __name__=='__main__':main()
