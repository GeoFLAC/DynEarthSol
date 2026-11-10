#!/usr/bin/env python3
"""Run the upstream PT analytical/restart checks on serial and compact MPI meshes.

All process commands, logs, configs and field comparisons are retained. Serial
and MPI builds must enable HDF5 (and MMG for --cases mmg-remesh). Comparisons use
an explicit tolerance because ADR scalar reductions can change floating-point
summation order; exactness is measured separately and is never assumed.
"""
import argparse
import configparser
import hashlib
import json
import os
from pathlib import Path
import shlex
import sys
import tempfile

import h5py
import numpy as np

from run_mpi import execute

ROOT = Path(__file__).resolve().parents[3]
FUNCTIONAL = ROOT / 'tests/functional'
CASES = ('elastic', 'maxwell', 'hydrostatic', 'remesh', 'mmg-remesh', 'dynamic', 'oblique', 'running-remesh', 'failure')


def compare_outputs(reference, current, rtol, atol):
    files = sorted(p.name for p in reference.glob('*.vtkhdf'))
    assert files and files == sorted(p.name for p in current.glob('*.vtkhdf')), 'Output file sets differ'
    exact = True
    worst = 0.0
    checked = 0
    worst_field = None
    for filename in files:
        with h5py.File(reference / filename) as left, h5py.File(current / filename) as right:
            a, b = {}, {}
            for source, target in ((left, a), (right, b)):
                def visit(name, obj):
                    if isinstance(obj, h5py.Dataset) and name.split('/')[-1] not in ('walltime_sec', 'provenance'):
                        target[name] = obj[()]
                source.visititems(visit)
            assert a.keys() == b.keys(), filename
            pressure_scale = float(np.max(np.abs(a.get('pore pressure', np.array([]))), initial=0))
            if '.chkpt.' in filename:
                save = reference/filename.replace('.chkpt.', '.save.')
                if save.is_file():
                    with h5py.File(save) as saved:
                        if 'pore pressure' in saved:
                            pressure_scale = float(np.max(np.abs(saved['pore pressure'][()]), initial=0))
            for name, value in a.items():
                other = b[name]
                assert value.shape == other.shape and value.dtype == other.dtype, (filename, name)
                equal = value.tobytes() == other.tobytes()
                exact &= equal
                if value.dtype.kind in 'fc':
                    assert np.isfinite(value).all() and np.isfinite(other).all(), (filename, name)
                    scale = float(np.max(np.abs(value), initial=0))
                    # The projected residual is a cancellation of loaded nodal
                    # forces. Compare it in force units; independent upstream
                    # tests still require its RMS to meet the requested tolerance.
                    if name.split('/')[-1] == 'PT residual force' and 'force' in a:
                        scale = max(scale, float(np.max(np.abs(a['force']), initial=0)))
                    if name.split('/')[-1] == 'pending pore pressure increment':
                        scale = max(scale, pressure_scale)
                    error = float(np.max(np.abs(value-other), initial=0))
                    allowed = atol + rtol*scale
                    assert error <= allowed, (filename, name, error, allowed)
                    fraction = error/allowed if allowed else 0
                    if fraction > worst:
                        worst = fraction
                        worst_field = dict(file=filename, field=name, absolute_error=error,
                                           reference_scale=scale, allowed_error=allowed)
                else:
                    assert equal, (filename, name)
                checked += 1
    return dict(bit_exact=bool(exact), datasets=checked, max_tolerance_fraction=worst, worst_field=worst_field, files=len(files))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--serial', type=Path, required=True)
    parser.add_argument('--mpi', type=Path, required=True)
    parser.add_argument('--serial3d', type=Path, required=True)
    parser.add_argument('--mpi3d', type=Path, required=True)
    parser.add_argument('--launcher', default='mpiexec --bind-to none')
    parser.add_argument('--ranks', nargs='+', type=int, default=[3, 4])
    parser.add_argument('--dimensions', nargs='+', type=int, choices=(2, 3), default=[2, 3])
    parser.add_argument('--options', nargs='+', type=int, choices=(0, 1), default=[0, 1])
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--workdir', type=Path, default=Path(tempfile.gettempdir()))
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--maxwell-steps', type=int, default=100)
    parser.add_argument('--rtol', type=float, default=1e-8)
    parser.add_argument('--atol', type=float, default=1e-12)
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if min(args.ranks) < 1 or args.threads < 1 or args.rtol < 0 or args.atol < 0:
        parser.error('Ranks/threads must be positive and tolerances nonnegative')
    args.workdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='pt-mpi-', dir=args.workdir.resolve()))
    print(work, flush=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), PYTHONPATH=str(ROOT))
    launcher = shlex.split(args.launcher)
    report = dict(completed=False, rtol=args.rtol, atol=args.atol, threads=args.threads,
                  partitioner=os.environ.get('DES_MPI_PARTITIONER', 'rcb'), checks=[])
    executables = {2: (args.serial.resolve(), args.mpi.resolve()),
                   3: (args.serial3d.resolve(), args.mpi3d.resolve())}
    report['executables'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                             for pair in executables.values() for p in pair}
    try:
        for dim in args.dimensions:
            for option in args.options:
                for case in args.cases:
                    if case == 'oblique' and dim == 3:
                        continue
                    fixture = 'pt-elastic-load.cfg' if case in ('elastic', 'maxwell', 'failure') else {
                        'hydrostatic': 'pt_hydrostatic_equilibrium.cfg',
                        'remesh': 'pt_remesh_restart.cfg', 'mmg-remesh': 'pt_remesh_restart.cfg',
                        'dynamic': 'pt_dynamic_restart.cfg', 'oblique': 'pt_oblique_restart.cfg',
                        'running-remesh': 'pt_remesh_restart.cfg'}[case]
                    source = ROOT/'examples'/fixture if case in ('elastic', 'maxwell', 'failure') else FUNCTIONAL/fixture
                    cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
                    cfg.optionxform = str
                    cfg.read(source)
                    cfg['control']['PT_option'] = str(option)
                    if 'remesh' in case:
                        # MMG-adapted meshes need more relaxation passes than
                        # the original non-MMG fixture; keep its tolerances.
                        cfg['control']['PT_max_iter'] = '50000'
                    if case == 'running-remesh':
                        cfg['control']['PT_use_running_scale'] = 'yes'
                    if case == 'failure':
                        cfg['control'].update(PT_max_iter='1', PT_use_running_scale='no',
                                              PT_relative_tolerance='1e-30', PT_absolute_tolerance='0')
                    if case == 'mmg-remesh':
                        cfg['mesh']['remeshing_option'] = '11'
                    if cfg.getint('mesh', 'meshing_option', fallback=1) in (90, 91):
                        cfg['mesh']['poly_filename'] = str((source.parent/cfg['mesh']['poly_filename']).resolve())
                    label = f'{case}-{dim}d-option{option}'
                    config = work/(label+'.cfg')
                    with config.open('w') as stream:
                        cfg.write(stream)
                    reference = None
                    for ranks in (None, *args.ranks):
                        tag = 'serial' if ranks is None else f'mpi{ranks}'
                        output = work/f'{label}-{tag}'
                        binary = executables[dim][ranks is not None]
                        if ranks is not None:
                            wrapper = work/f'launch-{dim}d-{ranks}'
                            wrapper.write_text('#!/usr/bin/python3\nimport os,sys\nos.execvp('+repr(launcher[0])+','+
                                               repr(launcher+['-n',str(ranks),str(binary)])+'+sys.argv[1:])\n')
                            wrapper.chmod(0o700)
                            binary = wrapper
                        if case == 'failure':
                            output.mkdir()
                            with (output/'model.cfg').open('w') as stream:
                                cfg.write(stream)
                            log, _ = execute([str(binary),'model.cfg'],output,env,expected_failure=True)
                            assert '[DES exit 54]' in log and 'status=max_iterations' in log
                            assert 'physical-start mechanical history restored' in log
                            if ranks is not None:
                                report['checks'].append(dict(case=case,dimension=dim,option=option,ranks=ranks,
                                                             expected_exit=54,passed=True))
                                print(label,tag,'collective failure passed',flush=True)
                            continue
                        if case in ('elastic', 'maxwell'):
                            script = FUNCTIONAL/('pt_residual_output.py' if case == 'elastic' else 'pt_maxwell_creep.py')
                            command = [sys.executable,str(script),str(binary),str(config),str(output),'--option',str(option)]
                            if case == 'maxwell':
                                command += ['--ndims',str(dim),'--steps',str(args.maxwell_steps)]
                        else:
                            command = [sys.executable,str(FUNCTIONAL/'pt_initial_checkpoint.py'),str(binary),str(config),str(output)]
                            if 'remesh' in case:
                                command += ['--restart-frame','1','--require-remesh']
                            elif case == 'hydrostatic':
                                command += ['--require-stationary']
                            else:
                                command += ['--restart-frame','1']
                        (work/f'{label}-{tag}-command.json').write_text(json.dumps(command,indent=2)+'\n')
                        execute(command,work,env,f'{label}-{tag}-check.log')
                        if ranks is None:
                            reference=output
                        else:
                            result=compare_outputs(reference,output,args.rtol,args.atol)
                            report['checks'].append(dict(case=case,dimension=dim,option=option,ranks=ranks,**result))
                            print(label,tag,'passed', 'exact' if result['bit_exact'] else 'within tolerance',flush=True)
        report['completed']=True
    finally:
        (work/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    if args.plot and 'maxwell' in args.cases:
        plot_creep(work, args.ranks, args.options, args.dimensions)
    print('Passed. Results:',work/'results.json',flush=True)


def plot_creep(work, ranks, options, dimensions=(2,3)):
    os.environ.setdefault('MPLCONFIGDIR', str(work/'.matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0,str(ROOT))
    from Dynearthsol import Dynearthsol
    fig, axes = plt.subplots(1,len(dimensions),figsize=(5*len(dimensions),4),squeeze=False)
    for ax, dim in zip(axes[0],dimensions):
        for option in options:
            for rank in (None,*ranks):
                tag='serial' if rank is None else f'mpi{rank}'
                model=Dynearthsol(str(work/f'maxwell-{dim}d-option{option}-{tag}'/'creep'))
                strain=[float(np.mean(model.read_field(f,'strain')[:,dim-1])) for f in model.frames]
                ax.plot(model.time,strain,label=f'PT {option}, {tag}',linestyle='-' if rank is None else '--')
        analytic = json.loads((work/f'maxwell-{dim}d-option{options[0]}-serial'/'analytic.json').read_text())
        ax.plot(model.time,[0]+[row['strain_zz'] for row in analytic['history']],
                color='black',linestyle=':',linewidth=1.5,label='Analytical reference')
        ax.set(title=f'{dim}D Maxwell creep',xlabel='Time (s)',ylabel='Mean vertical strain')
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(work/'pt-creep.png',dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    main()
