#!/usr/bin/env python3
"""Run a parent and the current DES build on the same finite-grain Biot column.

Output includes input files, raw DES saves, logs, hashes, CSV histories and a
JSON summary. --analyze-only reuses the recorded inputs and provenance.
"""

import argparse
import json
import math
from pathlib import Path

import consolidation as pf


import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/functional/hydromechanics'))
from des_case import digest, run_des, write_config, write_json

def run(output, binaries, biot, timeout, fixed_dt):
    provenance = {}
    for label, bindir in binaries.items():
        base = output / label
        base.mkdir()
        provenance[label] = {}
        for dim in (2, 3):
            case = base / f'terzaghi-{dim}d'
            case.mkdir()
            cfg = pf.make_config(dim)
            bulk = json.loads(cfg['mat']['bulk_modulus'])[0]
            cfg['mat']['biot_coeff'] = f'[{biot}]'
            if biot < 1:
                cfg['mat']['bulk_modulus_s'] = f'[{bulk / (1 - biot):.17g}]'
            if fixed_dt is not None:
                cfg['control']['fixed_dt'] = f'{fixed_dt:.17g}'
            write_config(case / 'input.cfg', cfg)
            exe = bindir / f'dynearthsol{dim}d'
            provenance[label][str(dim)] = dict(executable=str(exe), executable_sha256=digest(exe),
                                               input_sha256=digest(case / 'input.cfg'))
            # Persist before execution so later plotting never attributes an old
            # run to a binary that was rebuilt in the meantime.
            write_json(output / 'provenance.json', provenance)
            print(f'Running {label} {dim}D, alpha={biot}', flush=True)
            run_des(exe, case, 'input', None, timeout)
            if digest(exe) != provenance[label][str(dim)]['executable_sha256']:
                raise RuntimeError('Executable changed during the run; rebuild before benchmarking')


def analyze(output, biot):
    summary = {label: {str(dim): pf.load_case(output/label, dim)['metrics'] for dim in (2, 3)}
               for label in ('before', 'after')}
    summary['biot'] = biot
    summary['provenance'] = json.loads((output/'provenance.json').read_text())
    write_json(output/'summary.json', summary)
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--parent-bindir', type=Path)
    parser.add_argument('--current-bindir', type=Path, default=pf.ROOT)
    parser.add_argument('--biot', type=float, default=0.6)
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--fixed-dt-seconds', type=float,
                        help='Use the same explicit physical timestep in both versions')
    parser.add_argument('--analyze-only', dest='analyze_only', action='store_true', help='Check saved runs')
    args = parser.parse_args()
    output = args.output.resolve()
    if args.analyze_only:
        biot = json.loads((output / 'settings.json').read_text())['biot']
    else:
        if (args.parent_bindir is None or not 0.3 <= args.biot <= 1 or
                not (args.timeout > 0 and math.isfinite(args.timeout))):
            parser.error('Require --parent-bindir, 0.3 <= --biot <= 1 and positive --timeout')
        if args.fixed_dt_seconds is not None and not (args.fixed_dt_seconds > 0 and math.isfinite(args.fixed_dt_seconds)):
            parser.error('--fixed-dt-seconds must be positive and finite')
        output.mkdir(parents=True, exist_ok=False)
        biot = args.biot
        write_json(output / 'settings.json', dict(biot=biot, fixed_dt_seconds=args.fixed_dt_seconds))
        run(output, dict(before=args.parent_bindir.resolve(), after=args.current_bindir.resolve()),
            biot, args.timeout, args.fixed_dt_seconds)
    analyze(output, biot)


if __name__ == '__main__':
    main()
