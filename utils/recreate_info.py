#!/usr/bin/env python
'''Recreate <modelname>.info from DynEarthSol .save.NNNNNN.vtkhdf frame files.

Dynearthsol.py and restart fall back to the frame files when .info is missing,
so this script is only needed to materialize the file again (e.g. for tools
that parse .info themselves). It rebuilds it from the per-frame metadata that
output.cxx writes as FieldData in every frame. Dynearthsol.py does not read the
dt/walltime/nseg columns.

Usage: python utils/recreate_info.py <modelname> [-f]
       <modelname> may include a directory prefix, no extension.
'''
from __future__ import print_function
import argparse, os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Dynearthsol import scan_frames

# must match Output::write_info() in output.cxx
ROW_FMT = '%6d\t%10d\t%12.6e\t%12.4e\t%12.6e\t%8d\t%8d\t%8d\n'


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('modelname', help='model output path prefix, no extension')
    p.add_argument('-f', '--force', action='store_true',
                   help='overwrite an existing .info (backed up to .info.bak)')
    args = p.parse_args()

    info_path = args.modelname + '.info'
    if os.path.exists(info_path) and not args.force:
        print('Error: %s already exists; use -f to overwrite (backs it up to .bak)'
              % info_path, file=sys.stderr)
        sys.exit(1)

    try:
        frames = scan_frames(args.modelname)
    except (OSError, KeyError) as e:
        print('Error: %s' % e, file=sys.stderr)
        sys.exit(2)

    if not frames:
        print('Error: no %s.save.* frame files found' % args.modelname,
              file=sys.stderr)
        sys.exit(1)

    rows = [ROW_FMT % (d['frame'], d['steps'], d['time'], d['dt'], d['walltime'],
                       d['nnode'], d['nelem'], d['nseg'])
            for d in frames]

    if os.path.exists(info_path):
        os.replace(info_path, info_path + '.bak')
        print('Backed up existing file to %s.bak' % info_path)

    with open(info_path, 'w') as f:
        f.writelines(rows)
    print('Wrote %s (%d frames)' % (info_path, len(rows)))


if __name__ == '__main__':
    main()
