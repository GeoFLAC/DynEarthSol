#!/usr/bin/env python3
"""Compare all saved/checkpoint physical payloads of two matching DES runs.

This strict comparator is for behavior-preserving refactors. It excludes
runtime/provenance records, but permits no numerical tolerance or missing file.
The same binary/HDF5 reader is used by the existing PT restart checks.
"""

import argparse
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tests' / 'functional'))
from pt_initial_checkpoint import fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('current', type=Path)
    parser.add_argument('--model', default='terzaghi')
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    pattern = re.compile(re.escape(args.model) + r'\.(save|chkpt)\.\d{6}(?:\.vtkhdf)?$')

    def filenames(directory):
        return {p.name for p in directory.iterdir() if pattern.fullmatch(p.name)}

    old_files, new_files = filenames(args.reference), filenames(args.current)
    if not old_files or old_files != new_files:
        parser.error(f'Output file sets differ: missing={sorted(old_files-new_files)}, '
                     f'extra={sorted(new_files-old_files)}; reference files={len(old_files)}')
    ignored = {'walltime_sec', 'provenance'}

    def physical(key):
        if key.startswith('@') and key.endswith('.layout'):
            key = key[1:-len('.layout')]
        return key.split('/')[-1] not in ignored

    differences, payloads = {}, 0
    for name in sorted(old_files):
        # The shared reader takes a format-neutral base path and appends the
        # HDF5 suffix itself when no binary file exists.
        base = name[:-len('.vtkhdf')] if name.endswith('.vtkhdf') else name
        a, b = fields(args.reference / base), fields(args.current / base)
        keys = {key for key in a.keys() | b.keys() if physical(key)}
        payloads += len(keys)
        changed = sorted(key for key in keys if key not in a or key not in b or a[key] != b[key])
        if changed:
            differences[name] = changed
    report = dict(reference=str(args.reference.resolve()), current=str(args.current.resolve()),
                  files_compared=len(old_files), physical_payloads_compared=payloads,
                  ignored_fields=sorted(ignored), differences=differences)
    args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return 1 if differences else 0


if __name__ == '__main__':
    sys.exit(main())
