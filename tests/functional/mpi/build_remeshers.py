#!/usr/bin/env python3
"""Build pinned standalone ParMmg compatibility tools without installing globally."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


REVISIONS = {'2d': ('https://github.com/MmgTools/parMmg2D.git', 'd7e568f79d126bafc61a94fdb81adc1349f0c1b5'),
             '3d': ('https://github.com/MmgTools/ParMmg.git', 'd2eddc502cae46084ca2947859f333d7a7b23943')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--metis-prefix', type=Path)
    parser.add_argument('--dimensions', choices=['2d', '3d'], nargs='+', default=['2d', '3d'])
    parser.add_argument('--patch-2d-merge', action='store_true',
                        help='Apply the retained experimental metric export/indexing fix; boundary labels remain incompatible')
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1: parser.error('--jobs must be positive')
    root = args.workdir.resolve(); root.mkdir(parents=True, exist_ok=True)
    manifest = dict(completed=False, libraries={})
    def call(command, log):
        with log.open('a') as output:
            subprocess.run(list(map(str, command)), stdout=output, stderr=subprocess.STDOUT, check=True)
    for dim in args.dimensions:
        url, revision = REVISIONS[dim]
        source, build, log = root/f'source-{dim}', root/f'build-{dim}', root/f'build-{dim}.log'
        if not source.exists():
            call(['git', 'clone', '--no-checkout', url, source], log)
            call(['git', '-C', source, 'checkout', '--detach', revision], log)
        actual = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
        if actual != revision: raise RuntimeError(f'Unexpected revision in {source}')
        patch = Path(__file__).with_name('parmmg2d-merge.patch').resolve()
        if dim == '2d' and args.patch_2d_merge:
            result = subprocess.run(['git', '-C', str(source), 'apply', '--reverse', '--check', str(patch)], capture_output=True)
            if result.returncode:
                call(['git', '-C', source, 'apply', '--check', patch], log)
                call(['git', '-C', source, 'apply', patch], log)
        elif subprocess.check_output(['git', '-C', str(source), 'diff', '--name-only'], text=True).strip():
            raise RuntimeError(f'Unrequested source changes in {source}; use a fresh workdir')
        command = ['cmake', '-S', source, '-B', build, '-DCMAKE_BUILD_TYPE=Release',
                   '-DDOWNLOAD_METIS=OFF', '-DUSE_VTK=OFF', f'-DCMAKE_INSTALL_PREFIX={root}/install-{dim}']
        if args.metis_prefix:
            prefix = args.metis_prefix.resolve()
            candidates = [prefix/'lib', prefix/'lib64']+list((prefix/'lib').glob('*-linux-gnu'))
            library = next((p for p in candidates if list(p.glob('libmetis.*'))), None)
            if library is None: raise FileNotFoundError('No METIS library below --metis-prefix')
            command += [f'-DMETIS_DIR={prefix}', f'-DMETIS_LIBDIR={library}']
        call(command, log)
        call(['cmake', '--build', build, f'-j{args.jobs}'], log)
        executable = build/'bin'/('parmmg2d_O3' if dim == '2d' else 'parmmg_O3')
        diff = subprocess.check_output(['git', '-C', str(source), 'diff'], text=True)
        (root/f'changes-{dim}.patch').write_text(diff)
        manifest['libraries'][dim] = dict(repository=url, revision=actual, executable=str(executable),
            executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
            patch_sha256=hashlib.sha256(diff.encode()).hexdigest(), configuration=list(map(str, command)))
        (root/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print(dim, executable, flush=True)
    manifest['completed'] = True
    (root/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__': main()
