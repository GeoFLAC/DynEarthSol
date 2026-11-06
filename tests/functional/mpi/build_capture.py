#!/usr/bin/env python3
"""Relink existing serial DES objects with a test-only MMG capture wrapper (Linux).

The production adapter and remeshing policy are unchanged. GNU ld --wrap also
works with static MMG libraries, unlike LD_PRELOAD. Run make first in --build-dir.
"""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--make-options', nargs='*', default=[])
    p.add_argument('--dimensions', nargs='+', type=int, choices=[2, 3], default=[2, 3])
    a = p.parse_args()
    if not sys.platform.startswith('linux'): p.error('The capture relinker currently requires Linux/GNU ld')
    build = a.build_dir.resolve(); out = a.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    manifest = dict(completed=False, executables={})
    for dim in a.dimensions:
        obj = out/f'capture-{dim}d.o'
        command = ['cc', '-DDES_CAPTURE_LINK_WRAP', '-Wall', '-Wextra',
                   '-I'+str(build/'mmg/build/include'), '-c', str(Path(__file__).with_name('capture_mmg.c').resolve()),
                   '-o', str(obj)]
        if dim == 3: command.insert(1, '-DTHREED')
        subprocess.run(command, check=True)
        variables = [f'ndims={dim}', 'usemmg=1', 'hdf5=1', 'usempi=0']+a.make_options
        # Read the real build's resolved objects and libraries; do not duplicate
        # compiler/dependency selection in this test utility.
        rule = 'des_capture_link:; @echo $(CXX) $(M_OBJS) $(OBJS) $(LDFLAGS) -L$(C3X3_DIR) -l$(C3X3_LIBNAME)'
        link = subprocess.check_output(['make', '-s', '--eval', rule, 'des_capture_link']+variables, cwd=build, text=True)
        args = shlex.split(link.strip()); args[1:1] = [str(obj), f'-Wl,--wrap=MMG{dim}D_mmg{dim}dlib']
        # MMG's generic save dispatcher can pull optional VTK dependencies from
        # a static archive, even though the capture writes only Medit files.
        native = build/f'mmg/build/CMakeFiles/mmg{dim}d.dir/link.txt'
        if native.exists():
            shared = [s for s in shlex.split(native.read_text()) if s.startswith('/') and '.so' in s]
            args += shared+[f'-Wl,-rpath,{path}' for path in sorted({str(Path(s).parent) for s in shared})]
        target = out/f'dynearthsol{dim}d.capture'
        args += ['-o', str(target)]
        with (out/f'link-{dim}d.log').open('w') as log:
            subprocess.run(args, cwd=build, check=True, stdout=log, stderr=subprocess.STDOUT)
        manifest['executables'][str(dim)] = dict(path=str(target), command=args)
        print(target, flush=True)
    manifest['completed'] = True
    (out/'capture-build.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__': main()
