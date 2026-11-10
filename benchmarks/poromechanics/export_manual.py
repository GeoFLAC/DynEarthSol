#!/usr/bin/env python3
"""Export a sealed review into the DES3D benchmark and tutorial chapters.

The destination is a new upload folder, not a live website checkout.
Figures and inputs are static downloads; chapters use native Markdown/math.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
from urllib.parse import unquote, urlsplit

LINK = re.compile(r'(!?)\[([^\]]*)\]\(([^)\n]+)\)')
PRIVATE_PATH = re.compile(r'(?:(?<![\w:/])|(?<=-I)|(?<=-L)|(?<==))/(?:home|mnt|data|Users)/[^\s\"\'<>),\]`]+')
PUBLIC_PROVENANCE = {'code_rev', 'code_dirty', 'build_os', 'os', 'cpu_model',
                     'omp_threads', 'kernel', 'gpu_model', 'gpu_cuda_driver'}
SECTIONS = dict(groundwater='groundwater-flow', unsaturated='unsaturated-flow',
                consolidation='poroelastic-consolidation', injection='injection-induced-stress',
                dynamics='poroelastic-dynamics', mechanics='material-strength',
                tunnels='circular-openings', heat='heat-and-groundwater')
APPLICATIONS = {'regional', 'embankment', 'drawdown', 'slope', 'rainfall'}
ASSET_URL = '/benchmarks/poromechanics'


def chapter_paths(review, output, provenance):
    """Match the site's existing sections without importing analysis modules."""
    docs = output/'docs/benchmarks/poromechanics'
    mapping = {p: docs/'background'/p.relative_to(review/'appendix')
               for p in (review/'appendix').rglob('*.md')}
    for case in provenance['cases']:
        source = review/Path(case['page']).with_suffix('.md')
        if case['key'] in APPLICATIONS:
            target = output/'docs/tutorial/groundwater-slopes'/source.name
        elif case['key'] in ('rsf', 'fault_heat'):
            name = 'rsf-velocity-step.md' if case['key'] == 'rsf' else source.name
            target = output/'docs/benchmarks/faults-friction'/name
        else:
            target = docs/SECTIONS[case['theme']]/source.name
        assert source.is_file(), source
        mapping[source] = target
    mapping[review/'README.md'] = docs/'index.md'
    mapping[review/'REPRODUCE.md'] = docs/'running.md'
    expected = {p for p in review.rglob('*.md')
                if 'data' not in p.relative_to(review).parts and p.name != 'WORK-QUEUE.md'}
    assert set(mapping) == expected, ('unmapped chapters', expected-set(mapping))
    assert len(set(mapping.values())) == len(mapping), 'duplicate chapter destination'
    return mapping


def portable_path(match):
    parts = Path(match[0]).parts
    for anchor in ('results', 'benchmarks', 'tests', 'doc'):
        if anchor in parts:
            return '/'.join(parts[parts.index(anchor):])
    identity = hashlib.sha256(match[0].encode()).hexdigest()[:12]
    return 'external/'+identity+'/'+parts[-1]


def public_metadata(value):
    """Retain scientific values and hashes without workstation identities."""
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in ('host', 'hostname', 'build_host', 'remote_host', 'machine',
                       'builder', 'username', 'user', 'gpu_uuid', 'gpu_visible_devices'):
                continue
            if key == 'runner' and isinstance(item, str) and '@' in item:
                continue
            if key in ('origin', 'code_origin') and isinstance(item, str) and \
                    (item == 'local' or '@' in item or '://' in item or PRIVATE_PATH.search(item)):
                continue
            name = PRIVATE_PATH.sub(portable_path, key)
            assert name not in result, (key, name)
            result[name] = public_metadata(item)
        return result
    if isinstance(value, list):
        return [public_metadata(item) for item in value]
    return PRIVATE_PATH.sub(portable_path, value) if isinstance(value, str) else value


def public_frame(path):
    """Change runtime provenance only; every physical payload stays identical."""
    if path.suffix == '.vtkhdf':
        import h5py
        import numpy as np
        temporary = path.with_name(path.name+'.public')
        with h5py.File(path, 'r') as source, h5py.File(temporary, 'w') as target:
            for key, value in source.attrs.items(): target.attrs[key] = value
            for key in source:
                if key != 'provenance': source.copy(key, target)
            if 'provenance' in source:
                group = target.create_group('provenance')
                for key, value in source['provenance'].attrs.items():
                    if key in PUBLIC_PROVENANCE: group.attrs[key] = value
            def verify(name, item):
                if name == 'provenance' or name.startswith('provenance/'): return
                other = target[name]
                assert set(item.attrs) == set(other.attrs)
                for key in item.attrs: np.testing.assert_array_equal(item.attrs[key], other.attrs[key])
                if isinstance(item, h5py.Dataset):
                    assert item.shape == other.shape and item.dtype == other.dtype
                    np.testing.assert_array_equal(item[...], other[...])
            source.visititems(verify)
        temporary.replace(path)
        return
    data = path.read_bytes()
    assert data.startswith(b'# DynEarthSol'), path
    entries = [line.split(b'\t') for line in data[:4096].split(b'\0', 1)[0].splitlines()[1:]]
    for i, (name, position) in enumerate(entries):
        if name != b'provenance': continue
        start = int(position); end = int(entries[i+1][1]) if i+1 < len(entries) else len(data)
        lines = data[start:end].decode().splitlines()
        text = '\n'.join(line for line in lines if line.partition('=')[0] in PUBLIC_PROVENANCE)+'\n'
        replacement = text.encode()
        assert len(replacement) <= end-start
        changed = data[:start]+replacement.ljust(end-start, b' ')+data[end:]
        assert changed[:start] == data[:start] and changed[end:] == data[end:]
        assert len(changed) == len(data) and changed[:4096] == data[:4096]
        path.write_bytes(changed)
        return


def portable_metadata(output):
    changes = {}
    for path in output.rglob('*'):
        if not path.is_file(): continue
        before = hashlib.sha256(path.read_bytes()).hexdigest()
        if path.suffix == '.json':
            value = public_metadata(json.loads(path.read_text()))
            path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')
        elif path.suffix in ('.md','.html','.txt','.cfg','.csv'):
            path.write_text(PRIVATE_PATH.sub(portable_path,path.read_text()))
        elif path.suffix in ('.py','.sh') and PRIVATE_PATH.search(path.read_text()):
            raise ValueError(f'Portable source contains an absolute workstation path: {path.relative_to(output)}')
        elif re.search(r'\.(?:save|chkpt)\.\d+(?:\.vtkhdf)?$', path.name):
            public_frame(path)
        after = hashlib.sha256(path.read_bytes()).hexdigest()
        if before != after:
            changes[str(path.relative_to(output))] = dict(original_sha256=before, portable_sha256=after)
    (output/'PORTABLE-METADATA.json').write_text(json.dumps(dict(
        policy='Runtime identities and absolute workstation paths are omitted. Checkpoint physical payloads are verified unchanged.',
        files=changes), indent=2)+'\n')


def export(review, output):
    provenance = json.loads((review/'provenance.json').read_text())
    output.mkdir(parents=True, exist_ok=False)
    docs = output/'docs/benchmarks/poromechanics'
    assets = output/'static/benchmarks/poromechanics'
    docs.mkdir(parents=True); assets.mkdir(parents=True)
    mapping = chapter_paths(review, output, provenance)
    for source in review.iterdir():
        if source.name in ('fig', 'data', 'katex'):
            shutil.copytree(source, assets/source.name)
        elif source.name in ('provenance.json', 'VERSION.json', 'MANIFEST.sha256',
                             'review-checks.json', 'style.css'):
            name = 'REVIEW-MANIFEST.sha256' if source.name == 'MANIFEST.sha256' else source.name
            shutil.copy2(source, assets/name)
    for source in (review/'appendix').rglob('*'):
        if source.is_file() and source.suffix not in ('.md', '.py', '.sh') and source.name != 'analytical-solutions.html':
            target = assets/source.relative_to(review)
            target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, target)
    evidence = (review/'evidence.html').read_text().replace('href="index.html"',
        'href="https://geoflac.github.io/des3d/docs/benchmarks/poromechanics"')
    (assets/'evidence.html').write_text(evidence)

    for source, target in mapping.items():
        def link(match):
            bang, label, value = match.groups()
            url = urlsplit(value)
            if url.scheme or url.netloc or not url.path: return match.group(0)
            original = (source.parent/unquote(url.path)).resolve()
            assert original.is_relative_to(review), (source, value)
            if original == review/'WORK-QUEUE.md':
                return ''
            if original.suffix == '.html' and original.with_suffix('.md') in mapping:
                original = original.with_suffix('.md')
            if original in mapping:
                destination = Path(os.path.relpath(mapping[original], target.parent)).as_posix()
            else:
                if original.is_dir(): original /= 'input.cfg'
                relative = original.relative_to(review)
                if relative.as_posix() == 'MANIFEST.sha256': relative = Path('REVIEW-MANIFEST.sha256')
                assert (assets/relative).is_file(), (source, value, relative)
                destination = ASSET_URL+'/'+relative.as_posix()
                if original.suffix == '.html':
                    destination = 'pathname://'+destination
            if url.fragment:
                fragment = url.fragment.replace('_', '-') if original.stem == 'analytical-solutions' else url.fragment
                destination += '#'+fragment
            return f'{bang}[{label}]({destination})'
        text = source.read_text()
        text = text.replace('[Full work queue](WORK-QUEUE.md). ', '')
        # Implementation history belongs in the downloadable provenance, as
        # required by the manual's policy; the chapter keeps physical sources.
        text = text.replace('Integration of the nonelastic work proposed in Eh Tan’s PR #88; thermal fault weakening',
                            'Thermal fault weakening')
        text = re.sub(r"This integrates Eh Tan's\n\[existing PR #88\]\([^)]+\) and its\ncorrected follow-ups into the current poromechanics owners\.",
                      'The work calculation uses the current poromechanics constitutive routines.', text)
        # Docusaurus indexes heading IDs for anchor checks; raw HTML anchors
        # render but are absent from its document metadata.
        text = re.sub(r'<a id="([^"]+)"></a>\s*\n\s*\[([^\]]+)\]\(([^)]+)\)\.',
                      lambda m: f'## {m[2]} {{#{m[1]}}}\n\n[Worked example]({m[3]}).', text)
        text = LINK.sub(link, text)
        if source == review/'README.md':
            text = '---\ntitle: Poromechanics worked examples\nslug: /benchmarks/poromechanics\nsidebar_position: 0\n---\n\n'+text
        elif source.relative_to(review).as_posix() == 'appendix/evidence/README.md':
            text = '---\ntitle: Supporting verification records\n---\n\n'+text
        elif not text.startswith('---\n'):
            title = next((line[2:] for line in text.splitlines() if line.startswith('# ')), target.stem)
            text = '---\ntitle: '+json.dumps(title, ensure_ascii=False)+'\n---\n\n'+text
        target.parent.mkdir(parents=True, exist_ok=True); target.write_text(text)
    (docs/'_category_.json').write_text(json.dumps(dict(label='Poromechanics', position=6,
        link=dict(type='doc', id='index')), indent=2)+'\n')
    (output/'README.txt').write_text('Upload candidate for GeoFLAC/des3d, under the existing benchmark and tutorial sections.\n'
        'Copy docs/ and static/ into a branch of that repository, then run yarn build with Node >=20.\n'
        'The existing site and released documentation versions are not modified by this exporter.\n'
        'This artifact is prepared locally; it is not a website publication.\n'
        'MANIFEST.sha256 and EXPORT.json verify the exported bytes. REVIEW-MANIFEST.sha256 records the original review.\n'
        'PORTABLE-METADATA.json relates changed runtime metadata to the original hashes; physical checkpoint payloads are unchanged.\n')
    portable_metadata(output)
    files = {str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in output.rglob('*') if p.is_file()}
    checksum = output/'MANIFEST.sha256'
    checksum.write_text(''.join(f'{sha}  {name}\n' for name,sha in sorted(files.items())))
    files[checksum.name] = hashlib.sha256(checksum.read_bytes()).hexdigest()
    manifest = dict(review_commit=provenance['review_commit'], problem_pages=len(provenance['cases']),
        chapters={str(p.relative_to(review)):str(t.relative_to(output)) for p,t in mapping.items()},
        exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), files_sha256=files)
    (output/'EXPORT.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('review', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = export(args.review.resolve(), args.output.resolve())
    print(f"Exported {result['problem_pages']} problem pages and {len(result['files_sha256'])} files.")
