#!/usr/bin/env python3
"""Validate local review links, hashes and mandatory case material; make a ZIP."""

import argparse
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
import zipfile


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links = []; self.ids = set()

    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if 'id' in d:
            self.ids.add(d['id'])
        self.links.extend(d[k] for k in ('href', 'src') if k in d)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package(review, archive):
    count = 0
    def check(source, link):
        nonlocal count
        u = urlsplit(link)
        if u.scheme or u.netloc:
            return
        target = (source.parent/unquote(u.path)).resolve() if u.path else source
        # This package must not depend on the raw-results directory.
        assert target.is_relative_to(review), (source, link)
        assert target.is_file(), (source, link)
        if u.fragment:
            if target.suffix == '.html':
                parsed = Links(); parsed.feed(target.read_text()); ids = parsed.ids
            elif target.suffix == '.md':
                text = target.read_text()
                ids = set(re.findall(r'<a id="([^"]+)"', text))
                ids.update(re.sub(r'[^\w -]', '', title.lower()).replace(' ', '-')
                           for title in re.findall(r'^#+ (.+)$', text, flags=re.M))
            else:
                ids = set()
            assert unquote(u.fragment) in ids, (source, link)
        count += 1
    for path in review.rglob('*.html'):
        parsed = Links(); parsed.feed(path.read_text())
        for link in parsed.links:
            check(path, link)
    for path in review.rglob('*.md'):
        for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            check(path, link)
    provenance = json.loads((review/'provenance.json').read_text())
    assert provenance['cases']
    assert len({c['page'] for c in provenance['cases']}) == len(provenance['cases'])
    for name, sha in provenance['source_files'].items():
        assert digest(review.parent/name) == sha, name
    for case in provenance['cases']:
        page = review/case['page']
        assert page.is_file()
        assets = Path(case['schematic']).parent
        for stem in ('setup', 'result'):
            for ext in ('png', 'svg', 'pdf'):
                assert (review/assets/f'{stem}.{ext}').is_file(), stem
        assert case['equation_files']
        assert all((review/f).is_file() for f in case['equation_files'])
        assert case['analytical_solution'] == case['equation_files'][0]
        assert case['derivation'] == case['page']+'#derivation'
    for key in provenance.get('reference_only_cases', []):
        page = review/provenance['reference_pages'][key]
        assert page.is_file()
        assets = review/'fig'/page.stem
        assert (assets/'setup.svg').is_file() and (assets/'e1.svg').is_file()
    assets = {str(p.relative_to(review)): digest(p) for p in sorted(review.rglob('*'))
              if p.is_file() and p.name != 'review-checks.json' and p.suffix != '.log'}
    checks = dict(selected_cases=len(provenance['cases']),
        reference_only_cases=len(provenance.get('reference_only_cases', [])),
        source_files_hashed=len(provenance['source_files']), local_links_checked=count,
        all_local_HTML_and_Markdown_links_valid=True,
        every_case_has_setup_equations_and_derivation=True, parameters_fitted=False,
        artifacts_sha256=assets)
    (review/'review-checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    temp = archive.with_suffix('.zip.tmp')
    files = [p for p in sorted(review.rglob('*')) if p.is_file() and p.suffix != '.log']
    with zipfile.ZipFile(temp, 'w', zipfile.ZIP_DEFLATED) as z:
        for path in files:
            z.write(path, Path('hydromechanics-review')/path.relative_to(review))
    with zipfile.ZipFile(temp) as z:
        assert z.testzip() is None
        assert len(z.namelist()) == len(files)
    os.replace(temp, archive)
    sha = digest(archive)
    archive.with_suffix('.zip.sha256').write_text(sha+'  '+archive.name+'\n')
    return dict(selected_cases=len(provenance['cases']),
                reference_only_cases=len(provenance.get('reference_only_cases', [])), local_links_checked=count,
                source_files_hashed=len(provenance['source_files']), archive=str(archive),
                archive_files=len(files), archive_bytes=archive.stat().st_size, sha256=sha)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('review', type=Path); parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    print(json.dumps(package(args.review.resolve(), args.archive.resolve()), indent=2))
