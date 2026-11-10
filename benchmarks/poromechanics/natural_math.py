#!/usr/bin/env python3
"""Typeset the review's equations with KaTeX instead of equation images.

review.py writes each equation as an SVG image whose alt text is its LaTeX.
This pass renders that LaTeX once, at build time, with KaTeX (server-side,
no browser JavaScript), and replaces the image with the typeset HTML. The
KaTeX stylesheet and fonts are copied next to the pages, so the review stays
readable offline. The Markdown pages already carry $$...$$ for Docusaurus.

    python3 natural_math.py REVIEW --katex NODE_MODULES/katex

The shared build_review.py runs this pass after writing all problem pages.
"""
import argparse
import html
import json
from pathlib import Path
import re
import shutil
import subprocess

EQUATION = re.compile(r'<div class="equation"><img src="[^"]+\.svg" alt="([^"]*)">')
HEAD = '<link rel="stylesheet" href="{prefix}katex/katex.min.css">'
# Body text in the bundled KaTeX serif, so prose and mathematics match as on a
# printed page; headings, navigation and captions stay in the sans face.
STYLE = ('<style>body{font-family:KaTeX_Main,"Latin Modern Roman",Georgia,serif;font-size:17.5px}'
         'h1,h2,h3,nav,.kicker,figcaption,small,th,footer,summary{font-family:system-ui,sans-serif}'
         '.equation{display:flex;align-items:center;gap:16px;margin:8px 0;padding:.35em 0;'
         'overflow-x:auto;overflow-y:hidden}.equation .katex-display{margin:0;text-align:left}'
         '.equation .katex{font-size:1.1em}.katex{font-size:1.05em}</style>')
INLINE = re.compile(r'(?<![\\$])\$([^$<>\n]+?)\$(?!\$)')
RENDER = r"""
const katex = require(process.argv[1]);
let input = ''; process.stdin.on('data', d => input += d);
process.stdin.on('end', () => {
  const out = JSON.parse(input).map(([tex, display]) => katex.renderToString(tex,
      {displayMode: display, throwOnError: true, strict: 'error', output: 'html'}));
  process.stdout.write(JSON.stringify(out));
});
"""


def render(katex, items):
    """items: [(latex, display)], rendered in one node process."""
    result = subprocess.run(['node', '-e', RENDER, str(katex)], input=json.dumps(items),
                            capture_output=True, text=True, check=True)
    return json.loads(result.stdout)


def typeset(review, katex):
    sources = {p: p.read_text() for p in review.rglob('*.html')}
    sources = {p: t for p, t in sources.items()
               if 'katex.min.css' not in t and (EQUATION.search(t) or INLINE.search(t))}
    display = {html.unescape(m) for t in sources.values() for m in EQUATION.findall(t)}
    inline = {html.unescape(m) for t in sources.values() for m in INLINE.findall(t)}
    items = sorted((tex, True) for tex in display) + sorted((tex, False) for tex in inline)
    rendered = dict(zip(items, render(katex, items)))
    for page, text in sources.items():
        prefix = '../'*(len(page.relative_to(review).parts)-1)
        text = EQUATION.sub(lambda m: '<div class="equation">'+rendered[(html.unescape(m.group(1)), True)], text)
        text = INLINE.sub(lambda m: rendered[(html.unescape(m.group(1)), False)], text)
        if 'katex.min.css' not in text:
            text = text.replace('</head>', HEAD.format(prefix=prefix)+STYLE+'</head>', 1)
        page.write_text(text)
    target = review/'katex'
    if target.exists(): shutil.rmtree(target)
    target.mkdir()
    shutil.copy2(katex/'dist/katex.min.css', target/'katex.min.css')
    shutil.copy2(katex/'LICENSE', target/'LICENSE')
    shutil.copytree(katex/'dist/fonts', target/'fonts')
    return len(sources), len(items)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('review', type=Path)
    p.add_argument('--katex', type=Path, required=True, help='folder of the katex npm package')
    a = p.parse_args()
    pages, count = typeset(a.review.resolve(), a.katex.resolve())
    print(f'typeset {count} distinct equations on {pages} pages')
