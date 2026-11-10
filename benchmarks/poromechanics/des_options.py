"""DES's input options, read from the executable's own help, and the settings an input file changes.

`dynearthsol2d -h` prints every option as

      --section.name TYPE (=default)
          Help text, possibly on several lines.

catalogue(exe) parses that; changed(cfg, options) lists the settings of an
input file that differ from DES's defaults, with their help text, for the
review's "Input" sections.

    des_options.py EXECUTABLE INPUT.cfg
"""
import configparser
import re
import subprocess
import sys

OPTION = re.compile(r'^\s+--(\S+)\s+(.*?)\s*(?:\(=(.*)\))?\s*$')
TRUE, FALSE = {'1', 'yes', 'true', 'on'}, {'0', 'no', 'false', 'off'}


def catalogue(exe):
    """{'section.name': (type, default or None, help)} from `exe -h`."""
    text = subprocess.run([str(exe), '-h'], capture_output=True, text=True, timeout=60).stdout
    return parse_help(text)


def parse_help(text):
    """Parse a captured help listing from the qualified executable."""
    options, key = {}, None
    for line in text.split('Build snapshot:')[0].splitlines():
        m = OPTION.match(line)
        if m and line.lstrip().startswith('--'):
            key = m.group(1); options[key] = [m.group(2), m.group(3), '']
        elif key and line.strip():
            options[key][2] = (options[key][2] + ' ' + line.strip()).strip()
    return {k: tuple(v) for k, v in options.items()}


def _numbers(text):
    try:
        return [float(x) for x in text.strip().strip('[]').replace(',', ' ').split()]
    except ValueError:
        return None


def same(value, default):
    """True when a configured value equals DES's default (booleans, numbers and lists compared by value)."""
    if default is None: return False
    v, d = value.strip().strip('"\'').lower(), default.strip().strip('"\'').lower()
    if v == d: return True
    if {v, d} <= TRUE or {v, d} <= FALSE: return True
    a, b = _numbers(v), _numbers(d)
    return a is not None and b is not None and len(a) == len(b) and all(abs(x - y) <= 1e-12*max(abs(x), abs(y)) for x, y in zip(a, b))


def read(cfg_path):
    cfg = configparser.ConfigParser(inline_comment_prefixes=('#',), strict=False); cfg.optionxform = str
    text = open(cfg_path).read()
    if not text.lstrip().startswith('['): text = '[top]\n' + text
    cfg.read_string(text)
    return {(f'{s}.{k}' if s != 'top' else k): v.strip() for s in cfg.sections() for k, v in cfg[s].items()}


def changed(cfg_path, options):
    """[(key, value, default, help)] for every setting that differs from DES's default, and keys DES does not know."""
    rows, unknown = [], []
    for key, value in read(cfg_path).items():
        if key not in options: unknown.append(key); continue
        _, default, help_text = options[key]
        if not same(value, default): rows.append((key, value, default, help_text))
    return rows, unknown


if __name__ == '__main__':
    rows, unknown = changed(sys.argv[2], catalogue(sys.argv[1]))
    for key, value, default, help_text in rows: print(f'{key} = {value}   (default {default})  {help_text}')
    if unknown: print('unknown to this DES:', ', '.join(unknown))
