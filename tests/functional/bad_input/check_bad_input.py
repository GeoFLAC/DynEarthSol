#!/usr/bin/env python3
"""Check that bad config files are rejected, and with the right exit code.

Each case edits a valid 2-D config (../plane_strain_elastic/plane_strain_elastic.cfg)
into something the parser must refuse, runs the executable on it, and checks the
exit code and a piece of the message. A rejection for some other reason does not
count. The exit codes are those of utils.hpp: 10 is a config error that is not a
bad value (missing key, missing file, malformed line), 11 is a bad value or an
unknown or repeated key.

A few cases must be accepted. They guard against the parser getting stricter
than the old one did: a string value is taken whole, with its commas, spaces and
quotes, and a decimal integer may carry leading zeros.

Usage: check_bad_input.py [--exe PATH]
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
BASE_CFG = HERE.parent / "plane_strain_elastic" / "plane_strain_elastic.cfg"

EXIT_CONFIG = 10
EXIT_CONFIG_VALUE = 11


def replace(old: str, new: str):
    """An edit that swaps one line; it fails loudly if the base cfg changed."""
    def edit(text: str) -> str:
        if old not in text:
            raise AssertionError(f"base cfg has no {old!r}; update the test")
        return text.replace(old, new, 1)
    return edit


def delete(old: str):
    return replace(old + "\n", "")


def prepend(extra: str):
    return lambda text: extra + "\n" + text


def append(extra: str):
    return lambda text: text + "\n" + extra + "\n"


# (name, edit, expected exit code, text the message must contain)
REJECTED = [
    ("unknown key",
     replace("[mesh]\n", "[mesh]\nbogus_key = 1\n"), EXIT_CONFIG_VALUE, "misspelled or retired"),
    ("unknown section",
     append("[nosuch]\nfoo = 1"), EXIT_CONFIG_VALUE, "misspelled or retired"),
    ("repeated scalar key",
     replace("resolution = 0.25", "resolution = 0.25\nresolution = 0.5"),
     EXIT_CONFIG_VALUE, "more than once"),
    ("repeated list key",
     replace("rho0 = [1]", "rho0 = [1]\nrho0 = [2]"), EXIT_CONFIG_VALUE, "more than once"),
    ("repeated key, once at the top of the file with its section",
     prepend("mat.rho0 = [5]"), EXIT_CONFIG_VALUE, "more than once"),
    ("key in double quotes",
     replace("modelname = plane-strain-elastic", '"modelname" = plane-strain-elastic'),
     EXIT_CONFIG_VALUE, "cannot be quoted"),
    ("key in single quotes",
     replace("modelname = plane-strain-elastic", "'modelname' = plane-strain-elastic"),
     EXIT_CONFIG_VALUE, "cannot be quoted"),
    ("section name in quotes",
     replace("[sim]", '["sim"]'), EXIT_CONFIG_VALUE, "cannot be quoted"),
    ("quoted key that repeats a list key",
     replace("rho0 = [1]", 'rho0 = [1]\n"rho0" = [2]'), EXIT_CONFIG_VALUE, "cannot be quoted"),
    ("text for a number",
     replace("resolution = 0.25", "resolution = abc"), EXIT_CONFIG_VALUE, "decimal number"),
    ("empty number",
     replace("resolution = 0.25", "resolution ="), EXIT_CONFIG_VALUE, "decimal number"),
    ("hexadecimal integer",
     replace("max_steps = 2", "max_steps = 0x10"), EXIT_CONFIG_VALUE, "decimal integer"),
    ("hexadecimal float",
     replace("xlength = 1", "xlength = 0x10"), EXIT_CONFIG_VALUE, "decimal number"),
    ("integer written as a float",
     replace("max_steps = 2", "max_steps = 1e3"), EXIT_CONFIG_VALUE, "decimal integer"),
    ("bool that is not a bool",
     replace("has_PT = no", "has_PT = maybe"), EXIT_CONFIG_VALUE, "true/false"),
    ("bool abbreviation",
     replace("has_PT = no", "has_PT = t"), EXIT_CONFIG_VALUE, "true/false"),
    ("missing required key",
     delete("xlength = 1"), EXIT_CONFIG, "required but missing"),
    ("list without brackets",
     replace("rho0 = [1]", "rho0 = 1"), EXIT_CONFIG_VALUE, "must be a list"),
    ("empty list value",
     replace("rho0 = [1]", "rho0 ="), EXIT_CONFIG_VALUE, "must be a list"),
    ("list not closed",
     replace("rho0 = [1]", "rho0 = [1"), EXIT_CONFIG_VALUE, "must be a list"),
    ("list continued on the next line",
     replace("rho0 = [1]", "rho0 = [1,\n        2]"), EXIT_CONFIG, "invalid line"),
    ("list continued after a comment",
     replace("rho0 = [1]", "rho0 = [1,  # crust\n 2]\nalpha = [2e-5]"),
     EXIT_CONFIG, "invalid line"),
    ("list of the wrong length",
     replace("rho0 = [1]", "rho0 = [1, 2]"), EXIT_CONFIG_VALUE, "must list 1 number"),
    ("text in a list",
     replace("rho0 = [1]", "rho0 = [abc]"), EXIT_CONFIG_VALUE, "must list 1 number"),
    ("list inside a list",
     replace("rho0 = [1]", "rho0 = [[1]]"), EXIT_CONFIG_VALUE, "must be a list"),
    ("quoted list element",
     replace("rho0 = [1]", 'rho0 = ["1"]'), EXIT_CONFIG_VALUE, "must be a list"),
    ("empty list element",
     replace("rho0 = [1]", "rho0 = [1,,]"), EXIT_CONFIG_VALUE, "must be a list"),
    ("list element with a blank in it",
     replace("rho0 = [1]", "rho0 = [1 2]"), EXIT_CONFIG_VALUE, "must be a list"),
    ("list of only a comma",
     replace("rho0 = [1]", "rho0 = [,]"), EXIT_CONFIG_VALUE, "must be a list"),
    ("empty bool",
     replace("has_PT = no", "has_PT ="), EXIT_CONFIG_VALUE, "true/false"),
    ("string value opening with three double quotes",
     replace("modelname = plane-strain-elastic", 'modelname = """abc'),
     EXIT_CONFIG_VALUE, "three quotes"),
    ("string value opening with three single quotes",
     replace("modelname = plane-strain-elastic", "modelname = '''abc"),
     EXIT_CONFIG_VALUE, "three quotes"),
    ("string value opening with [",
     replace("modelname = plane-strain-elastic", "modelname = [abc"),
     EXIT_CONFIG_VALUE, "not a list"),
    ("line without an equals sign",
     replace("[bc]\n", "[bc]\nvbc_x0\n"), EXIT_CONFIG, "invalid line"),
    ("section header not closed",
     replace("[bc]", "[bc"), EXIT_CONFIG, "invalid line"),
    ("key with no name",
     replace("[bc]\n", "[bc]\n= 1\n"), EXIT_CONFIG_VALUE, "unrecognised option"),
]

# (name, edit, check on the files the run leaves in its directory)
ACCEPTED = [
    ("unchanged config", lambda t: t, None),
    ("comment after a value and after a list",
     replace("rho0 = [1]", "rho0 = [1]   # density"), None),
    ("leading zeros on a decimal integer",
     replace("max_steps = 2", "max_steps = 02"), None),
    ("list with a trailing comma and blanks",
     replace("rho0 = [1]", "rho0 = [ 1 , ]"), None),
    ("a single quote that is never closed",
     replace("modelname = plane-strain-elastic", "modelname = 'abc"), "'abc"),
    ("string kept whole: commas, spaces and quotes",
     replace("modelname = plane-strain-elastic", 'modelname = "a,b  c"'), '"a,b  c"'),
]


def run(exe: Path, cfg_text: str, workdir: Path, filename: str = "case.cfg"):
    cfg = workdir / filename
    cfg.write_text(cfg_text)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    return subprocess.run([str(exe), str(cfg)], cwd=workdir, env=env, timeout=120,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--exe", type=Path, default=REPO_ROOT / "dynearthsol2d")
    exe = parser.parse_args().exe.resolve()
    if not exe.is_file():
        raise FileNotFoundError(f"executable not found: {exe}")
    base = BASE_CFG.read_text()

    failures = []

    def check(name, ok, detail):
        print(f"{'ok  ' if ok else 'FAIL'} {name}")
        if not ok:
            failures.append(f"{name}: {detail}")

    for name, edit, code, text in REJECTED:
        with tempfile.TemporaryDirectory(prefix="des-bad-input-") as tmp:
            result = run(exe, edit(base), Path(tmp))
            tail = result.stdout.strip().splitlines()[-4:]
            check(f"rejects: {name}", result.returncode == code and text in result.stdout,
                  f"want exit {code} and {text!r}, got exit {result.returncode}: {tail}")

    with tempfile.TemporaryDirectory(prefix="des-bad-input-") as tmp:
        missing = subprocess.run([str(exe), str(Path(tmp) / "nosuch.cfg")], cwd=tmp,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        check("rejects: file that does not exist", missing.returncode == EXIT_CONFIG,
              f"want exit {EXIT_CONFIG}, got {missing.returncode}")
        folder = subprocess.run([str(exe), tmp], cwd=tmp,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        check("rejects: a directory in place of the file", folder.returncode == EXIT_CONFIG,
              f"want exit {EXIT_CONFIG}, got {folder.returncode}")

    for name, edit, prefix in ACCEPTED:
        with tempfile.TemporaryDirectory(prefix="des-bad-input-") as tmp:
            result = run(exe, edit(base), Path(tmp))
            ok = result.returncode == 0
            detail = f"want exit 0, got {result.returncode}: {result.stdout.strip().splitlines()[-4:]}"
            if ok and prefix:
                names = [p.name for p in Path(tmp).iterdir() if p.name.startswith(prefix)]
                ok = bool(names)
                detail = f"no output file starts with {prefix!r}"
            check(f"accepts: {name}", ok, detail)

    if failures:
        print(f"\n{len(failures)} case(s) failed:", file=sys.stderr)
        for line in failures:
            print(f"  {line}", file=sys.stderr)
        return 1
    print(f"\nall {len(REJECTED) + 2 + len(ACCEPTED)} cases behave as required")
    return 0


if __name__ == "__main__":
    sys.exit(main())
