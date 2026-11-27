#!/usr/bin/env python3
"""Check that every per-mesh field is distributed by the MPI backend or listed here.

Compute ranks hold only the fields that mpi_local_fields.hpp lists; any other field
stays unallocated there, so a kernel that reads it under MPI silently does nothing or
reads stale data. This check reads the sources (no build needed): every field that a
.cxx file allocates on a Variables with `new` must be in mpi_local_fields.hpp, or in
OUTSIDE below with the reason it never needs to reach the compute ranks.
"""
import argparse
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[3]

# Fields that stay off the compute ranks' field list, and why.
OUTSIDE = {
    # Each rank builds its own copy of the mesh topology when the mesh is partitioned.
    "connectivity": "per-rank topology", "bcflag": "per-rank topology",
    "bnormals": "per-rank topology", "neighbor": "per-rank topology",
    "top_elems": "per-rank topology",
    # Rank zero's mesh management: generation, remeshing and output.
    "connectivity_surface": "rank zero", "segment": "rank zero", "segflag": "rank zero",
    "regattr": "rank zero",
    # Markers stay on rank zero; compute ranks receive the material counts.
    "elemmarkers": "markers on rank zero", "markers_in_elem": "markers on rank zero",
    "hydrous_markers_in_elem": "markers on rank zero",
    # Scratch arrays that each mesh allocates for itself.
    "ntmp": "scratch", "ctmp": "scratch", "etmp": "scratch", "etmp_int": "scratch",
    "tmp_result": "scratch",
    # Lookup tables, identical everywhere.
    "log_table": "table", "sin_table": "table", "tan_table": "table",
    # Five accepted-step sums; the metadata carries them and every rank holds the reduced values.
    "fluid_balance": "reduced over the ranks",
}

ALLOCATION = re.compile(r"\b(?:var|master)\s*(?:\.|->)\s*(\w+)\s*=\s*(?:[^;]*[?:]\s*)?new\s+"
                        r"(?:double_vec|array_t|tensor_t|int_vec|uint_vec|int_vec2D|conn_t|"
                        r"segment_t|segflag_t|regattr_t|elem_cache)\b")


def without_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="repository root")
    root = parser.parse_args().root.resolve()

    listed = set(re.findall(r"\bF\((\w+)\)", (root / "mpi_local_fields.hpp").read_text()))
    allocated = {}
    for source in sorted(root.glob("*.cxx")):
        for name in ALLOCATION.findall(without_comments(source.read_text())):
            allocated.setdefault(name, source.name)

    missing = sorted(set(allocated) - listed - set(OUTSIDE))
    stale = sorted(set(OUTSIDE) - set(allocated))
    for name in missing:
        print(f"{allocated[name]}: '{name}' is neither in mpi_local_fields.hpp nor in OUTSIDE "
              f"of {Path(__file__).name}")
    for name in stale:
        print(f"'{name}' is in OUTSIDE but no longer allocated; remove it there")
    if missing or stale:
        print("Add a field the kernels read to mpi_local_fields.hpp (doc/mpi.md, "
              "\"Adding a distributed kernel\"); list one that never leaves rank zero in OUTSIDE.")
        return 1
    print(f"{len(allocated)} allocated per-mesh fields: "
          f"{len(set(allocated) & listed)} distributed, {len(set(allocated) & set(OUTSIDE))} outside")
    return 0


if __name__ == "__main__":
    sys.exit(main())
