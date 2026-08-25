#!/usr/bin/env python3

"""Check monitor CSV history preservation and fail-fast restart handling."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_EXE = REPO_ROOT / "dynearthsol2d"
DEFAULT_CFG = HERE / "monitor_restart_base.cfg"
POINT_FILES = tuple(f"monitor_history_point_{point}.csv" for point in range(2))


def render_cfg(
    template: str,
    *,
    restarting: bool,
    write_header: bool,
    restart_frame: int = 1,
) -> str:
    replacements = {
        "__IS_RESTARTING__": "yes" if restarting else "no",
        "__WRITE_HEADER__": "yes" if write_header else "no",
        "__RESTART_FRAME__": str(restart_frame),
    }
    rendered = template
    for token, value in replacements.items():
        if token not in rendered:
            raise AssertionError(f"missing config token {token}")
        rendered = rendered.replace(token, value)
    if "__" in rendered:
        raise AssertionError("unexpanded config token remains")
    return rendered


def run_config(
    exe: Path,
    cfg: str,
    case_dir: Path,
    name: str,
    *,
    threads: int,
    expected_error: str | None = None,
) -> str:
    case_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = case_dir / f"{name}.cfg"
    cfg_path.write_text(cfg, encoding="ascii")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    result = subprocess.run(
        [str(exe), str(cfg_path)],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    (case_dir / f"{name}.log").write_text(result.stdout, encoding="utf-8")
    if expected_error is None:
        if result.returncode != 0:
            print(result.stdout, file=sys.stderr)
            raise RuntimeError(f"{name}: DynEarthSol exited with {result.returncode}")
    else:
        if result.returncode != 22:
            print(result.stdout, file=sys.stderr)
            raise AssertionError(
                f"{name}: expected restart-mismatch status 22, got {result.returncode}"
            )
        if expected_error not in result.stdout:
            print(result.stdout, file=sys.stderr)
            raise AssertionError(
                f"{name}: missing expected diagnostic {expected_error!r}"
            )
    return result.stdout


def point_paths(case_dir: Path) -> tuple[Path, ...]:
    return tuple(case_dir / filename for filename in POINT_FILES)


def read_histories(case_dir: Path) -> tuple[bytes, ...]:
    paths = point_paths(case_dir)
    for path in paths:
        if not path.is_file():
            raise AssertionError(f"missing monitor history: {path}")
    return tuple(path.read_bytes() for path in paths)


def read_steps(path: Path, *, has_header: bool) -> list[int]:
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.reader(stream))
    if has_header:
        if not rows or rows[0][:2] != ["step", "time_s"]:
            raise AssertionError(f"{path}: missing expected header")
        rows = rows[1:]
    return [int(row[0]) for row in rows]


def assert_steps(case_dir: Path, expected: list[int], *, has_header: bool) -> None:
    for path in point_paths(case_dir):
        steps = read_steps(path, has_header=has_header)
        if steps != expected:
            raise AssertionError(f"{path}: expected steps {expected}, got {steps}")


def clone_baseline(baseline: Path, destination: Path) -> None:
    shutil.copytree(baseline, destination)


def insert_checkpoint_duplicate(path: Path, *, has_header: bool) -> None:
    lines = path.read_bytes().splitlines(keepends=True)
    row_offset = 1 if has_header else 0
    checkpoint_index = row_offset + 2
    lines.insert(checkpoint_index + 1, lines[checkpoint_index])
    path.write_bytes(b"".join(lines))


def insert_internal_duplicate(path: Path, *, has_header: bool) -> None:
    lines = path.read_bytes().splitlines(keepends=True)
    row_offset = 1 if has_header else 0
    internal_index = row_offset + 1
    lines.insert(internal_index + 1, lines[internal_index])
    path.write_bytes(b"".join(lines))


def remove_checkpoint_row(path: Path, *, has_header: bool) -> None:
    lines = path.read_bytes().splitlines(keepends=True)
    row_offset = 1 if has_header else 0
    del lines[row_offset + 2]
    path.write_bytes(b"".join(lines))


def corrupt_last_row_column_count(path: Path) -> None:
    lines = path.read_bytes().splitlines(keepends=True)
    last = lines[-1].rstrip(b"\r\n")
    newline = b"\r\n" if lines[-1].endswith(b"\r\n") else b"\n"
    fields = last.split(b",")
    if len(fields) < 3:
        raise AssertionError("fixture row has too few fields to corrupt")
    lines[-1] = b",".join(fields[:-1]) + newline
    path.write_bytes(b"".join(lines))


def corrupt_last_row_value(path: Path) -> None:
    lines = path.read_bytes().splitlines(keepends=True)
    last = lines[-1].rstrip(b"\r\n")
    newline = b"\r\n" if lines[-1].endswith(b"\r\n") else b"\n"
    fields = last.split(b",")
    fields[-1] = b"not-a-number"
    lines[-1] = b",".join(fields) + newline
    path.write_bytes(b"".join(lines))


def check_history_replacement(
    exe: Path,
    template: str,
    run_root: Path,
    *,
    write_header: bool,
    threads: int,
) -> None:
    suffix = "header" if write_header else "no_header"
    baseline = run_root / f"baseline_{suffix}"
    fresh_cfg = render_cfg(template, restarting=False, write_header=write_header)
    restart_cfg = render_cfg(template, restarting=True, write_header=write_header)
    run_config(exe, fresh_cfg, baseline, "fresh", threads=threads)
    assert_steps(baseline, [0, 1, 2, 3, 4], has_header=write_header)
    expected = read_histories(baseline)

    ordinary = run_root / f"ordinary_{suffix}"
    clone_baseline(baseline, ordinary)
    run_config(exe, restart_cfg, ordinary, "restart", threads=threads)
    assert_steps(ordinary, [0, 1, 2, 3, 4], has_header=write_header)
    if read_histories(ordinary) != expected:
        raise AssertionError(
            f"{suffix}: restart did not reproduce the uninterrupted CSV histories"
        )

    frame_zero = run_root / f"frame_zero_{suffix}"
    clone_baseline(baseline, frame_zero)
    frame_zero_cfg = render_cfg(
        template,
        restarting=True,
        write_header=write_header,
        restart_frame=0,
    )
    run_config(exe, frame_zero_cfg, frame_zero, "restart", threads=threads)
    assert_steps(frame_zero, [0, 1, 2, 3, 4], has_header=write_header)
    if read_histories(frame_zero) != expected:
        raise AssertionError(
            f"{suffix}: frame-0 restart duplicated or changed the initial row"
        )

    duplicate = run_root / f"checkpoint_duplicate_{suffix}"
    clone_baseline(baseline, duplicate)
    for path in point_paths(duplicate):
        insert_checkpoint_duplicate(path, has_header=write_header)
    run_config(exe, restart_cfg, duplicate, "restart", threads=threads)
    assert_steps(duplicate, [0, 1, 2, 3, 4], has_header=write_header)
    if read_histories(duplicate) != expected:
        raise AssertionError(
            f"{suffix}: exact checkpoint duplicates were not canonicalized"
        )

    new_history = run_root / f"new_history_{suffix}"
    clone_baseline(baseline, new_history)
    for path in point_paths(new_history):
        path.unlink()
    run_config(exe, restart_cfg, new_history, "restart", threads=threads)
    assert_steps(new_history, [3, 4], has_header=write_header)


def check_fail_fast(
    exe: Path, template: str, run_root: Path, *, threads: int
) -> None:
    fresh_cfg = render_cfg(template, restarting=False, write_header=True)
    restart_cfg = render_cfg(template, restarting=True, write_header=True)
    baseline = run_root / "baseline_fail_fast"
    run_config(exe, fresh_cfg, baseline, "fresh", threads=threads)

    malformed = run_root / "malformed_tail"
    clone_baseline(baseline, malformed)
    corrupt_last_row_column_count(point_paths(malformed)[1])
    before = read_histories(malformed)
    run_config(
        exe,
        restart_cfg,
        malformed,
        "restart",
        threads=threads,
        expected_error="wrong number of columns",
    )
    if read_histories(malformed) != before:
        raise AssertionError("malformed point file partially changed another history")

    non_numeric = run_root / "non_numeric_tail"
    clone_baseline(baseline, non_numeric)
    corrupt_last_row_value(point_paths(non_numeric)[1])
    before = read_histories(non_numeric)
    run_config(
        exe,
        restart_cfg,
        non_numeric,
        "restart",
        threads=threads,
        expected_error="non-numeric data column",
    )
    if read_histories(non_numeric) != before:
        raise AssertionError("non-numeric point file was modified before failure")

    internal_duplicate = run_root / "internal_duplicate"
    clone_baseline(baseline, internal_duplicate)
    for path in point_paths(internal_duplicate):
        insert_internal_duplicate(path, has_header=True)
    before = read_histories(internal_duplicate)
    run_config(
        exe,
        restart_cfg,
        internal_duplicate,
        "restart",
        threads=threads,
        expected_error="duplicate row before the restart cutoff",
    )
    if read_histories(internal_duplicate) != before:
        raise AssertionError("ambiguous duplicate history was modified before failure")

    missing_checkpoint = run_root / "missing_checkpoint_row"
    clone_baseline(baseline, missing_checkpoint)
    for path in point_paths(missing_checkpoint):
        remove_checkpoint_row(path, has_header=True)
    before = read_histories(missing_checkpoint)
    run_config(
        exe,
        restart_cfg,
        missing_checkpoint,
        "restart",
        threads=threads,
        expected_error="missing the row corresponding to the restart checkpoint",
    )
    if read_histories(missing_checkpoint) != before:
        raise AssertionError("history missing its checkpoint row was modified")

    partial = run_root / "partial_point_set"
    clone_baseline(baseline, partial)
    point_paths(partial)[1].unlink()
    first_before = point_paths(partial)[0].read_bytes()
    run_config(
        exe,
        restart_cfg,
        partial,
        "restart",
        threads=threads,
        expected_error="same history as the other monitor points",
    )
    if point_paths(partial)[0].read_bytes() != first_before:
        raise AssertionError("complete point history was modified before partial-set failure")
    if point_paths(partial)[1].exists():
        raise AssertionError("missing point history was silently recreated before failure")


def check_fresh_byte_parity(
    exe: Path,
    reference_exe: Path,
    template: str,
    run_root: Path,
    *,
    threads: int,
) -> None:
    for write_header in (True, False):
        suffix = "header" if write_header else "no_header"
        cfg = render_cfg(template, restarting=False, write_header=write_header)
        reference = run_root / f"fresh_reference_{suffix}"
        candidate = run_root / f"fresh_candidate_{suffix}"
        run_config(reference_exe, cfg, reference, "fresh", threads=threads)
        run_config(exe, cfg, candidate, "fresh", threads=threads)
        if read_histories(candidate) != read_histories(reference):
            raise AssertionError(
                f"{suffix}: fresh monitor CSV changed relative to the reference executable"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run monitor CSV checkpoint-restart regressions."
    )
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
    parser.add_argument("--reference-exe", type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Keep outputs below this directory (default: temporary directory).",
    )
    return parser.parse_args()


def run_all(
    exe: Path,
    template: str,
    run_root: Path,
    *,
    threads: int,
    reference_exe: Path | None,
) -> None:
    check_history_replacement(
        exe, template, run_root, write_header=True, threads=threads
    )
    check_history_replacement(
        exe, template, run_root, write_header=False, threads=threads
    )
    check_fail_fast(exe, template, run_root, threads=threads)
    if reference_exe is not None:
        check_fresh_byte_parity(
            exe, reference_exe, template, run_root, threads=threads
        )


def main() -> None:
    args = parse_args()
    exe = args.exe.expanduser().resolve()
    cfg = args.cfg.expanduser().resolve()
    if not exe.is_file():
        raise FileNotFoundError(f"executable not found: {exe}")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    reference_exe = None
    if args.reference_exe is not None:
        reference_exe = args.reference_exe.expanduser().resolve()
        if not reference_exe.is_file():
            raise FileNotFoundError(f"reference executable not found: {reference_exe}")
    template = cfg.read_text(encoding="ascii")

    if args.run_dir:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(
            exe,
            template,
            run_root,
            threads=args.threads,
            reference_exe=reference_exe,
        )
    else:
        with tempfile.TemporaryDirectory(prefix="des-monitor-restart-") as tmp:
            run_all(
                exe,
                template,
                Path(tmp),
                threads=args.threads,
                reference_exe=reference_exe,
            )

    print("monitor restart regression: PASS")


if __name__ == "__main__":
    main()
