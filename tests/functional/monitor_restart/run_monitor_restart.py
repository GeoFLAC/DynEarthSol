#!/usr/bin/env python3

"""Check monitor CSV history preservation and fail-fast restart handling."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT / "tests/functional"))
from des_test_utils import run_des

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
    result = run_des(exe, case_dir, name, cfg, threads=threads, check=False)
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


def corrupt_kept_row_column_count(path: Path, *, has_header: bool) -> None:
    lines = path.read_bytes().splitlines(keepends=True)
    index = (1 if has_header else 0) + 1
    fields = lines[index].rstrip(b"\r\n").split(b",")
    if len(fields) < 3:
        raise AssertionError("fixture row has too few fields to corrupt")
    lines[index] = b",".join(fields[:-1]) + b"\n"
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

    cases = (
        ("malformed_row", (1,),
         lambda path: corrupt_kept_row_column_count(path, has_header=True),
         "wrong number of columns"),
        ("internal_duplicate", (0, 1),
         lambda path: insert_internal_duplicate(path, has_header=True),
         "duplicate or out-of-order rows for one step"),
        ("missing_checkpoint_row", (0, 1),
         lambda path: remove_checkpoint_row(path, has_header=True),
         "missing the row corresponding to the restart checkpoint"),
    )
    for name, point_ids, corrupt, diagnostic in cases:
        folder = run_root / name
        clone_baseline(baseline, folder)
        for point in point_ids:
            corrupt(point_paths(folder)[point])
        before = read_histories(folder)
        run_config(exe, restart_cfg, folder, "restart", threads=threads,
                   expected_error=diagnostic)
        if read_histories(folder) != before:
            raise AssertionError(f"{name}: monitor history changed before failure")

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

    for name, cfg in (("fresh_open_failure", fresh_cfg), ("restart_open_failure", restart_cfg)):
        folder = run_root / name
        clone_baseline(baseline, folder)
        cfg = cfg.replace("output_prefix = monitor_history", "output_prefix = missing/monitor_history")
        result = run_des(exe, folder, name, cfg, threads=threads, check=False)
        assert result.returncode == 20, (name, result.returncode, result.stdout)
        assert "cannot open" in result.stdout and "monitor file" in result.stdout


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run monitor CSV checkpoint-restart regressions."
    )
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE)
    parser.add_argument("--cfg", type=Path, default=DEFAULT_CFG)
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
) -> None:
    check_history_replacement(
        exe, template, run_root, write_header=True, threads=threads
    )
    check_history_replacement(
        exe, template, run_root, write_header=False, threads=threads
    )
    check_fail_fast(exe, template, run_root, threads=threads)


def main() -> None:
    args = parse_args()
    exe = args.exe.expanduser().resolve()
    cfg = args.cfg.expanduser().resolve()
    if not exe.is_file():
        raise FileNotFoundError(f"executable not found: {exe}")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    template = cfg.read_text(encoding="ascii")

    if args.run_dir:
        run_root = args.run_dir.expanduser().resolve()
        run_root.mkdir(parents=True, exist_ok=True)
        run_all(
            exe,
            template,
            run_root,
            threads=args.threads,
        )
    else:
        with tempfile.TemporaryDirectory(prefix="des-monitor-restart-") as tmp:
            run_all(
                exe,
                template,
                Path(tmp),
                threads=args.threads,
            )

    print("monitor restart regression: PASS")


if __name__ == "__main__":
    main()
