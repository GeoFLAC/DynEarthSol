#!/usr/bin/env python3
"""Plot and check the homogeneous PT tutorial runs, using DES's output reader.

Run from the directory containing the model outputs:
    python analyze_pt.py elastic elastic-0 elastic-1 --source ..
    python analyze_pt.py maxwell maxwell-0 maxwell-1 --source ..
"""

import argparse
import csv
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def save_figure(figure, stem):
    for extension in ("svg", "png"):
        path = Path(f"{stem}.{extension}")
        figure.savefig(path, dpi=160,
                       metadata={"Date": None} if extension == "svg" else None)
        if extension == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")


def reference(kind, count):
    """Confined elastic response or an independent Maxwell CN recurrence."""
    bulk = shear = 1e7
    load, dt, viscosity = -1e5, 0.01, 1e6
    strain = np.zeros(count)
    stress = np.zeros(count)
    increments = np.zeros(count)
    if kind == "elastic":
        strain[1:] = load / (bulk + 4 * shear / 3)
        stress[1:] = load
        increments[1] = strain[1]
    else:
        a = shear * dt / (2 * viscosity)
        decay = (1 - a) / (1 + a)
        modulus = 4 * shear / (3 * (1 + a))
        mean = deviatoric = 0.0
        for step in range(1, count):
            increment = (load - mean - decay * deviatoric) / (bulk + modulus)
            mean += bulk * increment
            deviatoric = decay * deviatoric + modulus * increment
            increments[step] = increment
            strain[step] = strain[step - 1] + increment
            stress[step] = mean + deviatoric
    return strain, stress, increments


def read_model(Dynearthsol, prefix, kind):
    model = Dynearthsol(str(prefix))
    expected_steps = 3 if kind == "elastic" else 100
    if list(model.steps) != list(range(expected_steps + 1)):
        raise ValueError(f"{prefix}: expected every step from 0 to {expected_steps}")
    # The .info time column is seconds; fixed_dt is 0.01 s in all four inputs.
    time = np.asarray(model.time)
    np.testing.assert_allclose(time, np.arange(expected_steps + 1) * 0.01,
                               rtol=1e-12, atol=1e-14)
    ref_strain, ref_stress, increments = reference(kind, len(model.frames))
    strains, stresses, displacement = [], [], []
    top_displacement = 0.0
    max_error = 0.0
    for index, frame in enumerate(model.frames):
        fields = {name: model.read_field(frame, name)
                  for name in ("strain", "stress", "coordinate", "velocity")}
        for name, field in fields.items():
            if not np.isfinite(field).all():
                raise ValueError(f"{prefix}: nonfinite {name} at frame {frame}")
        strain = fields["strain"][:, model.ndims - 1]
        stress = fields["stress"][:, model.ndims - 1]
        coordinate = fields["coordinate"]
        velocity = fields["velocity"][:, model.ndims - 1]
        reference_velocity = (coordinate[:, -1] + 1) * increments[index] / 0.01
        max_error = max(max_error,
                        np.max(np.abs(strain - ref_strain[index])) / 0.01,
                        np.max(np.abs(stress - ref_stress[index])) / 1e5,
                        np.max(np.abs(velocity - reference_velocity)) / 0.5)
        if index:
            top = np.isclose(coordinate[:, -1], 0)
            top_displacement += float(np.mean(velocity[top])) * (time[index] - time[index - 1])
        strains.append(float(np.mean(strain)))
        stresses.append(float(np.mean(stress)))
        displacement.append(top_displacement)
    if max_error > 5e-5:
        raise ValueError(f"{prefix}: normalized field error {max_error:.6g} exceeds 5e-5")
    records = []
    pattern = (r"PT step=(\d+) status=converged iterations=(\d+) residual=(\S+) "
               r"initial=(\S+) threshold=(\S+) warm_initial=(\S+) evaluations=(\d+) "
               r"stage=before_physical_updates")
    for match in re.finditer(pattern, prefix.with_suffix(".log").read_text()):
        step, iterations, residual, initial, threshold, warm, evaluations = match.groups()
        record = [int(step), int(iterations), float(residual), float(threshold),
                  int(evaluations), float(initial), float(warm)]
        if not np.isfinite(record).all() or record[2] > record[3]:
            raise ValueError(f"{prefix}: invalid converged record at step {step}")
        records.append(record)
    if [record[0] for record in records] != list(range(1, expected_steps + 1)):
        raise ValueError(f"{prefix}: missing PT diagnostic records; keep the solver .log")
    if kind == "elastic" and any(record[1] != 1 for record in records[1:]):
        raise ValueError(f"{prefix}: the held elastic steps should need no relaxation update")
    with prefix.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["step", "time_s", "mean_strain_zz", "mean_stress_zz_Pa",
                         "equivalent_top_displacement_m", "reference_strain_zz"])
        writer.writerows(zip(model.steps, time, strains, stresses, displacement, ref_strain))
    print(f"{prefix}: {model.ndims}D, max normalized field error={max_error:.6g}; "
          f"final strain_zz={strains[-1]:.9g}, top displacement={displacement[-1]:.9g} m")
    return time, np.array(strains), np.array(stresses), ref_strain, np.array(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("elastic", "maxwell"))
    parser.add_argument("models", nargs="+", type=Path)
    parser.add_argument("--source", required=True, type=Path,
                        help="DynEarthSol directory containing Dynearthsol.py")
    args = parser.parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    from Dynearthsol import Dynearthsol

    plt.rcParams.update({"font.size": 11, "svg.hashsalt": "des-pt-tutorial",
                         "svg.fonttype": "none"})
    figure, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    all_records = []
    for index, prefix in enumerate(args.models):
        time, strain, stress, ref_strain, records = read_model(Dynearthsol, prefix, args.kind)
        axes[0].plot(time, -100 * strain, "o-", ms=3, markevery=max(1, len(time) // 10), label=str(prefix))
        axes[1].plot(time, -stress / 1000, "o-", ms=3, markevery=max(1, len(time) // 10), label=str(prefix))
        if index == len(args.models) - 1:
            axes[0].plot(time, -100 * ref_strain, "k--", lw=1, label="Elastic reference" if args.kind == "elastic" else "Discrete Maxwell reference")
        all_records.append((prefix, records))
    axes[0].set(ylabel="Axial compression (%)", xlabel="Physical time (s)")
    axes[1].set(ylabel="Compressive axial stress (kPa)", xlabel="Physical time (s)")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    save_figure(figure, f"{args.kind}-response")
    plt.close(figure)

    if args.kind == "elastic":
        figure, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
        for prefix, records in all_records:
            axes[0].plot(records[:, 0], records[:, 1], "o-", label=str(prefix))
            axes[1].semilogy(records[:, 0], records[:, 2] / records[:, 3], "o-", label=str(prefix))
        axes[0].set(ylabel="PT loop passes", xlabel="Physical step", xticks=[1, 2, 3])
        axes[1].axhline(1, color="black", ls="--", lw=1, label="Acceptance limit")
        axes[1].set(ylabel="Accepted residual / threshold", xlabel="Physical step", xticks=[1, 2, 3])
        for axis in axes:
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
        save_figure(figure, "elastic-convergence")
        plt.close(figure)


if __name__ == "__main__":
    main()
