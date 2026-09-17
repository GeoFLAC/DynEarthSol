#!/usr/bin/env python3
"""Reference solutions and monitor readers for the local EP/RSF benchmarks."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

from run_simple_shear_benchmark import BenchmarkCase


VX_TOP_M_S = 1.0e-5
HEIGHT_M = 1.0
SHEAR_MODULUS_PA = 200.0e6
COHESION_PA = 1.0e6
SHEAR_DT_S = 0.01
VELOCITY_FLOOR_M_S = 1.0e-16


@dataclass(frozen=True)
class MonitorData:
    time_s: list[float]
    stress_by_element: list[list[float]]
    friction_by_element: list[list[float] | None]
    state_by_element: list[list[float] | None]

    @property
    def mean_abs_stress(self) -> list[float]:
        return [
            sum(values) / len(values)
            for values in zip(*self.stress_by_element)
        ]


@dataclass(frozen=True)
class ErrorMetrics:
    samples: int
    mean_fraction: float
    max_fraction: float
    final_fraction: float


def slip_rate(time_s: float) -> float:
    """Return 2 w strain_rate_II for the deforming two-triangle cell."""
    gamma = VX_TOP_M_S * time_s / HEIGHT_M
    return VX_TOP_M_S / math.sqrt((1.0 - gamma) ** 2 + 1.0)


def _integrating_factor(time_s: float, dc_m: float) -> float:
    gamma = VX_TOP_M_S * time_s / HEIGHT_M
    return (HEIGHT_M / dc_m) * (
        math.asinh(1.0) - math.asinh(1.0 - gamma)
    )


def aging_state(
    time_s: list[float],
    dc_m: float,
    v0_m_s: float,
) -> list[float]:
    """Solve the aging law by the paper's integrating-factor quadrature."""
    if not time_s:
        return []
    if any(right < left for left, right in zip(time_s, time_s[1:])):
        raise ValueError("Reference times must be nondecreasing.")

    max_time_s = time_s[-1]
    grid_steps = int(math.ceil(max_time_s / SHEAR_DT_S - 1.0e-12))
    integrals = [0.0] * (grid_steps + 1)
    previous_exp = 1.0
    for step in range(1, grid_steps + 1):
        grid_time_s = step * SHEAR_DT_S
        current_exp = math.exp(_integrating_factor(grid_time_s, dc_m))
        integrals[step] = (
            integrals[step - 1]
            + 0.5 * (previous_exp + current_exp) * SHEAR_DT_S
        )
        previous_exp = current_exp

    theta0_s = dc_m / v0_m_s
    result: list[float] = []
    for sample_time_s in time_s:
        position = sample_time_s / SHEAR_DT_S
        lower = min(int(math.floor(position)), grid_steps)
        upper = min(lower + 1, grid_steps)
        if upper == lower:
            integral = integrals[lower]
        else:
            fraction = position - lower
            integral = integrals[lower] + fraction * (
                integrals[upper] - integrals[lower]
            )
        phi = _integrating_factor(sample_time_s, dc_m)
        result.append(math.exp(-phi) * (theta0_s + integral))
    return result


def yield_stress_from_mu(mu: float) -> float:
    """Return the zero-pressure Mohr-Coulomb shear strength."""
    phi = math.atan(max(mu, 1.0e-6))
    sin_phi = math.sin(phi)
    n_phi = (1.0 + sin_phi) / (1.0 - sin_phi)
    return 2.0 * COHESION_PA * math.sqrt(n_phi) / (1.0 + n_phi)


def analytical_stress(
    case: BenchmarkCase,
    time_s: list[float],
) -> list[float]:
    """Return the paper's small-strain reference stress."""
    elastic = [
        SHEAR_MODULUS_PA * VX_TOP_M_S * value / HEIGHT_M
        for value in time_s
    ]
    mu0 = math.tan(math.radians(case.friction_angle_deg))

    if case.group == "ep":
        strength = yield_stress_from_mu(mu0)
        return [min(value, strength) for value in elastic]

    velocities = [slip_rate(value) for value in time_s]
    if case.state_var_model == 0:
        friction = [
            mu0
            + (case.direct_a - case.evolution_b)
            * math.log(velocity / case.characteristic_velocity)
            for velocity in velocities
        ]
    else:
        states = aging_state(
            time_s,
            case.characteristic_distance,
            case.characteristic_velocity,
        )
        friction = [
            mu0
            + case.direct_a
            * math.log(velocity / case.characteristic_velocity)
            + case.evolution_b
            * math.log(
                case.characteristic_velocity
                * state
                / case.characteristic_distance
            )
            for velocity, state in zip(velocities, states)
        ]

    return [
        min(value, yield_stress_from_mu(mu))
        for value, mu in zip(elastic, friction)
    ]


def analytical_state_ratio(
    case: BenchmarkCase,
    time_s: list[float],
) -> list[float]:
    states = aging_state(
        time_s,
        case.characteristic_distance,
        case.characteristic_velocity,
    )
    return [
        state / (case.characteristic_distance / slip_rate(value))
        for value, state in zip(time_s, states)
    ]


def load_monitor_case(case_dir: Path) -> MonitorData:
    paths = sorted(case_dir.glob("monitor_point_*.csv"))
    if len(paths) != 2:
        raise FileNotFoundError(f"Expected two monitor CSVs in {case_dir}")

    times: list[list[float]] = []
    stresses: list[list[float]] = []
    frictions: list[list[float] | None] = []
    states: list[list[float] | None] = []
    for path in paths:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = sorted(
                csv.DictReader(handle),
                key=lambda row: float(row["time_s"]),
            )
        if not rows:
            raise ValueError(f"No monitor rows found in {path}")
        header = rows[0].keys()
        times.append([float(row["time_s"]) for row in rows])
        stresses.append([abs(float(row["stress_2"])) for row in rows])
        frictions.append(
            [float(row["dynamic_friction"]) for row in rows]
            if "dynamic_friction" in header
            else None
        )
        states.append(
            [float(row["state_variable"]) for row in rows]
            if "state_variable" in header
            else None
        )

    reference_time = times[0]
    if len(times[1]) != len(reference_time):
        raise ValueError(f"Monitor time-axis length mismatch in {case_dir}")
    for left, right in zip(reference_time, times[1]):
        if abs(left - right) > 1.0e-10:
            raise ValueError(f"Monitor time-axis mismatch in {case_dir}")

    return MonitorData(reference_time, stresses, frictions, states)


def pointwise_error_metrics(
    numerical: list[float],
    reference: list[float],
) -> ErrorMetrics:
    if len(numerical) != len(reference):
        raise ValueError("Numerical and reference series have different lengths.")
    errors = [
        abs(num - ref) / abs(ref)
        for num, ref in zip(numerical[1:], reference[1:])
    ]
    if not errors:
        raise ValueError("The t=0 sample is the only available sample.")
    return ErrorMetrics(
        len(errors),
        sum(errors) / len(errors),
        max(errors),
        errors[-1],
    )


def recovered_rate_errors(
    case: BenchmarkCase,
    data: MonitorData,
) -> list[float]:
    """Recover the constitutive rate from monitored friction and state."""
    if case.group == "ep":
        return []
    mu0 = math.tan(math.radians(case.friction_angle_deg))
    errors: list[float] = []
    for friction, state in zip(
        data.friction_by_element,
        data.state_by_element,
    ):
        if friction is None:
            raise ValueError("Dynamic friction is required to recover slip rate.")
        element_errors: list[float] = []
        for index in range(1, len(data.time_s)):
            if case.state_var_model == 0:
                velocity = case.characteristic_velocity * math.exp(
                    (friction[index] - mu0)
                    / (case.direct_a - case.evolution_b)
                )
            else:
                if state is None:
                    raise ValueError("State variable is required to recover slip rate.")
                evolution = case.evolution_b * math.log(
                    case.characteristic_velocity
                    * state[index]
                    / case.characteristic_distance
                )
                velocity = case.characteristic_velocity * math.exp(
                    (friction[index] - mu0 - evolution) / case.direct_a
                )
            reference = slip_rate(data.time_s[index])
            element_errors.append(abs(velocity - reference) / reference)
        errors.append(max(element_errors))
    return errors


def healing_reference(max_step: int) -> list[float]:
    dc_m = 1.0e-2
    theta_s = dc_m / 4.0e-9
    result = [theta_s]
    for _ in range(max_step):
        theta_s += (
            1.0 - VELOCITY_FLOOR_M_S * theta_s / dc_m
        ) * 259_200.0
        theta_s = min(max(theta_s, 1.0e-12), 1.0e12)
        result.append(theta_s)
    return result
