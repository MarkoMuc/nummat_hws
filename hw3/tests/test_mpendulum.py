from math import exp, pi, sin, sqrt

import numpy as np
import pytest

from mpendulum import (
    DOPRI5,
    InitialValueProblem,
    PendulumParameters,
    harmonic_angle,
    harmonic_period,
    pendulum_angle,
    pendulum_period,
    solve,
)


def test_dopri5_solves_an_equation_with_known_exponential_solution():
    problem = InitialValueProblem(
        lambda _t, y, _p: y,
        np.array([1.0]),
        (0.0, 1.0),
    )
    solution = solve(
        problem,
        DOPRI5(rtol=1e-11, atol=1e-13, initial_step=0.2, max_step=0.3),
    )

    points = np.linspace(0.0, 1.0, 21)
    errors = np.abs(solution(points)[:, 0] - np.exp(points))

    assert solution.times[-1] == 1.0
    assert solution.rejected_steps >= 1
    assert np.max(errors) < 2e-9
    assert solution.values[-1, 0] == pytest.approx(exp(1.0), rel=2e-11)


def _complete_elliptic_k(modulus):
    arithmetic = 1.0
    geometric = sqrt(1.0 - modulus ** 2)
    while abs(arithmetic - geometric) > 2e-16 * arithmetic:
        arithmetic, geometric = (
            (arithmetic + geometric) / 2.0,
            sqrt(arithmetic * geometric),
        )
    return pi / (2.0 * arithmetic)


def test_period_matches_the_independent_elliptic_integral_formula():
    parameters = PendulumParameters(gravity=1.0, length=1.0)
    theta0 = pi / 2.0
    reference = 4.0 * _complete_elliptic_k(sin(theta0 / 2.0))

    result = pendulum_period(theta0, parameters=parameters)

    assert result.period == pytest.approx(reference, rel=8e-11)
    assert result.crossings[0].direction == result.crossings[1].direction == -1
    assert abs(result.crossings[0].state[0]) < 2e-12
    assert abs(result.crossings[1].state[0]) < 2e-12


def test_harmonic_formula_and_angle_interface_respect_initial_data():
    parameters = PendulumParameters(gravity=4.0, length=1.0)
    theta0 = 0.3
    omega0 = -0.2

    assert harmonic_angle(0.0, theta0, omega0, parameters=parameters) == theta0
    assert harmonic_angle(
        harmonic_period(parameters), theta0, omega0, parameters=parameters
    ) == pytest.approx(theta0, abs=2e-15)
    assert pendulum_angle(0.0, theta0, omega0, parameters=parameters) == theta0
    assert abs(
        pendulum_angle(
            harmonic_period(parameters) / 4.0,
            1e-3,
            parameters=parameters,
            method=DOPRI5(rtol=1e-11, atol=1e-13),
        )
    ) < 2e-9
