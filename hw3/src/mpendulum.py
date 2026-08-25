from dataclasses import dataclass
from math import isfinite, pi, sqrt

import numpy as np


# Parts of the code were adapted from Vaje 16

class InitialValueProblem:
    """Initial value problem u'(t) = f(t, u(t),p), u(t0) = u0.

        Where t is the independent variable, u(t) is the dependent variable, and p are the parameters.
    """

    def __init__(self, function, initial_value, time_step, parameters=None):
        # Just validate the inputs
        if initial_value.ndim != 1 or initial_value.size == 0:
            raise ValueError(f"Initial value array should be an one dimensional non-empty array.")

        t0, t1 = (float(value) for value in time_step)
        if not t0 < t1:
            raise ValueError("DOPRI 5 demands t0 < t1.")

        self.function = function
        self.initial_value = np.array(initial_value, dtype=float, copy=True)
        self.time_span = (t0, t1)
        self.parameters = parameters


@dataclass(frozen=True)
class DOPRI5:
    """General setting for the adaptive DOPRI5 method we use"""

    rtol: float = 1e-10
    atol: float = 1e-12
    initial_step: float | None = None
    max_step: float = float("inf")
    min_step: float = 0.0
    safety: float = 0.9
    min_factor: float = 0.2
    max_factor: float = 5.0
    max_steps: int = 1_000_000

    def __post_init__(self):
        if self.rtol < 0.0 or self.atol < 0.0:
            raise ValueError("Tolerance must be positive.")

        if self.min_step > self.max_step:
            raise ValueError("min_step < max_step.")


# We will also use the cubic Hermite interpolation as in Vaje 16 so we can plot better
class ODESolution:
    """Stores the numerical result and allows to interpolate on the solutino for plotting"""

    def __init__(
            self,
            problem,
            times,
            values,
            derivatives,
            accepted_steps,
            rejected_steps,
            function_evaluations,
    ):
        times = np.asarray(times, dtype=float)
        values = np.asarray(values, dtype=float)
        derivatives = np.asarray(derivatives, dtype=float)
        if values.shape != derivatives.shape:
            raise ValueError("Values and derivatives must match")
        if values.ndim != 2 or values.shape[0] != times.size:
            raise ValueError("Valyes and times must be the same size.")
        if values.shape[1] != problem.initial_value.size:
            raise ValueError("Values vector has changed in size.")
        if not np.all(np.diff(times) > 0.0):
            raise ValueError("Times must be strictly increasing.")

        times = times.copy()
        values = values.copy()
        derivatives = derivatives.copy()
        times.setflags(write=False)
        values.setflags(write=False)
        derivatives.setflags(write=False)

        self.problem = problem
        self.times = times
        self.values = values
        self.derivatives = derivatives
        self.accepted_steps = accepted_steps
        self.rejected_steps = rejected_steps
        self.function_evaluations = function_evaluations

    @property
    def step_sizes(self):
        """Returns the step size"""

        return np.diff(self.times)

    def _evaluate_scalar(self, time):
        # Last time in the interval
        if time == self.times[-1]:
            return self.values[-1].copy()

        # Find where it will be inserted
        index = int(np.searchsorted(self.times, time, side="right") - 1)
        left_time = self.times[index]
        right_time = self.times[index + 1]
        h = right_time - left_time
        s = (time - left_time) / h

        # As in Vaje 16
        h00 = 2.0 * s ** 3 - 3.0 * s ** 2 + 1.0
        h10 = s ** 3 - 2.0 * s ** 2 + s
        h01 = -2.0 * s ** 3 + 3.0 * s ** 2
        h11 = s ** 3 - s ** 2
        return (
                h00 * self.values[index]
                + h * h10 * self.derivatives[index]
                + h01 * self.values[index + 1]
                + h * h11 * self.derivatives[index + 1]
        )

    def __call__(self, time):
        """Calculate using Hermite interpolation."""

        points = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(points)):
            raise ValueError("Time must be finite")
        if points.ndim == 0:
            return self._evaluate_scalar(float(points))

        rows = np.stack([self._evaluate_scalar(float(point)) for point in points.flat])
        return rows.reshape(points.shape + (self.values.shape[1],))


# Helper function calculates the right-hand side of the ODE.
def _right_hand_side(problem, time, state):
    derivative = np.asarray(
        problem.function(time, state, problem.parameters), dtype=float
    )
    if derivative.shape != state.shape:
        raise ValueError("Vectors should match.")
    return derivative


def _dopri5_step(
        problem,
        time,
        state,
        step,
        first_stage,
):
    """Performs a single DOPRI5-step
    Calculates up to k7, then returns the fifth order and fourth order stages.
    It also returns the last stage, which is needed for the FSAL method.
    """

    k1 = _right_hand_side(problem, time, state) if first_stage is None \
        else np.asarray(first_stage, dtype=float)
    k2 = _right_hand_side(
        problem,
        time + step / 5.0,
        state + step * (k1 / 5.0),
    )
    k3 = _right_hand_side(
        problem,
        time + 3.0 * step / 10.0,
        state + step * (3.0 * k1 / 40.0 + 9.0 * k2 / 40.0),
    )
    k4 = _right_hand_side(
        problem,
        time + 4.0 * step / 5.0,
        state
        + step
        * (44.0 * k1 / 45.0 - 56.0 * k2 / 15.0 + 32.0 * k3 / 9.0),
    )
    k5 = _right_hand_side(
        problem,
        time + 8.0 * step / 9.0,
        state
        + step
        * (
                19372.0 * k1 / 6561.0
                - 25360.0 * k2 / 2187.0
                + 64448.0 * k3 / 6561.0
                - 212.0 * k4 / 729.0
        ),
    )
    k6 = _right_hand_side(
        problem,
        time + step,
        state
        + step
        * (
                9017.0 * k1 / 3168.0
                - 355.0 * k2 / 33.0
                + 46732.0 * k3 / 5247.0
                + 49.0 * k4 / 176.0
                - 5103.0 * k5 / 18656.0
        ),
    )
    fifth_order = state + step * (
            35.0 * k1 / 384.0
            + 500.0 * k3 / 1113.0
            + 125.0 * k4 / 192.0
            - 2187.0 * k5 / 6784.0
            + 11.0 * k6 / 84.0
    )
    k7 = _right_hand_side(problem, time + step, fifth_order)
    fourth_order = state + step * (
            5179.0 * k1 / 57600.0
            + 7571.0 * k3 / 16695.0
            + 393.0 * k4 / 640.0
            - 92097.0 * k5 / 339200.0
            + 187.0 * k6 / 2100.0
            + k7 / 40.0
    )
    return fifth_order, fourth_order, k7


def solve(problem, method):
    """Solves the Initial problem using DOPRI5"""

    method = DOPRI5() if method is None else method
    start, end = problem.time_span
    interval_length = end - start

    # Sets up the initial time steps
    if method.initial_step is None:
        step = interval_length / 100.0
    else:
        step = method.initial_step
    step = min(step, method.max_step, interval_length)

    # Initial state
    time = start
    state = problem.initial_value.copy()
    derivative = _right_hand_side(problem, time, state)

    times = [time]
    values = [state.copy()]
    derivatives = [derivative.copy()]

    accepted_steps = 0
    rejected_steps = 0
    function_evaluations = 1

    for _ in range(method.max_steps):
        if time >= end:  # Stop when meeting the end of the interval
            break
        step = min(step, end - time)

        # If the step size is too small, we can end up never moving due to loss of floating point precision
        roundoff_floor = 10.0 * np.finfo(float).eps * max(1.0, abs(time))
        effective_minimum = max(method.min_step, roundoff_floor)
        if step < effective_minimum:
            raise RuntimeError("This catches the problem, where the floating point loses precision.")

        # Calculate a single step of DOPRI5
        fifth_order, fourth_order, last_stage = _dopri5_step(
            problem, time, state, step, derivative
        )

        function_evaluations += 6

        # Precision takes into account the fact that omega and theta might have different magnitutes
        scale = method.atol + method.rtol * np.maximum(
            np.abs(state), np.abs(fifth_order)
        )
        scale = np.maximum(scale, np.finfo(float).tiny)  # Get max error
        error = fifth_order - fourth_order
        error_norm = float(np.sqrt(np.mean((error / scale) ** 2)))  # Mean sqrt error

        if error_norm == 0.0:
            # No error, make large step
            factor = method.max_factor
        elif isfinite(error_norm):
            # This adaption corresponds to h_new = h * 0.9q^-1.5, i got it on the web
            factor = method.safety * error_norm ** (-1.0 / 5.0)
            factor = min(method.max_factor, max(method.min_factor, factor))
        else:  # Error like state, reset
            factor = method.min_factor

        if error_norm <= 1.0:
            # Accept when in tolerance range and iterates further
            time += step
            state = fifth_order
            derivative = last_stage  # FSAL, i.e., First Same As Last
            accepted_steps += 1
            times.append(time)
            values.append(state.copy())
            derivatives.append(derivative.copy())
            step = min(method.max_step, step * factor)
        else:
            rejected_steps += 1
            step *= min(1.0, factor)
    else:
        raise RuntimeError("Max steps reached.")

    return ODESolution(
        problem=problem,
        times=np.asarray(times),
        values=np.asarray(values),
        derivatives=np.asarray(derivatives),
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
        function_evaluations=function_evaluations,
    )


@dataclass(frozen=True)
class PendulumParameters:
    """Parameters of the pendulum."""

    gravity: float = 9.81
    length: float = 1.0
    mass: float = 1.0

    def __post_init__(self):
        if not all(
                isfinite(value) and value > 0.0
                for value in (self.gravity, self.length, self.mass)
        ):
            raise ValueError("All parameters must be positive finite numbers.")


def pendulum_rhs(_time, state, parameters):
    """theta' = omega = - (g/l)sin(theta)"""

    theta, omega = state
    return np.array([omega, -(parameters.gravity / parameters.length) * np.sin(theta)])


def pendulum_problem(theta0,  # Initial value of theta
                     omega0,  # Initial value of omega
                     time_span,
                     parameters=None,
                     ):
    """Creates the initial problem for the pendulum."""

    parameters = PendulumParameters() if parameters is None else parameters
    return InitialValueProblem(
        pendulum_rhs,
        np.array([theta0, omega0], dtype=float),
        time_span,
        parameters,
    )


def solve_pendulum(
        theta0,
        omega0,
        end_time,
        *,
        parameters,
        method,  # DOPRI5
):
    """Estimate the pendulum on the interval [0, end_time] """

    problem = pendulum_problem(theta0, omega0, (0.0, end_time), parameters)
    return solve(problem, method)


def pendulum_angle(
        time,
        theta0,
        omega0=0.0,
        *,
        parameters=None,
        method=None,
):
    """Returns the angle of the pendulum at time t. """

    if time == 0.0:
        return float(theta0)
    if not isfinite(time) or time < 0.0:
        raise ValueError("Non negative and finite time")
    return float(
        solve_pendulum(
            theta0,
            omega0,
            time,
            parameters=parameters,
            method=method,
        ).values[-1, 0]
    )


def harmonic_angle(
        time,
        theta0,
        omega0=0.0,
        *,
        parameters=None,
):
    """Harmonic angle can be directly calculated
        as theta_h(t) = theta_0 cos(sqrt(g/l)t) + omega_0 sin(sqrt(g/l)t)
    """

    parameters = PendulumParameters() if parameters is None else parameters
    frequency = sqrt(parameters.gravity / parameters.length)
    points = np.asarray(time, dtype=float)
    values = theta0 * np.cos(frequency * points) + (
            omega0 / frequency
    ) * np.sin(frequency * points)
    if points.ndim == 0:
        return float(values)
    return values


def harmonic_period(parameters):
    """Harmonic period 2*pi*sqrt(l/g)"""

    parameters = PendulumParameters() if parameters is None else parameters
    return 2.0 * pi * sqrt(parameters.length / parameters.gravity)


def pendulum_energy(
        theta,
        omega,
        parameters,
):
    """Calculates the energy of the pendulum."""

    parameters = PendulumParameters() if parameters is None else parameters
    theta_values = np.asarray(theta, dtype=float)
    omega_values = np.asarray(omega, dtype=float)
    energy = (
            0.5 * parameters.mass * parameters.length ** 2 * omega_values ** 2
            + parameters.mass
            * parameters.gravity
            * parameters.length
            * (1.0 - np.cos(theta_values))
    )
    if energy.ndim == 0:  # If input is a scalar
        return float(energy)
    return energy


# Adapted from vaje 16
class ZeroCrossing:
    """Ničla izbrane komponente rešitve in smer prehoda skozi ničlo."""

    def __init__(self, time, state, direction):
        state = np.array(state, dtype=float, copy=True)
        # -1 if theta decressing through zero, +1 if increasing, 0 if the derivative is zero
        if direction not in (-1, 0, 1):
            raise ValueError("Direction must be -1, 0 or 1.")
        self.time = time
        self.state = state
        self.direction = direction


def _state_from_step_start(solution, index, time):
    """Calculates the state at a time inside one accepted DOPRI5 interval """
    # The left endpoint t_i
    left_time = solution.times[index]
    if time == left_time:
        return solution.values[index].copy()
    # Perform one DOPRI5 step
    state, _, _ = _dopri5_step(
        solution.problem,
        left_time,
        solution.values[index],
        time - left_time,
        solution.derivatives[index],
    )
    return state


def _refine_zero(
        solution,
        index,
        component,
        *,
        rtol,
        atol,
):
    """Refines the zero by bisection because we are working inside a bracket, only need sign change"""
    left = float(solution.times[index])
    right = float(solution.times[index + 1])
    left_value = float(solution.values[index, component])
    right_value = float(solution.values[index + 1, component])
    if left_value == 0.0:
        return left, solution.values[index].copy()
    if right_value == 0.0:
        return right, solution.values[index + 1].copy()
    if left_value * right_value > 0.0:
        raise ValueError("Check we are really clamping around zero")

    state = solution.values[index].copy()
    for _ in range(100):
        middle = left + (right - left) / 2.0
        state = _state_from_step_start(solution, index, middle)
        middle_value = float(state[component])
        # Tolerance
        if middle_value == 0.0 or right - left <= max(atol, rtol * max(1.0, abs(middle))):
            return middle, state
        if left_value * middle_value <= 0.0:
            right = middle
        else:
            left = middle
            left_value = middle_value
    raise RuntimeError("Could not find zero crossing.")


def find_zero_crossings(
        solution,
        component,
        *,
        direction=0,
        rtol=1e-13,
        atol=1e-14,
):
    """Finds zero crossings by paying attention to the sign of the derivative."""

    crossings = []
    # FOr every accepted interval read the left and right values
    for index in range(solution.times.size - 1):
        left_value = float(solution.values[index, component])
        right_value = float(solution.values[index + 1, component])

        if index == 0 and left_value == 0.0:
            # Initial state could be 0, which we dont count
            continue
        if left_value * right_value > 0.0:
            # Skip interval if both values have the same sign
            continue
        if left_value == 0.0 and crossings and (
                abs(solution.times[index] - crossings[-1].time)
                <= max(atol, rtol * max(1.0, abs(solution.times[index])))
        ):
            continue

        time, state = _refine_zero(solution, index, component, rtol=rtol, atol=atol)
        # Derivate of selected component
        derivative = _right_hand_side(solution.problem, time, state)[component]
        # Crossing direction
        crossing_direction = int(np.sign(derivative))

        # Remove duplicate events
        if direction == 0 or crossing_direction == direction:
            if not crossings or abs(time - crossings[-1].time) > max(atol, rtol * max(1.0, abs(time))):
                crossings.append(ZeroCrossing(time, state, crossing_direction))
    return crossings


@dataclass(frozen=True)
class PeriodResult:
    """Stores all the information for our period counter used to plot"""

    period: float
    crossings: tuple[ZeroCrossing, ZeroCrossing]
    solution: ODESolution


def pendulum_period(
        theta0,
        omega0=0.0,
        *,
        parameters=None,
        method=None,
        max_doublings=7,
):
    """Perform energy validation + DOPRI5 + zero crossing for periods"""

    parameters = PendulumParameters() if parameters is None else parameters

    linear_period = harmonic_period(parameters)
    if method is None:
        method = DOPRI5(
            rtol=2e-12,
            atol=2e-14,
            initial_step=linear_period / 100.0,
            max_step=linear_period / 20.0,
        )

    # Select initial horizon
    horizon = 2.0 * linear_period
    for _ in range(max_doublings + 1):
        solution = solve_pendulum(
            theta0,
            omega0,
            horizon,
            parameters=parameters,
            method=method,
        )
        crossings = find_zero_crossings(solution, 0)
        # Group crossing by direction
        for crossing_direction in (-1, 1):
            # Seelct matching crossings
            same_direction = [
                crossing
                for crossing in crossings
                if crossing.direction == crossing_direction
            ]

            if len(same_direction) >= 2:
                # second time - first time
                first, second = same_direction[:2]
                return PeriodResult(second.time - first.time, (first, second), solution)
        horizon *= 2.0

    raise RuntimeError("Did not find two periods in the allowed range.")
