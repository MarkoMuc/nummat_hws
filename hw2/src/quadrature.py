from dataclasses import dataclass
from math import sqrt, isfinite, sin
from typing import Callable

# A lot of the implementations here are similar or identical to @Vaje13 for chapter 13 of the book


@dataclass(frozen=True)
class Interval:
    """"Interval structure [a,b]"""
    a: float
    b: float
    def __post_init__(self):
        if not isfinite(self.a) or not isfinite(self.b):
            raise ValueError("The end points must be finite numbers.")

    @property
    def length(self):
        return self.b - self.a

@dataclass(frozen=True)
class Integral:
    """Definite integral of the function on an interval [a,b]."""
    function: Callable[[float], float] # This type defines a scalar function
    interval: Interval


@dataclass(frozen=True)
class CompositeGaussLegendre2:
    """Composite two point Gauss-Legendre"""

    N: int # n

    def __post_init__(self):
        if self.N < 1:
            raise ValueError("Number of subintervals needs to be >0.")

    @property
    def evaluations(self):
        """Number of func. evaluations"""

        return 2 * self.N


@dataclass(frozen=True)
class StepDoublingResult:
    """Quadrature result using the step doubling estimate"""

    value: float #Q_2n
    richardson_value: float # I_n = Q_2n + (Q_2n - Q_n)/15
    error_estimate: float # |(Q_2n - q_n) / 15|
    N: int
    evaluations: int # Final cost of the quadrature
    total_evaluations: int # All earlier trials


def integrate(integral, rule):
    """Approximation of the integral using the two point Gauss-Legendre rule

        Q_N = h/2 sum^(n-1)_0 (f (m_i - h/(2sqrt(3))) + f (m_i + h/(2sqrt(3))))
        where
        h = (b - a) / N
        m_i = a + (i + 1/2) * h
    """

    a = float(integral.interval.a)
    h = integral.interval.length / rule.N # Interval length of the subintervals

    offset = h / (2.0 * sqrt(3.0)) # Mapping from [-1,1] to [a,b] so h/(2 sqrt(3))
    weighted_sum = 0.0 # Sum of the subintervals
    for i in range(rule.N): # Repeat over subintervals
        m_i = a + (i + 0.5) * h
        fst_value = float(integral.function(m_i - offset)) # f(m_i - h/(2sqrt(3))
        snd_value = float(integral.function(m_i + offset)) # f(m_i + h/(2sqrt(3))
        weighted_sum += fst_value + snd_value
    return h * weighted_sum / 2.0 # h/2 * weighted_sum(subintervals)


def estimate_by_step_doubling(
    integral,
    *,
    rtol= 1e-10,
    atol= 0.0,
    start_N = 1,
    max_N = 1 << 20):
    """Doubles the number of intervals until we reach the error estimate

    The error is estimated wtih |Q_2n - Q_n|/ 15.
    """

    if rtol == 0.0 and atol == 0.0:
        raise ValueError("At least one of the tolerances must be nonzero.")
    if start_N < 1:
        raise ValueError("Starting N must be a positive integer.")
    if max_N < 2 * start_N:
        raise ValueError("Maximum N must be at least 2N.")

    # Start with initial estimation
    curr_N = start_N
    coarse = integrate(integral, CompositeGaussLegendre2(curr_N)) # This is our Q_n
    total_evaluations = 2 * curr_N # Calculates the cost of the initial Q_n

    while 2 * curr_N <= max_N:
        curr_N *= 2 # n -> 2n
        fine = integrate(integral, CompositeGaussLegendre2(curr_N)) # This is Q_2n
        total_evaluations += 2 * curr_N # Number of evaluations for Q_2n

        # Check the tolerance
        difference = fine - coarse
        error_estimate = abs(difference) / 15.0 # |Q_2n - Q_n| / 15
        tolerance = max(atol, rtol * abs(fine))
        if error_estimate <= tolerance:
            return StepDoublingResult(
                value=fine,
                richardson_value=fine + difference / 15.0,
                error_estimate=error_estimate,
                N=curr_N,
                evaluations=2 * curr_N,
                total_evaluations=total_evaluations,
            )
        coarse = fine

    raise RuntimeError("We did not reach the desired tolerance.")
