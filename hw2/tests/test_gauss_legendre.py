import pytest

from quadrature import (
    CompositeGaussLegendre2,
    Integral,
    Interval,
    integrate,
)


def test_exact_for_deg3_polynomials():
    function = lambda x: x**3 - 2.0 * x + 1.0
    antiderivative = lambda x: x**4 / 4.0 - x**2 + x
    interval = Interval(-2.0, 3.0)

    approximation = integrate(
        Integral(function, interval),
        CompositeGaussLegendre2(1)
    )

    assert approximation == pytest.approx( antiderivative(interval.b) - antiderivative(interval.a), abs=2e-14 )


def test_quartic_error_matches_the_derived_remainder():
    height = 3.0
    exact = height**5 / 5.0
    approximation = integrate(
        Integral(lambda x: x**4, Interval(0.0, height)),
        CompositeGaussLegendre2(1),
    )

    # For f(x)=x^4 je f^(4)=24, thus R=h^5*24/4320=h^5/180.
    assert exact - approximation == pytest.approx(height**5 / 180.0, rel=1e-14)


def test_composite_rule_has_fourth_order_convergence():
    integral = Integral(lambda x: x**4, Interval(0.0, 1.0))
    error_4 = abs(1.0 / 5.0 - integrate(integral, CompositeGaussLegendre2(4)))
    error_8 = abs(1.0 / 5.0 - integrate(integral, CompositeGaussLegendre2(8)))

    assert error_4 / error_8 == pytest.approx(16.0, rel=2e-10)


def test_reversed_interval_changes_the_sign():
    forward = integrate(
        Integral(lambda x: x**2, Interval(-1.0, 2.0)),
        CompositeGaussLegendre2(3),
    )
    backward = integrate(
        Integral(lambda x: x**2, Interval(2.0, -1.0)),
        CompositeGaussLegendre2(3),
    )

    assert backward == pytest.approx(-forward, abs=1e-14)
