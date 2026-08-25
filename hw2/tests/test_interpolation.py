from math import exp

import numpy as np
import pytest

from interpolation import (
    BarycentricInterpolator,
    chebyshev_nodes,
)


def test_nodes_are_mapped_affinely_to_the_requested_interval():
    nodes = chebyshev_nodes(2, 2.0, 6.0)

    np.testing.assert_allclose(nodes, [6.0, 4.0, 2.0], atol=1e-15)


def test_interpolator_reproduces_a_polynomial_and_preserves_shape():
    polynomial = lambda x: x ** 4 - 2.0 * x ** 2 + 3.0 * x - 1.0
    interpolant = BarycentricInterpolator.from_function(polynomial, -2.0, 3.0, 4)
    points = np.array([[-2.0, -0.7, 0.25], [1.1, 2.4, 3.0]])

    values = interpolant(points)

    assert isinstance(values, np.ndarray)
    assert values.shape == points.shape
    np.testing.assert_allclose(values, polynomial(points), atol=2e-13)
    assert interpolant(0.25) == pytest.approx(polynomial(0.25), abs=2e-13)


def test_function_values_are_returned_exactly_at_nodes():
    interpolant = BarycentricInterpolator.from_function(exp, -1.0, 1.0, 8)

    np.testing.assert_array_equal(interpolant(interpolant.nodes), interpolant.values)
