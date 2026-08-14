import numpy as np
import pytest

from conjugate_gradient import conj_grad
from sparse_matrix import RedkaMatrika


def test_conjugate_gradient_recovers_known_solution():
    dense = np.array(
        [
            [4.0, -1.0, 0.0, 0.0],
            [-1.0, 4.0, -1.0, 0.0],
            [0.0, -1.0, 4.0, -1.0],
            [0.0, 0.0, -1.0, 3.0],
        ]
    ) #symmetric PD matrix
    expected = np.array([1.0, -2.0, 3.0, 0.5])
    matrix = RedkaMatrika.from_dense(dense)
    rhs = dense @ expected

    solution, iterations = conj_grad(matrix, rhs)

    np.testing.assert_allclose(solution, expected, atol=1e-11)
    assert 1 <= iterations <= len(rhs)
    assert np.linalg.norm(matrix @ solution - rhs) < 1e-10


def test_zero_residual_finishes_without_iteration():
    matrix = RedkaMatrika.from_dense(np.eye(3))
    rhs = np.array([1.0, 2.0, 3.0])

    solution, iterations = conj_grad(matrix, rhs)

    np.testing.assert_array_equal(solution, rhs)
    assert iterations == 0


def test_initial_approximation_and_iteration_limit_are_respected():
    matrix = RedkaMatrika.from_dense([[2.0, -1.0], [-1.0, 2.0]])
    rhs = np.array([1.0, 0.0])

    solution, iterations = conj_grad(matrix, rhs, x0=np.zeros(2), max_iter=1)

    assert iterations == 1
    assert np.linalg.norm(matrix @ solution - rhs) > 1e-10


def test_non_positive_definite_matrix_is_detected():
    matrix = RedkaMatrika.from_dense([[1.0, 0.0], [0.0, -1.0]])

    with pytest.raises(ValueError, match=r"P\.D\."):
        conj_grad(matrix, [1.0, 1.0], x0=[0.0, 0.0])


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"tol": 0.0}, "Tolerance"),
        ({"max_iter": -1}, "max_iter"),
        ({"x0": [0.0]}, "Dimensions"),
    ],
)
def test_invalid_solver_arguments(kwargs, message):
    matrix = RedkaMatrika.from_dense(np.eye(2))
    with pytest.raises(ValueError, match=message):
        conj_grad(matrix, [1.0, 2.0], **kwargs)
