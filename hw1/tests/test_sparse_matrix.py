import numpy as np
import pytest

from sparse_matrix import RedkaMatrika


def test_storage_indexing_and_shape():
    dense = np.array([[4.0, 0.0, -1.0], [0.0, 2.0, 0.0], [-1.0, 0.0, 3.0]])
    matrix = RedkaMatrika.from_dense(dense)

    assert matrix.shape == (3, 3)
    assert matrix.nnz == 5
    assert matrix.V == [[4.0, -1.0], [2.0], [-1.0, 3.0]]
    assert matrix.I == [[0, 2], [1], [0, 2]]
    assert matrix[0, 1] == 0.0
    assert matrix[2, 0] == -1.0
    np.testing.assert_array_equal(matrix.to_dense(), dense)


def test_setting_inserts_updates_and_removes_entry():
    matrix = RedkaMatrika.zeros(3)

    matrix[1, 2] = 5
    matrix[1, 0] = 2
    matrix[1, 2] = 7
    assert matrix.I[1] == [0, 2]
    assert matrix.V[1] == [2.0, 7.0]

    matrix[1, 0] = 0
    assert matrix.I[1] == [2]
    assert matrix.V[1] == [7.0]
    assert matrix.nnz == 1


def test_sparse_matrix_vector_product_matches_hand_calculation():
    matrix = RedkaMatrika(
        [[4, -1], [-1, 4, -1], [-1, 3]],
        [[0, 1], [0, 1, 2], [1, 2]],
    )

    np.testing.assert_allclose(matrix @ np.array([1.0, 2.0, 3.0]), [2.0, 4.0, 7.0])


def test_invalid_sparse_matrix_operations_are_rejected():
    with pytest.raises(ValueError):
        RedkaMatrika([[1]], [[], []])
    with pytest.raises(ValueError):
        RedkaMatrika([[1, 2]], [[0, 0]])
    with pytest.raises(ValueError):
        RedkaMatrika.from_dense(np.ones((2, 3)))

    matrix = RedkaMatrika.zeros(2)
    with pytest.raises(IndexError):
        _ = matrix[2, 0]
    with pytest.raises(ValueError):
        _ = matrix @ np.ones(3)
