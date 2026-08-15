import numpy as np

from sparse_matrix import RedkaMatrika


# Pretty much the implementation from Vaje06
def conj_grad(
        A: RedkaMatrika,
        b,
        *,
        tol=1e-10,
        max_iter=None,
        x0=None,
):
    """Solve the linear system A @ x = b with the conjugate gradient method

    Matrix A is symmetrical and positive definite.
    The stopping condition is ||b - A @ x||_2 < tol``.
    """

    right_side = np.asarray(b, dtype=float)
    n, m = A.shape

    # Some check are not needed, since our RedkaMatrika does similar checks
    # But we still explicitly check here
    # We do not check if A is P.D.
    if n != m or right_side.ndim != 1 or len(right_side) != n:
        raise ValueError("Check dimensions of matrices and vectors.")
    if not np.isfinite(right_side).all():
        raise ValueError("Finite numbers only")
    if tol <= 0:
        raise ValueError("Tolerance parameter must be positive.")

    if max_iter is None:
        max_iter = n
    if max_iter < 0:
        raise ValueError("max_iter must be positive")

    # If no initial approximation is given, just take the right side
    if x0 is None:
        x = right_side.copy()
    else:
        x = np.asarray(x0, dtype=float).copy()
        if x.ndim != 1 or len(x) != n:
            raise ValueError("Dimensions of init. approximation are not compatible.")
        if not np.isfinite(x).all():
            raise ValueError("Finite numbers only")

    residual = right_side - A @ x  # Measures how much current approx. satisfies the system
    direction = residual.copy()  # Move in direction opposite the gradient
    residual_norm_squared = float(residual @ residual)  # Dot product in numpy
    if np.sqrt(residual_norm_squared) < tol:
        return x, 0

    for iteration in range(1, max_iter + 1):  # From 1 since the first iteration was already done
        matrix_direction = A @ direction
        curvature = float(direction @ matrix_direction)
        if curvature <= 0 or not np.isfinite(curvature):
            raise ValueError("Matrix is not P.D.")

        alpha = residual_norm_squared / curvature  # Alpha represents the step size
        x += alpha * direction  # Updating the approximation
        residual -= alpha * matrix_direction  # Negative gradient, we want to descent
        new_residual_norm_squared = float(residual @ residual)

        if np.sqrt(new_residual_norm_squared) < tol:
            return x, iteration

        # Update the search direction so it is conjugate to the previous one.
        direction = residual + (new_residual_norm_squared / residual_norm_squared) * direction
        residual_norm_squared = new_residual_norm_squared

    return x, max_iter
