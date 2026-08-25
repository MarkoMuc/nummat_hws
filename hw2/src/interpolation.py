import numpy as np


def chebyshev_nodes(n, a=-1.0, b=1.0):
    """Returns chebyshev nodes on the interval [a,b] using the formula x_k=cos((2k-1)/2n pi) for k = 0,1,..., n-1
    The nodes are mapped from [-1,1] to [a,b]
    """

    a, b = float(a), float(b)

    k = np.arange(n + 1, dtype=float)  # Indices
    cnodes = np.cos(np.pi * k / n)

    # This is as described in the notes, but I think its a mistake
    # k = np.arange(n, dtype=float) # Indices
    # cnodes = np.cos((2.0 * k - 1.0) * np.pi / (2.0 * n)) # Chebyshev nodes generated with the formula

    # Perform the mapping [-1,1] -> [a, b]
    scale = (b - a) / 2.0  # 2.0 = -1.0 - 1.0, the target interval.
    return (a + b) / 2.0 + scale * cnodes  # Simplified form of: scale * (cnodes - (-1.0)) + a


def chebyshev_weights(n):
    """Returns chebyshev weights using the formula lambda_k=(-1)^k case(1 if 0 < i < n, 1/2 if i=0, n else"""

    weights = np.ones(n + 1, dtype=float)
    weights[1::2] = -1.0
    weights[[0, -1]] *= 0.5

    # This is as described in the notes, but I think its a mistake
    # indices = np.arange(n, dtype=int)  # k = 0, 1, ..., n - 1
    # factors = np.where(
    #    (0 < indices) & (indices < n),
    #    1.0,
    #    np.where(indices == 0, 0.5, float(n))) # Calculate the cases
    # weights = (-1.0) ** indices * factors

    return weights


class BarycentricInterpolator:
    """Barycentric interpolation on the interval [a,b]

    Represents the interpolation polynomial with chebyshev nodes and weight,
    and function values in the corresponding nodes.
    """

    def __init__(self, nodes, values, weights, interval):
        self.nodes = np.asarray(nodes, dtype=float)
        self.values = np.asarray(values, dtype=float)
        self.weights = np.asarray(weights, dtype=float)
        self.interval = interval

        if nodes.ndim != 1 or values.ndim != 1 or weights.ndim != 1:
            raise ValueError("Must be 1D.")
        if not (len(nodes) == len(values) == len(weights)):
            raise ValueError("Must be of same length.")
        if len(nodes) < 2:
            raise ValueError("At least two nodes are required.")
        if np.any(weights == 0.0):
            raise ValueError("Weights must be nonzero.")
        if len(np.unique(nodes)) != len(nodes):
            raise ValueError("Interpolation nodes must be unique.")

    @property
    def degree(self):
        """Degree of the polynomial"""

        return len(self.nodes) - 1

    @classmethod
    def from_function(cls, function, a, b, n):
        """Static function for creating an instance of the class"""

        # Generate the chebyshev nodes (and map them onto to the [a,b] interval)
        nodes = chebyshev_nodes(n, a, b)

        # Evaluate the function in the chebyshev nodes
        values = np.fromiter(
            (float(function(float(node))) for node in nodes),
            dtype=float,
            count=n + 1,
        )
        # Generate the chebyshev weights
        weights = chebyshev_weights(n)
        return cls(nodes, values, weights, (float(a), float(b)))

    def _evaluate_scalar(self, point):
        """Returns the value of the barycentric Lagrande interpolation using the provided formula"""
        diff = point - self.nodes  # x - x_j
        matches = np.flatnonzero(diff == 0.0)  # Returns indices that are zero
        if matches.size:  # Evaluates l(x) = f(x_j) if x = x_j
            return float(self.values[matches[0]])

        # Evaluates l(x) = sum(f(x_j)weight_j/diff)/sum(weight/diff)
        quot = self.weights / diff
        return float(np.dot(quot, self.values) / np.sum(quot))

    def __call__(self, x):
        """Interpolates over the provided points"""

        points = np.asarray(x, dtype=float)

        if points.ndim == 0:  # In a single point
            return self._evaluate_scalar(float(points))

        # Interpolate over all points
        result = np.fromiter(
            (self._evaluate_scalar(float(point)) for point in points.flat),
            dtype=float,
            count=points.size,
        )
        # TODO: Why is reshape called here?
        return result.reshape(points.shape)


def max_error_on_grid(function, interpolant, points):
    """Finds the max absolute error, used for finding the first passing degree"""

    grid = np.asarray(points, dtype=float)  # Points over which we interpolate
    if grid.ndim != 1 or grid.size == 0:
        raise ValueError("The point grid should be nonzero and one dimensional.")

    actual = np.fromiter(
        (float(function(float(point))) for point in grid),  # Direct values from the function we are interpolating
        dtype=float,
        count=grid.size,
    )
    errors = np.abs(np.asarray(interpolant(grid)) - actual)  # Calculates the errors
    index = int(np.argmax(errors))  # Finds the index of max error
    return float(errors[index]), float(grid[index])  # Returns the error and the corresponding point
