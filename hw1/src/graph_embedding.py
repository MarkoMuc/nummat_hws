import numpy as np

from .conjugate_gradient import conj_grad
from .sparse_matrix import RedkaMatrika


# The Graph embedding functions are from Vaje06 + The new sparse matrix and conjugate gradient implementations

def _add_edge(graph, first, second):
    graph[first].add(second)
    graph[second].add(first)


def circular_ladder(n):
    """ Circular ladder graph with 2n nodes"""

    graph = [set() for _ in range(2 * n)]  # Adjacency list
    for vertex in range(n):
        next_vertex = (vertex + 1) % n

        _add_edge(graph, vertex, next_vertex)  # Outer circle edge
        _add_edge(graph, vertex + n, next_vertex + n)  # Inner circle edge
        _add_edge(graph, vertex, vertex + n)  # The ladder edges
    return graph


def physical_system(graph, free_vertices):
    """ Creates a P.D. matrix of the physical system

    Note that the real system matrix A would be negative definite, so we are returning -A!
    We do this by not negating the deg beforehand and insertin -1 if they are neighbors.
    """

    free = list(free_vertices)

    position = {vertex: index for index, vertex in enumerate(free)}
    matrix = RedkaMatrika.zeros(len(free))
    # The starting matrix is already the zero matrix
    # This way we already have the a_ij that equal zero (not neighbors).
    for row, vertex in enumerate(free):
        matrix[row, row] = len(graph[vertex])  # The diagonal has degrees
        for neighbor in graph[vertex]:
            if neighbor in position:
                matrix[row, position[neighbor]] = -1.0  # If neighbors add -1 otherwise 0
    return matrix


def right_hand_side(graph, free_vertices, coordinates):
    """Right hand side of the system of coordinates of the fixed vertices"""

    free = list(free_vertices)
    coordinate = np.asarray(coordinates, dtype=float)
    if coordinate.ndim != 1 or len(coordinate) != len(graph):
        raise ValueError("Single coordinate for each vertex.")

    free_set = set(free)
    return np.array(
        [
            sum(coordinate[neighbor] for neighbor in graph[vertex] if neighbor not in free_set)  # Fixed neighbors only
            for vertex in free  # For each free vertex
        ],
        dtype=float,
    )


def embed_graph(graph, fixed_vertices, points, *, tol=1e-10, ):
    """ Embed the graph """

    coordinates = np.asarray(points, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != len(graph):
        raise ValueError("Coordinate of the graph must match the nodes of the graph.")

    fixed = list(fixed_vertices)
    if not fixed:
        raise ValueError("At least one node needs to be fixed")

    fixed_set = set(fixed)
    free = [vertex for vertex in range(len(graph)) if vertex not in fixed_set]  # Get the free set
    result = coordinates.copy()
    if not free:
        return result, [0] * coordinates.shape[0]

    matrix = physical_system(graph, free)  # Get matrix A which is the physical system
    iterations = []
    for dimension in range(coordinates.shape[0]):
        rhs = right_hand_side(graph, free, coordinates[dimension])
        solution, count = conj_grad(matrix, rhs, tol=tol)
        result[dimension, free] = solution
        iterations.append(count)
    return result, iterations
