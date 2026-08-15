from pathlib import Path

import matplotlib.pyplot as plt

from src.graph_embedding import *


def draw_graph(axis, graph, points, fixed, title):
    """Draws the graph
        **This plotting function was generated with an LLM and checked by hand**
    """

    for first, neighbors in enumerate(graph):
        for second in neighbors:
            # Draws the edges
            if first < second:
                axis.plot(
                    points[0, [first, second]],  # x coordinate
                    points[1, [first, second]],  # y coordinate
                    color="#758195",
                    linewidth=1.4,
                    zorder=1,
                )
    vertices = np.arange(len(graph))  # Enumerate all the vertices
    fixed_mask = np.isin(vertices, list(fixed))  # Get the fixed vertices as a boolean mask
    axis.scatter(
        points[0, ~fixed_mask],  # Draw the free vertices
        points[1, ~fixed_mask],
        color="#2e77bb",
        s=42,
        label="Free vertices",
        zorder=2,
    )
    axis.scatter(
        points[0, fixed_mask],  # Draw the fixed vertices
        points[1, fixed_mask],
        color="#d1495b",
        marker="s",
        s=46,
        label="Fixed vertices",
        zorder=3,
    )
    axis.set(title=title, aspect="equal")
    axis.axis("off")


if __name__ == "__main__":
    # Ladder graph with 8 edges and 16 vertices
    # 2 cycles, from 0...7 and 8...15
    n = 8
    graph = circular_ladder(n)
    fixed = range(n)  # Fix first 8 vertices
    angles = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)  # Place the fixed vertices on a circle
    initial = np.column_stack((
        np.vstack((np.cos(angles), np.sin(angles))),
        np.zeros((2, n)),))  # Initial coordinate of fixed and free vertices (matrix= fixed | free)
    embedded, iterations = embed_graph(graph, fixed, initial)  # Compute the embedding

    # Here we reconstruct the system matrix A, so we analyse the number of zero elements
    # We also check if the solution is correct, by using the fact that the residual = Ax - b
    # If the embedding is correct, its norm should be close to zero (better the approximation)
    matrix = physical_system(graph, list(range(n, 2 * n)))
    residual_norms = []
    for dimension in range(2):
        rhs = right_hand_side(graph, list(range(n, 2 * n)), embedded[dimension])
        residual_norms.append(
            float(np.linalg.norm(matrix @ embedded[dimension, n:] - rhs))
        )

    image_directory = Path("img")
    image_directory.mkdir(exist_ok=True)

    figure, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), constrained_layout=True)
    draw_graph(axes[0], graph, initial, fixed, "Starting assignment")
    draw_graph(axes[1], graph, embedded, fixed, "Physical assignment")
    handles, labels = axes[1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
    figure.savefig(image_directory / "circular_ladder_embedding.svg", bbox_inches="tight")
    plt.close(figure)

    print("Graph:", "circular ladder")
    print("Vertices:", len(graph))
    print("Fixed vertices:", n)
    print("System shape:", matrix.shape)
    print("Nonzero elements:", matrix.nnz)
    print("Iterations:", iterations)
    print("Residual norms:", residual_norms)
