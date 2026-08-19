from math import exp, sin
from pathlib import Path

import json
import matplotlib.pyplot as plt

from src.interpolation import *

MAX_ERROR = 1e-6

def gaussian(x):
    return exp(-(x ** 2))


def sinc(x: float) -> float:
    if x == 0.0:
        return 1.0

    return sin(x) / x


def f3(x):
    return abs(x ** 2 - 2.0 * x)


def first_passing_degree(function, interval, maximum_degree, grid_size):
    """Search for the first passing polynomial degrees, so the error is < 10^-6"""
    grid = np.linspace(*interval, grid_size)  # Where we will test it out
    history = []
    for degree in range(1, maximum_degree + 1):
        interpolant = BarycentricInterpolator.from_function(function, *interval, degree)  # Creates the interpolant

        error, location = max_error_on_grid(function, interpolant, grid)  # Find the maximum error
        history.append((degree, error, location))  # Store degree, max error and location for graphing
        if error < MAX_ERROR:  # We've found the first passing degree
            return interpolant, history
    raise RuntimeError("No passing degree found on the provided interval.")


def interpolation_experiments():
    experiments = [
        ("exp(-x^2)", gaussian, (-1.0, 1.0), 20, 20_001),
        ("sin(x)/x", sinc, (0.0, 10.0), 25, 40_001),
    ]
    results = {}
    interpolants = {}

    for name, function, interval, maximum_degree, grid_size in experiments:
        # Find the first passing degree (depending on the max error)
        interpolant, history = first_passing_degree(function, interval, maximum_degree, grid_size)
        degree, error, location = history[-1]  # Last entry is the first passing degree
        results[name] = {
            "interval": interval,
            "degree": degree,
            "observed_max_error": error,
            "error_location": location,
            "control_points": grid_size,
            "previous_degree_error": history[-2][1],
        }  # Save the results
        interpolants[name] = interpolant

    # The degree for the polynomial of the third function is > 1million,
    # Since the derivative from the left differs from the right one, corner in x=2.
    # That's why we process it separately
    # Interval is [1.0, 3.0]
    nonsmooth_degrees = [20, 50, 100, 200, 500, 1_000, 2_000, 5_000]
    convergence = []
    global_grid = np.linspace(1.0, 3.0, 20_001)
    for degree in nonsmooth_degrees:
        interpolant = BarycentricInterpolator.from_function(
            f3, 1.0, 3.0, degree
        ) # Construct the interpolant
        local_grid = 2.0 + np.linspace(-2.5, 2.5, 401) / degree # Error grows around x=2, don't want to miss it
        control_grid = np.unique(np.concatenate((global_grid, local_grid))) # Local + Global grid
        error, location = max_error_on_grid(
            f3, interpolant, control_grid
        ) # Find max error
        convergence.append(
            {"degree": degree, "error": error, "location": location}
        ) # Add to convergence

    final_degree = 1_200_000
    comparison_degree = 1_000_000 # Test a large degree
    comparison_interpolant = BarycentricInterpolator.from_function(f3, 1.0, 3.0, comparison_degree)
    comparison_grid = np.unique(
        np.concatenate(
            (
                np.linspace(1.0, 3.0, 81),
                2.0 + np.linspace(-2.5, 2.5, 401) / comparison_degree,
            )
        )
    )
    comparison_error, comparison_location = max_error_on_grid(
        f3, comparison_interpolant, comparison_grid
    ) # Max error with degree 1 million
    final_interpolant = BarycentricInterpolator.from_function(
        f3, 1.0, 3.0, final_degree
    ) # Max error with the final degree in 1.2 million range, space is interval we take is smaller
    local_grid = 2.0 + np.linspace(-2.5, 2.5, 401) / final_degree
    verification_grid = np.unique(np.concatenate((np.linspace(1.0, 3.0, 81), local_grid)))
    # Compare solution with 1.2 million and the real function
    final_error, final_location = max_error_on_grid(f3, final_interpolant, verification_grid)

    results["abs(x^2-2x)"] = {
        "interval": (1.0, 3.0),
        "degree": final_degree,
        "observed_max_error": final_error,
        "error_location": final_location,
        "control_points": len(verification_grid),
        "comparison": {
            "degree": comparison_degree,
            "observed_max_error": comparison_error,
            "error_location": comparison_location,
        },
        "moderate_degree_convergence": convergence,
    }

    interpolants["abs(x^2-2x)"] = final_interpolant
    return results, interpolants

# This plotting function was written by an LLM and then hand-checked by me
def draw_interpolation_figures(results, interpolants, image_directory):
    comparison_cases = [
        ("exp(-x^2)", gaussian, (-1.0, 1.0), 31),
        ("sin(x)/x", sinc, (0.0, 10.0), 41),
        ("abs(x^2-2x)", f3, (1.0, 3.0), 41),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.6), constrained_layout=True) # 3 graphs
    for axis, (name, function, interval, sample_count) in zip(
        axes, comparison_cases
    ):
        curve_grid = np.linspace(*interval, 1_001) #1001 spaced points
        exact = np.fromiter(
            (function(float(x)) for x in curve_grid),
            dtype=float,
            count=curve_grid.size,
        ) # Exact function evaluated
        sample_points = np.linspace(*interval, sample_count)
        approximation = interpolants[name](sample_points) # Evaluate multiple points with the interpolant

        axis.plot(
            curve_grid,
            exact,
            color="#2e77bb",
            linewidth=1.8,
            label="Function $f(x)$",
        ) # Draw the real function curve
        axis.scatter(
            sample_points,
            approximation,
            color="#d1495b",
            edgecolor="white",
            linewidth=0.45,
            s=25,
            zorder=3,
            label="Interpolant $p_n(x)$",
        ) # Draw the interpolant as nodes
        axis.set(
            title=f"{name}, n={results[name]['degree']}",
            xlabel="x",
            ylabel="Value",
        ) # Add in legend and the degree of the polynomial
        axis.grid(alpha=0.25)

    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    figure.savefig(image_directory / "hw2_interpolants.svg", bbox_inches="tight")
    plt.close(figure) # Save

    # Plot the absolute errors
    figure, axes = plt.subplots(1, 3, figsize=(12.0, 3.7), constrained_layout=True)
    cases = [
        ("exp(-x^2)", gaussian, (-1.0, 1.0)),
        ("sin(x)/x", sinc, (0.0, 10.0)),
    ]
    for axis, (name, function, interval) in zip(axes, cases):
        grid = np.linspace(*interval, 2_001) # Create the grid
        exact = np.fromiter((function(float(x)) for x in grid), float) # Evaluate exact function
        error = np.abs(interpolants[name](grid) - exact) # Evaluate the error (p(x) - f(x))
        axis.semilogy(grid, np.maximum(error, 1e-18), color="#2e77bb")
        axis.axhline(1e-6, color="#d1495b", linestyle="--", linewidth=1.1)
        axis.set(
            title=f"{name}, n={results[name]['degree']}",
            xlabel="x",
            ylabel="Absolute error",
        )
        axis.grid(alpha=0.25)

    # Absolute error of the third function
    name = "abs(x^2-2x)"
    degree = results[name]["degree"]
    scaled_grid = np.linspace(-2.5, 2.5, 801)
    corner_grid = 2.0 + scaled_grid / degree # Scale the coordinates, so we can capture what happens around x=2
    corner_exact = np.fromiter(
        (f3(float(x)) for x in corner_grid),
        dtype=float,
        count=corner_grid.size,) # Evaluate the exact values
    corner_error = np.abs(interpolants[name](corner_grid) - corner_exact) # Evaluate the error (p(x) - f(x))
    axes[2].semilogy(
        scaled_grid,
        np.maximum(corner_error, 1e-18),
        color="#2e77bb",
    )
    axes[2].axhline(1e-6, color="#d1495b", linestyle="--", linewidth=1.1)
    axes[2].set(
        title=f"{name}, n={degree} (locally)",
        xlabel=r"$n(x-2)$",
        ylabel="Absolute error",
    )
    axes[2].grid(alpha=0.25)
    figure.savefig(image_directory / "hw2_interpolation_errors.svg", bbox_inches="tight")
    plt.close(figure)

    # Plot convergence for third function
    convergence = results["abs(x^2-2x)"]["moderate_degree_convergence"]
    degrees = np.array([row["degree"] for row in convergence], dtype=float) # Get the full history of the degree
    errors = np.array([row["error"] for row in convergence]) # History of the errors
    final_degree = results["abs(x^2-2x)"]["degree"] # End degree
    final_error = results["abs(x^2-2x)"]["observed_max_error"] # Max error found at the end
    figure, axis = plt.subplots(figsize=(6.0, 3.8), constrained_layout=True)
    axis.loglog(degrees, errors, "o-", color="#2e77bb", label="Measured error")
    axis.loglog(
        [degrees[-1], final_degree],
        [errors[-1], errors[-1] * degrees[-1] / final_degree],
        "--",
        color="#758195",
        label=r"Extrapolation of $C/n$",
    )
    axis.loglog(
        [final_degree], [final_error], "s", color="#d1495b", label="Last check"
    )
    axis.axhline(1e-6, color="#d1495b", linestyle=":", linewidth=1.1)
    axis.set(xlabel="Degree $n$", ylabel="Max error")
    axis.grid(alpha=0.25, which="both")
    axis.legend(frameon=False)
    figure.savefig(image_directory / "hw2_convergence.svg", bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    image_directory = Path("img")
    result_directory = Path("results")
    image_directory.mkdir(exist_ok=True)
    result_directory.mkdir(exist_ok=True)

    interpolation_results, interpolants = interpolation_experiments()  # Perform the interpolation experiments
    results = {"interpolation": interpolation_results, }

    draw_interpolation_figures(interpolation_results, interpolants, image_directory)  # Draw the figures for interp.

    output_path = result_directory / "hw2_part1.json"
    output_path.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))