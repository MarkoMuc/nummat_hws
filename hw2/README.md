# Homework 2: Barycentric Interpolation and Gauss–Legendre Quadrature

**Author:** Marko Zupančič Muc

Implement polynomial interpolation using the barycentric Lagrange formula and Chebyshev--Lobatto nodes.
The nodes can be mapped from $[-1,1]$ to any interval $[a,b]$.

The implementation is used to interpolate three example functions and examine
how their maximum interpolation errors change as the polynomial degree
increases.

The second part implements the composite two-point Gauss–Legendre rule.
It also estimates the integration error by doubling the number of subintervals until the requested tolerance is reached.

## Implementation

- `src/interpolation.py` contains the implementation of the barycentric interpolation algorithm.
- `src/quadrature.py` contains the implementation of Gauss-Legendre quadrature.
- `scripts/part1_interpolate.py` runs the interpolation experiments, stores the numerical results, 
  and generates the plots used in the report.
- `scripts/part2_quadrature.py` approximates the integral of $\sin(x)/x$ on $[0,5]$ and plots the results.

## Example use

### Interpolation

```python
from math import exp

from src.interpolation import BarycentricInterpolator

interpolant = BarycentricInterpolator.from_function(
    lambda x: exp(-x**2),
    -1.0,
    1.0,
    10,
)

value = interpolant(0.3)
```

### Quadrature

```python
from src.quadrature import CompositeGaussLegendre2, Integral, Interval, integrate

problem = Integral(lambda x: x**3 - 2*x + 1, Interval(-2.0, 3.0))
approx = integrate(problem, CompositeGaussLegendre2(8))
```

## Tests

Run all tests:

```bash
uv run pytest
```

To run only the interpolation tests:

```bash
uv run pytest tests/test_interpolation.py
```

To run only the Gauss–Legendre quadrature tests:

```bash
uv run pytest tests/test_gauss_legendre.py
```

## Generate the results and plots

```bash
uv run python -m scripts.part1_interpolate
uv run python -m scripts.part2_quadrature
```

## Compile the report

```bash
typst compile hw2.typ
```
