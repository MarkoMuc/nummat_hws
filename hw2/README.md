# Homework 2: Barycentric Interpolation

**Author:** Marko Zupančič Muc

Implement polynomial interpolation using the barycentric Lagrange formula and Chebyshev--Lobatto nodes.
The nodes can be mapped from $[-1,1]$ to any interval $[a,b]$.

The implementation is used to interpolate three example functions and examine
how their maximum interpolation errors change as the polynomial degree
increases.

## Implementation

- `src/interpolation.py` contains the implementation of the barycentric interpolation algorithm.
- `scripts/part1_interpolate.py` runs the interpolation experiments, stores the numerical results, 
  and generates the plots used in the report.

**Example use**

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

## Tests

Run all tests:

```bash
uv run pytest
```

To run only the interpolation tests:

```bash
uv run pytest tests/test_interpolation.py
```

## Generate the results and plots

```bash
uv run python -m scripts.part1_interpolate
```

## Compile the report

```bash
typst compile hw2.typ
```
