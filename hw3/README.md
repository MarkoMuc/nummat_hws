# Homework 3: Mathematical Pendulum

**Author:** Marko Zupančič Muc

Model the mathematical pendulum as a first-order system and solve it with the adaptive DOPRI5 method.

The numerical solution is compared with the harmonic approximation for several initial angles.
The implementation is also used to measure how the pendulum's period changes with its initial energy.

## Implementation

- `src/mpendulum.py` the DOPRI5 solver, the pendulum model, etc.
- `scripts/hw3_mpendulum.py` runs the numerical experiments and generates the plots used in the report.

## Example use

```python
from src.mpendulum import (
    DOPRI5,
    PendulumParameters,
    harmonic_period,
    pendulum_period,
    solve_pendulum,
)

parameters = PendulumParameters(gravity=9.81, length=1.0, mass=1.0)
method = DOPRI5(rtol=1e-10, atol=1e-12)

solution = solve_pendulum(
    theta0=1.0,
    omega0=0.0,
    end_time=3 * harmonic_period(parameters),
    parameters=parameters,
    method=method,
)

angle_at_one_second = solution(1.0)[0]
period = pendulum_period(1.0, parameters=parameters, method=method).period
```

## Tests

Run all tests (technically a single one):

```bash
uv run pytest
```

## Generate the results and plots

```bash
uv run python -m scripts.hw3_mpendulum
```

## Compile the report

```bash
typst compile hw3.typ
```
