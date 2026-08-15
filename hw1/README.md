# Homework 1: Conjugate Gradient method for Sparse matrices

**Author:** Marko Zupančič Muc

Implement the conjugate gradient method for sparse matrices.
The sparse matrix is represented with a special Data type `RedkaMatrika`.
Conjugate gradient method can then use the sparse matrix type to solve the linear system.

At the end we also use conjugate gradient method to embed a ledder graph same as in lab `Vaje 06`.

## Implementation

- `src/sparse_matrix.py` contains the type implementation of the sparse matrix `RedkaMatrika`.
- `src/conjugate_gradient.py` contains the implementation of the conjugate gradient method.
- `src/graph_embedding.py` contains the implementation of the ladder graph generating function and embedding logic.
  - Most of the logic is from Chapter 6 of the `nummat_jl.pdf` book.

## Tests

Run all tests:
```bash
uv run pytest
```

To only run a specific test:

```bash
uv run pytest tests/test_conjugate_gradient.py
```

## Generate the plots

```bash
uv run python -m scripts.demo
```

## Compile the report

```bash
typst compile hw1.typ
```