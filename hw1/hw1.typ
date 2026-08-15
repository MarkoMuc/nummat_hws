#set document(
  title: "Conjugate gradient method for sparse matrices",
  author: "Marko Zupančič Muc",
)
#set page(paper: "a4", margin: (x: 2.4cm, y: 1.8cm), numbering: "1")
#set text(lang: "eng", size: 10.5pt)
#set par(justify: true, leading: 0.7em)
#set heading(numbering: "1.")
#show raw.where(block: true): set block(
  fill: luma(245),
  inset: 8pt,
  radius: 4pt,
  width: 100%,
)
#show link: set text(fill: rgb("2457a6"))

#align(center)[
  #text(17pt, weight: "bold")[Conjugate gradient method for sparse matrices]
  #v(0.5em)
  #text(11pt)[Homework 1]
  #v(0.8em)
  Marko Zupančič Muc
]

#v(1em)

= Conjugate gradient method

The conjugate gradient method is an iterative method for solving certain systems of linear equations.
Unlike Gaussian elimination, which eliminates variables, it improves an initial guess by moving along carefully chosen directions.
This produces a sequence of approximations $ x_0, x_1, x_2, ... $ that approaches the exact solution $x^*$ of $ A x = b. $

The method is designed for matrices $A$ that are symmetric, $A^T = A$, and positive definite, meaning that $x^T A x > 0$ for every nonzero $x$.

We start off at an initial approximation $x_0$, iteratively improving.
Thus we need some sort of metric to check when we are close enough to $x^*$.

To measure the quality of an approximation $x_k$, we use its residual
$ r_k = b - A x_k. $
A small residual means that $x_k$ nearly satisfies the original system.

Another way to view the method is through the quadratic function
$ f(x) = 1/2 x^T A x - b^T x. $
Solving $A x = b$ is equivalent to minimizing this function.

Because $A$ is symmetric, $gradient f(x) = A x - b$.
Therefore, $ gradient f(x) = 0 <=> A x = b. $

Since $gradient f(x_k) = A x_k - b$, we have $ r_k = - gradient f(x_k). $
The residual therefore points in the direction of the steepest decrease of $f$.

Two directions $p_i$ and $p_j$ are conjugate with respect to $A$ if $ p_i^T A p_j = 0. $
A set of $n$ nonzero, mutually conjugate directions forms a basis of $RR^n$.
This property allows the method to improve the approximation in one direction without undoing the progress made in the previous directions.

== The iterative method

We begin with an initial guess $x_0$ and compute the first residual $ r_0 = b - A x_0. $
The first search direction is $ p_0 = r_0. $

#pagebreak()
The approximations are updated according to the recursive expression
$ x_(k+1) = x_k + alpha_k p_k. $

We choose $alpha_k$ to minimize $f$ along the line through $x_k$ in the direction $p_k$:
$ alpha_k = (r_k^T r_k)/(p_k^T A p_k). $
This step is called a line search.

We then calculate the next residual:
$ r_(k+1) = r_k - alpha_k A p_k. $

Instead of using this residual directly as the next direction, we combine it with the previous direction:
$ p_(k+1) = r_(k+1) + beta_k p_k, $
where $beta_k$ is chosen so that $p_(k+1)^T A p_k = 0.$
This gives $ beta_k = (r_(k+1)^T r_(k+1))/(r_k^T r_k). $

We repeat these steps until the residual is sufficiently small.

== Sparse matrix

A sparse matrix contains mostly zero elements.
By storing only its nonzero elements, we can reduce both memory use and computation time of some operations.
For each row $i$, our representation stores a list of values $V[i]$ and the corresponding column indices $I[i]$, so that 
$ V[i][j] = a_(i,I[i][j]). $
If $A$ has $m$ nonzero elements, this representation requires $O(m)$ space as opposed to $O(n^2)$.

This representation also makes the matrix vector product $A p_k$, the main operation in each conjugate gradient iteration, more efficient.
Because only the nonzero elements take part in the calculation, a product can be computed in $O(m)$ operations:
$ (A x)[i] = sum_j V[i,j] x_(I[i,j]). $

For a dense matrix, the same operation would require $O(n^2)$ operations.

The product $A p_k$ needs to be computed only once per iteration.
The remaining work consists of dot products and vector updates, each requiring $O(n)$ operations.

The dominant step in an iteration is therefore the product $A p_k$, which requires $O(m)$ operations.
Since a typical sparse matrix has at least as many nonzero elements as rows, the overall cost of one iteration is $O(m)$.

= Graph embedding

The last part of the task is to use the new sparse matrix representation and
the conjugate gradient method to the graph embedding problem from Chapter 6 of the lab book.
The graph related functions were adapt from the lab work and now use the new matrix representation and solver where appropriate.

The embedding uses the physical model described in Chapter 6.
Neighboring vertices are connected by springs, and we seek a position in which the forces on every free vertex are balanced:
$ sum_(i in N(j)) F_(i j) = 0. $
Some vertices must be fixed at known coordinates to make the problem well-defined.
The equilibrium conditions for the remaining vertices form a linear system $A x = b$, where $x$ contains their unknown coordinates.

The matrix $A$ contains the degrees of the free vertices on its diagonal.
Its off-diagonal elements describe which free vertices are neighbors.
After the fixed vertices are taken into account, the matrix is symmetric, positive definite, and usually sparse,
so it is well suited to the conjugate gradient method.

We apply the method to a circular ladder graph with 16 vertices. The result is shown in @figure1.

#figure(
  image("./img/circular_ladder_embedding.svg"),
  caption: [
    The left image shows the initial assignment. The outer ring is fixed, while all vertices in the inner ring initially lie at $(0,0)$
    and are considered free.
    The right image shows the graph after we've embedded it in the plane with the conjugate gradient method.
  ],
)<figure1>

To verify the result, we reconstruct the system matrix $A$ and calculate the residual $r = A x - b$.
If the embedding satisfies the equilibrium conditions, the residual norm should be close to zero.

The computed values confirm this:

```text
Vertices: 16
Fixed vertices: 8
System shape: (8, 8)
Nonzero elements: 24
Iterations: [1, 1]
Residual norms: [4.991976076032729e-16, 8.97892199319327e-16]
```
