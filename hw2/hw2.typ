#set document(
  title: "Barycentric interpolation and Gaussian quadrature",
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
  #text(17pt, weight: "bold")[Barycentric interpolation and Gaussian quadrature]
  #v(0.5em)
  #text(11pt)[Homework 2]
  #v(0.8em)
  Marko Zupančič Muc
]

#v(1em)

= Interpolation with the barycentric formula

Polynomial interpolation constructs a polynomial $p(x_i)$ that matches a given function at selected nodes.
The interpolating polynomial $p_n$ is the unique polynomial of degree at most $n$ that satisfies
$ p_n (x_k) = y_k, $ for nodes $x_0, ..., x_k$ and sampled values $y_k = f(x_k)$.
Interpolation in the standard basis uses the Vandermonde matrix.
However, this can be computationally expensive and numerically unstable, especially for high-degree polynomials or poorly chosen nodes.
Interpolation costs approximately $2/3n^3 O(n^2)$ operations.

To solve these issues, we move from the standard basis to the more stable Lagrange basis.
The Lagrange form describes the same polynomial directly in terms of the sampled valyes.
For each node $x_k$, we define a Lagrange basis polynomial as
$ l_k (x) = product^n_(j=0, j!=k) (x - x_j)/(x_k - x_j). $

It satisfies:
$ l_k (x_j) = cases( 1 "if" j = k ",", 0 "else.") $

The interpolating polynomial for $j = k$ is then:
$ p_n (x) = sum^n_(k=0)y_k l_k (x). $

For each node $x_i$, all terms except the $i$-th vanish: $ p_n (x_i) = sum^n_(k=0)y_k l_k (x_i) = y_i. $

Polynomial interpolation constructs a polynomial that matches a given function at selected nodes.
The Lagrange form provides an explicit expression for this polynomial, but evaluating it directly repeats many calculations.
The barycentric form reorganizes the same polynomial so that reusable weights can be computed in advance,
making repeated evaluation more efficient.

The cost of evaluating one $l_k (x)$ consists of multiplying approximately $n$ factors.
There are $n+1$ basis polynomials, so directly evaluating the Lagrange form costs $O(n^2)$ operations per argument.

The barycentric form avoids explicitly calculating every Lagrange basis polynomial $l_k (x)$.
We use the following formula for barycentric Lagrange interpolation:
$
  l(x) = cases(
      sum((f(x_j) lambda_j)/x-x_j)/sum(lambda_j/(x - x_j)) "if" x!=x_j ",",
      f(x_j) "else."
    ).
$

== Chebyshev nodes

Polynomial interpolation using Chebyshev nodes minimizes the effect of Runge's phenomenon.
Runge's phenomenon refers to oscillations near the edges of an interval that can occur when high-degree polynomial interpolation
uses a set of equidistant interpolation points.

In our case, we interpolate using Chebyshev nodes generated on the interval $[-1, 1]$ with the formula
$ x_k = cos((2k - 1)/(2n) pi), quad k = 0,1,...,n-1. $

We also need the Chebyshev weights $lambda_j$ for barycentric interpolation. The weights are given by the following formula:
$ 
  lambda_k = (-1)^k cases(
    1 "if" 0 < 1 < n ",",
    1/2 "if" i = 0 ",",
    n "else."
  )
$

We map the interval from $t in [-1,1]$ to $x in [a,b]$ using an affine transformation
#footnote[#link("https://en.wikipedia.org/wiki/Chebyshev_nodes#Definition")]:
$ x = (a + b)/2 + (b-a)/2 t. $

*NOTE*: There seems to have been a mismatch between the formulas for the Chebyshev weights and nodes.
The provided node formula uses Chebyshev points of the first kind,
whereas the provided weight formula uses Chebyshev weights of the second kind.
We therefore use formulas of the second kind for both:

$ 
  x_k = cos(k/n pi), quad k = 0,1,...,n.\

  lambda_k = (-1)^k cases(
    1/2 "if" k = 0 "," n",",
    1 "else."
  )
$

== Interpolation of the example functions

The goal was to determine a polynomial degree for which the maximum interpolation error is below $10^(-6)$.
For the first two functions, we tested consecutive degrees and found the first one that met the target.
The third function required a much larger degree,
so we estimated it from the observed convergence rate and verified the result at $n=1.2 dot 10^6$.
The results are shown in @tabela.

#figure(
  table(
    columns: (2.1fr, 1.2fr, 0.7fr, 1.4fr, 1.4fr),
    align: (left, center, center, right, right),
    table.header([$f(x)$], [interval], [$n$], [Error at $n-1$], [Error at $n$]),
    [$e^(-x^2)$], [$[-1,1]$], [10], [$2.0811 dot 10^(-5)$], [$8.2919 dot 10^(-7)$],
    [$sin(x)/x$], [$[0,10]$], [13], [$4.4543 dot 10^(-6)$], [$1.3770 dot 10^(-7)$],
    [$abs(x^2-2x)$], [$[1,3]$], [1200000], [$1.1938 dot 10^(-6)$], [$9.9487 dot 10^(-7)$],
  ),
  caption: [Degrees used to reach an error below $10^(-6)$. For the third function, the earlier error was measured at $n=10^6$.],
)<tabela>

@abs_err shows the absolute interpolation error at the degrees listed in @tabela.

#pagebreak()

#figure(
  image("./img/hw2_interpolation_errors.svg", width: 92%),
  caption: [The curve shows the absolute interpolation error, and the dotted line marks the target of $10^(-6)$.],
)<abs_err>

The functions and their interpolants are compared in @graphs.

#figure(
  image("./img/hw2_interpolants.svg", width: 92%),
  caption: [The blue curves show the original functions, while the red points are values of the corresponding interpolants.],
)<graphs>

For the first two smooth functions,
the error decreases quickly as the polynomial degree increases because both functions can be approximated well by polynomials.
However, the function $g(x)=|x^2-2x|$ has a corner at $x=2$, where its first derivative is discontinuous:
$ g'_r(2^-)= -2 quad g'_r(2^+)= 2. $
I.e., the slope jumps from $-2$ to $2$, which produces a corner.
This matters because polynomials are smooth and have a continuous derivative everywhere.
To reproduce the corner, the polynomial starts bending rapidly near $x=2$, resulting in a large increase in the required degree.

In the experiments, the maximum error here decreases approximately as
$ E_n approx C/n. $

The experiments give $C approx 1.194$#footnote[This can be seen in the convergence results].
To satisfy the required error $E_n < 10^(-6)$, we therefore need approximately
$ n > C / 10^(-6) approx 1.194 dot 10^6 approx 1.2 "million". $

This explains the large degree required for the third function.
In @conver, the measured errors approximately follow a line proportional to $1/n$, confirming the much slower convergence.

#figure(
  image("./img/hw2_convergence.svg", width: 82%),
  caption: [Convergence of the error as the polynomial degree increases.],
)<conver>
