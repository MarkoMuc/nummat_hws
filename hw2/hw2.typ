#import "@preview/ctheorems:1.1.3": *
#show: thmrules

#let definition = thmbox("definition", "Definition", inset: (x: 1.2em, top: 1em))

#let proof = thmproof("proof", "Proof")

#set math.equation(numbering: "(1)")


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

#outline()
#line(length: 100%)

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

The experiments give $C approx 1.194$#footnote[This can be seen in the convergence results found in `./results/hw2_part1.json`.].
To satisfy the required error $E_n < 10^(-6)$, we therefore need approximately
$ n > C / 10^(-6) approx 1.194 dot 10^6 approx 1.2 "million". $

This explains the large degree required for the third function.
In @conver, the measured errors approximately follow a line proportional to $1/n$, confirming the much slower convergence.

#figure(
  image("./img/hw2_convergence.svg", width: 82%),
  caption: [Convergence of the error as the polynomial degree increases.],
)<conver>

#pagebreak()

= Gauss–Legendre quadrature

An integral can be approximate as a weighted sum
$ integral_a^b f(x) d x approx sum^r_(k=1) w_k f(x_k), $ where $x_k$ are knots and $w_k$ are the weights.

Unlike trapez or samson rule, Gauss positions the knots and can approximate higher order functions (why?).

Deriving the Gauss-Legendre quadrature rule on two points:
$ integral^1_(-1) f(x) d x = A f(x_1) + B f(x_2) + R_f, $
including the derivation of the error $R_f$ formula.

We take polynomials, and pick $w_i$ and $x_i$ such that the formula is accurate for polynomials of higher order.

== Deriving the two point formula on $[-1,1]$

#definition[
  The two point formula for Gauss-Legendre quadrature on $[-1,1]$ is:
  $ integral_(-1)^1 f(x) d x = f(- 1/sqrt(3)) + f(1/sqrt(3)) + R_f. $

  Where the error formula is $ R_f = 1/135 f^(4)(xi), quad xi in (-1 ,1). $
]<deg1>

#proof([of @deg1])[
  Start off with the 2-point rule:
  $ integral_(-1)^1 f(x) d x approx w_1f(x_1) + w_2f(x_2). $

  There are four unknowns: $w_1, w_2, x_1, x_2.$
  We will insert successive polynomials to find the solution.

  For degree $0$ we have $f(x) = 1$:
  $ integral_(-1)^1 1 dif x = 2 = w_1 + w_2. $

  Hence $w_1 + w_2 = 2.$

  For degree $1$ we have $f(x) = x$:
  $ integral_(-1)^1 x dif x = 0 = w_1 x_1 + w_2 x_2. $

  Hence $w_1 x_1 + w_2 x_2 = 0.$

  For degree $2$ we have $f(x) = x^2$:
  $ integral_(-1)^1 x^2 dif x = 2/3 = w_1 x_1^2 + w_2 x_2^2. $

  Hence $ w_1 x_1^2 + w_2 x_2^2 = 2/3. $

  For degree $3$ we have $f(x) = x^3$:
  $ integral_(-1)^1 x^3 dif x = 0 = w_1 x_1^3 + w_2 x_2^3. $

  Hence $w_1 x_1^3 + w_2 x_2^3 = 0.$

  The full system of equations is:
  $ 
    w_1 + w_2 &= 2,\
    w_1 x_1 + w_2 x_2 &= 0,\
    w_1 x_1^2 + w_2 x_2^2 &= 2/3,\
    w_1 x_1^3 + w_2 x_2^3 &= 0.
  $

  We can now solve the system.
  Using the fact that there are two points on the interval $[-1,1]$ and our aim is to make it symmetric,
  we deduce that $x_2 = -x_1$ and thus also $w_1 = w_2$.

  Now use the $x^2$ equation:
  $ x_1^2 + x_2^2 = 2/3. $

  Insert $x_2 = -x_1$,
  $ 2 x_1^2 = 2/3\
  x_1^2 = 1/3.
  $

  Thus we can choose
  $ x_1 = -1/sqrt(3), quad x_2 = 1/sqrt(3) . $

  The complete 2-point formula is:
  $ integral_(-1)^1 f(x) dif x approx f(-1/sqrt(3)) + f(1/sqrt(3)) + R_f . $

  The error is calculated as $ R_f = integral_(-1)^1 f(x) dif x - f(-1/sqrt(3)) + f(1/sqrt(3)). $

  We already know that Gauss-Legendre is exact for every polynomial of degree $<= 3$:
  Therefore the first derivative that can contribute to the error is the fourth derivate.

  The first polynomial where the error can appear is in $f(x) = x^4$.
  We can once again insert it as before, yielding:
  $ integral_(-1)^1 x^4 dif x = 2/5 != 2/9 = (-1/sqrt(3))^4 + (1/sqrt(3))^4, $
  unlike for degrees 0 to 3, they are not equal.

  Their difference is the error:
  $ R_(x^4) = 2/5 - 2/9 = 8/45. $

  This is only the error for degree 4, but we need the general error formula for arbitrary $f$.
  The general error theorem has the form $ R_f = C f^((4)) (xi), quad "for some" xi in (-1, 1). $

  We can get $C$ by using $f(x) = x^4$, we already know $R_(x^4) = 8/45$ and $f^((4))(x) = 24.$
  Therefore
  $ 8/45 = C dot 24\
    C = 1/135.
  $

  Finally we get the general formula $ R_f = 1/135 f^((4)) (xi). $
]

== Deriving the two point formula on $[a,b]$

We can go directly from $[-1,1]$ to an arbitrary interval $[a,b]$ using the linear substitution

$ x(t) = (b-a)/2 t + (a+b)/2, quad t(x) = 2 (x-a)/(b-a). $

Integral on $[a,b]$ can be expressed with the integral on $[-1,1]$ as:
$ integral_a^b f(x) dif x = integral_(-1)^1 f(x(t)) x'(t) dif t = 1/2 (b-a) integral_(-1)^1 f(x(t)) dif t, $
where $x'(t) = 1/2 (b-a)$.

As such the Gauss-Legendre formula on $[a,b]$ is:
$ integral_a^b f(x) dif x = 1/2 (b-a) integral_(-1)^1 f(x(t)) dif t approx (b-a)/2 (f(x(-1/sqrt(3))) + f(x(1/sqrt(3)))). $

Applying the 2-point rule on $[-1,1]$, where $t_1 = -1 1/sqrt(3), quad t_2= 1/sqrt(3)$:
$ integral_a^b f(x) dif x  approx (b-a)/2 (f((a+b)/2 - (b-a)/(2 sqrt(3))) + f((a+b)/2 + (b-a)/(2 sqrt(3)))). $

Again we determine the error formula $R_f$ by testing degree $4$ polynomial.
We already have the rule
$ (b-a)/2 (f(m - (b-a)/(2 sqrt(3))) + f(m + (b-a)/(2 sqrt(3)))), $
where $m = (a+b)/2$.

In this case we choose the polynomial $f(x) = (x-m)^4$, because the interval is symmetric around $m$.

Now let $ d = (b-a)/2, quad  a=m -d, quad b = m + d. $

Then the interval is effectively $m-d$ to $m+d$.
We introduce the substitution $y = x - m,$ then the integration limits are $a -> y = (m - d) - m = -d$ and $b -> y = (m+d) -m = d$.

The exact integral is $ integral_a^b (x-m)^4 dif x = integral_(-d)^d y^4 dif y = [y^5/5]^d_(-d) = (2 d^5)/5. $

Now apply the Gaussian rule.
The two nodes are $m plus.minus d/sqrt(3)$, so
$ d ( (-d/sqrt(3))^4 + (d/sqrt(3))^4 ) = d (d^4/9 + d^4/9) = (2 d^5)/9. $

Therefore the error for this fourth-degree polynomial is $ R_f = (2 d^5)/5 - (2 d^5)/9 = (8 d^5)/45. $

Now using the general error form $R_f = C f^((4)) (xi)$ and the fourth derivative $f^((4)) (x) = 24$ we get:
$ (8 d^5)/45 = C dot 24,\
  C = d^5/135 = 1/135 dot ((b-a)/2)^5 = (b-a)^5/4320.
$

Thus $ R_f = (b-a)^5/4320 f^(4)(xi), quad xi in (a,b). $

== Deriving the two point formula on $[0,h]$

This is just a special case of $[a,b]$ where $ a =0, quad b = h $.

From $ integral_a^b f(x) dif x  approx (b-a)/2 (f((a+b)/2 - (b-a)/(2 sqrt(3))) + f((a+b)/2 + (b-a)/(2 sqrt(3)))), $ set $a=0, b=h$.

Then $ (a+b)/2 = h/2, quad (b-a)/2 = h/2. $

Therefore $ integral_a^b f(x) dif x approx h/2 (f(h/2 - h/(2 sqrt(3))) + f(h/2 + h/(2 sqrt(3)))). $
Which is equivalent to $ integral_a^b f(x) dif x approx h/2 (f(h/2 (1 - 1/sqrt(3))) + f(h/2 (1 + h/sqrt(3)))). $

The error follows the same way:
$ R_f = (b-a)^5/4320 f^((4)) (xi), $ using $a=0, b=h$: $ R_f = h^5/4320 f^((4)) (xi), quad xi in (0,h). $

#pagebreak()
== Composite two point formula

Here we will double the number of steps.

Divide $[a,b]$ into $n$ equal subintervals and let
$ h = (b-a)/n, quad x_i = a + i h. $

On each subinterval $[x_i, x_(i+1)]$, apply the two-point Gauss–Legendre rule.
Since the midpoint of the subinterval is $(x_i + x_(i+1))/2 = x_i + h/2$, we get
$ integral_a^b f(x) dif x approx Q_n = h/2 sum_(i=0)^(n-1) ( f(x_i + h/2 - h/(2 sqrt(3))) + f(x_i + h/2 + h/(2 sqrt(3))) ). $
So we split the entire integral into all $n$ subintervals, and used Gauss approximation on each subinterval, summing the results.

Hence $ integral_a^b f(x) dif x =  Q_n + R_n. $
For one subinterval, the error is $ R_i = h^5/4320 f^(4)(xi_i), $ for some $xi_i in (x_i, x_(i+1))$.

Adding the errors from all subintervals gives
$ R_n = h^5/4320 sum_(i=0)^(n-1) f^(4)(xi_i). $

If $f^(4)$ is continuous, then for some $xi in (a,b)$ we can write $sum_(i=0)^(n-1) f^(4)(xi_i) = n f^(4)(xi).$

Plugging this in into the formula for all subintervals:
$ R_n = n h^5/4320 f^(4)(xi). $

We can further simplify with $n h = b-a$ and the definition of $h$, we can rewrite it as
$ R_n = (b-a)^5/(4320 n^4) f^(4)(xi). $

Therefore the composite two-point Gauss–Legendre rule is fourth order:
$ R_n = O(h^4) = O(n^(-4)). $

=== Error estimation by doubling the number of subintervals

Since the composite rule is fourth order, its error behaves approximately as
$ R_n approx C h^4. $

If we double the number of subintervals from $n$ to $2n$, then the step size is halved $ h -> h/2$.

Therefore $ R_(2n) approx C (d/2)^4 = 1/16 C d^4 approx 1/16 R_n. $

Let $I$ denote the exact value of the integral.
Since $Q_n = I - R_n$ and $Q_(2n) = I - R_(2n),$ we can subtract the formulas:
$ Q_(2n) - Q_n = (I - R_(2n)) - (I - R_n) = R_n - R_(2n). $

Using $R_n approx 16 R_(2n)$, we obtain $ Q_(2n) - Q_n approx 16 R_(2n) - R_(2n) = 15 R_(2n). $

Therefore $ R_(2n) approx (Q_(2n) - Q_n)/15. $

Wrap absolute value, yielding $abs(R_(2n)) approx abs(Q_(2n) - Q_n)/15.$
Since $R_(2n) = I - Q_(2n)$, the error of the finer approximation is estimated by

$ abs(I - Q_(2n)) approx abs(Q_(2n) - Q_n)/15. $

== Example integral


