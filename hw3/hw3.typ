#import "@preview/ctheorems:1.1.3": *
#show: thmrules

#let definition = thmbox("definition", "Definition", inset: (x: 1.2em, top: 1em))

#let proof = thmproof("proof", "Proof")

#set math.equation(numbering: "(1)")


#set document(
  title: "Mathematical pendulum",
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
  #text(17pt, weight: "bold")[Mathematical pendulum]
  #v(0.5em)
  #text(11pt)[Homework 3]
  #v(0.8em)
  Marko Zupančič Muc
]

#v(1em)

#outline()
#line(length: 100%)

= Mathematical model of the pendulum

The angle movement of the mathematical pendulum is described with the equation:
$ theta''(t) + g/l sin(theta(t)) = 0, quad theta(0)=theta_0, quad theta'(0)=omega_0, $
where $l$ is the string length and $g$ is the gravitational acceleration.

To evaluate the function, we will use the DOPRI5 method, which is a version of the Rugga-Kutta methods.
But Runge-Kutta methods work on first order differential equations, while our equation is of the second order.

== Converting to first-order system

Let $ omega(t) = theta'(t) $ be the rate at which the angle changes (angular velocity).

Then $ w'(t) = dot.double(theta) = - g/l sin(theta) $

The result is the following first-order system 
$
  theta'(t) = omega(t),\
  omega'(t) = - g/l sin(omega(t)).
$

Combine the two variables $theta, omega$ into one vector $ u(t) = vec(theta(t), omega(t)). $

Runga-Kutta methods are designed for equations that look like $ u(t) = f(t, bold(u)). $
This system can also be represented in the following vector form
$ bold(u)'(t) = bold(f)(t, bold(u), bold(p)) = vec(omega, - q/l sin(theta)), quad bold(u)(0) = vec(theta_0, omega_0). $
where the parameter parameters combined into $bold(p) = (g,l)$.

== DOPRI5

Runge-Kutta methods, unlike e.g. Euler method, compute several intermediate slopes.
R-K methods, when approximation between points $x_n -> x_(n+1)$ takes into account multiple derivatives between the two points,
and creating a weighted sum (average): $ beta_1 f(x_n, y_n) + beta_2 f(x_n + c_2 h, y_n + d_2 k) + beta_3 f(x_n + c_3 h,y_n + d_3 k) + ... $

In general a single slope is calculated as
$ k_i = f( t_n + c_i h, bold(y)_n + h sum_(j < i) a_(i j) k_j), $
then the slopes are combined into a single weighted sum
$ bold(y)_(n+1) = bold(y)_n + h sum_i b_i k_i, $
where $b_i$ are the weights.

The parameters $c_i$ and $a_i$ are often stored together in a Butcher tableau.

DOPRI5 is an embedded seven stage Runge-Kutta method that calculates two approximations at the same time.
One approximation is of order 5 and the other is of order 4.
The difference between the approximation gives us an estimate of the numerical error:
$ e = y^(5)_(n+1) - y^(4)_(n+1). $

If the difference is too big, the step was too large and be repeated using a smaller $h$.
Thus DOPRI5 can adapt its time step.

The Butcher tableau for the DOPRI5 method can be found here #link("https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method#Butcher_tableau").

== Evaluating and comparing the mathematical and harmonic pendulum

The harmonic pendulum formula can be directly sovled using the following stuff:


To track the period and to interpolate for plotting, we will use the Hermitic interpolation and zero finding as in Vaje 16.


