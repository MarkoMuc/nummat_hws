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

= Numerical experiments

The experiments use $ g=9.81 "m/s"^2, quad l=1 "m", quad m=1 "kg"$
and initial angular velocity $omega_0=0$.
The harmonic period associated with these parameters is $ T_0=2 pi sqrt(l/g) approx 2.00607 "s". $

First we will compare the mathematical pendulum is approximated with DOPRI5.
But the harmonic approximation replaces $sin(theta)$ by $theta$:
$ theta_h'' + g/l theta_h=0. $

For $omega_0=0$, the harmonic solution is
#footnote[#link("https://en.wikipedia.org/wiki/Pendulum_(mechanics)#Small-angle_approximation")]
$ theta_h(t)=theta_0 cos(sqrt(g/l)t).$
Plus the previously mentioned harmonic period.

The mathematical pendulum is approximated for three initial angles $theta_0=0.1$, $1.0$, and $2.5$.
Every solution is computed over three harmonic periods.

DOPRI5 chooses its time steps adaptively, so different solutions do not generally contain values at the same times.
For plotting, we evaluate them on one common uniform grid.
Intermediate values are obtained with cubic Hermite interpolation, as in `Vaje 16`
#footnote[#link("https://gitlab.com/nummat/nummat-knjiga/-/blob/master/Vaja16/src/Vaja16.jl?ref_type=heads")].
The comparison is plotted on @img1.

#figure(
  image("img/hw3_pendulum_comparison.svg", width: 88%),
  caption: [
    Comparison of the nonlinear mathematical pendulum approximated with DOPRI5 and
    the harmonic approximation for three initial angles.
  ],
)<img1>

For the smallest angle $theta_0=0.1$, the two curves are almost identical.
This is expected because $sin(theta) approx theta$ near zero.

At $theta_0=1.0$, the difference in the periods becomes visible and gradually produces a shift.

At $theta_0=2.5$, the harmonic approximation is no longer accurate, as it oscillates considerably faster than the nonlinear pendulum.
For $0<theta<pi$, we have $sin(theta)<theta$,
so the restoring acceleration in the nonlinear model is weaker than in the harmonic model.
The nonlinear pendulum consequently needs more time to complete one oscillation.

The second experiment is the comparison between the period of the nonlinear pendulum and the energy of the system.
Note that the graph is not constructed from all time steps of one trajectory.
Instead, we repeat the complete simulation for many different initial angles.
Each initial angle produces one energy and one measured period, thus one point $(E,T)$ in the graph.

The mechanical energy is
#footnote[#link("https://en.wikipedia.org/wiki/Pendulum_(mechanics)") and #link("https://en.wikipedia.org/wiki/Spherical_pendulum")]
$ E(t)=1/2 m l^2 omega(t)^2 + m g l(1-cos(theta(t))). $

Because every pendulum is released from rest, $omega_0=0$, its initial energy is simplifies to $ E_0=m g l(1-cos(theta_0)). $

The model contains no damping, so this total energy should remain constant during one simulation.
Consequently, we calculate one energy value for each initial condition rather than treating the DOPRI5 time points as different energy levels.
The graph uses the dimensionless quantity $E/(m g l)$.

*Determining the period*

For every initial angle, DOPRI5 first produces a numerical trajectory $(theta(t),omega(t))$.
We determine its period by detecting zero crossings of $theta(t)$.
Starting from a positive turning point with $omega_0=0$, the pendulum first crosses $theta=0$ while moving downwards.
After one complete oscillation it crosses $theta=0$ downwards again.
The difference between the times of these two crossings is therefore the period: $ T=t_2-t_1. $

The zero search is performed in two short stages.
First, consecutive stored DOPRI5 states are inspected until the angle changes sign, i.e., we bracket a zero
The bracket is then narrowed with bisection.
At each trial time the state is reconstructed from the left endpoint with a DOPRI5 step.

The procedure is repeated for initial angles between $0.05$ and $0.95 pi$.
For every angle we store $ (E_0,T). $

The plotted graph can be seen on @img2.

#figure(
  image("img/hw3_period_energy.svg", width: 82%),
  caption: [
    Period of the mathematical pendulum as a function of normalized initial energy.
    The dashed line is the constant harmonic period $T_0$.
  ],
)<img2>

At low energy, the mathematical pendulum is close to the harmonic model, so its period is close to $T_0$.
As the initial energy and thus amplitude increase, the nonlinear period increases.
The harmonic line stays constant because the harmonic period is independent of amplitude and energy.
Near $E/(m g l)=2$, the initial angle approaches the unstable upright position $theta_0=pi$.
The pendulum then spends increasingly more time near the turning point, which explains the steep rise of the curve.
