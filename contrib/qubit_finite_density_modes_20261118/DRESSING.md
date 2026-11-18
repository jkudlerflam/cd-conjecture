# Raw density, dressed charge and the virtual-pair background

## Setting and order of the assertion

Use `n=(1-Z)/2`, `b=(X+iY)/2`, `N=sum n`,
`D_r=sum_j n_j n_(j+r)` and `T_r=sum_j b_(j+r)^dagger b_j`.
The orientation of T does not affect the real T2 combination below. Work in
the periodic bulk, with P the no-adjacent-particle subspace and Q=1-P.
The large term is `V D_1`, the pair perturbation is
`Delta(A+A^dagger)`, `A=sum_j b_j^dagger b_(j+1)^dagger`.
At the prescribed point Delta=1, V=40, epsilon=Delta/V=1/40.

The actual low-energy hopping and bulk chemical potential are O(epsilon^2 V).
Let `R=Q D_1^(-1) Q`, and choose the first canonical off-block generator
`S=epsilon (R A P - P A^dagger R)+...`, so `e^S H e^(-S)` has no first-order
P--Q pair matrix elements. Higher diagonal-block rotations change operator
representatives; this convention fixes the coefficient stated here.

**[proved, order by order]** Since A adds two particles and cannot annihilate
a pair from P,

\[
P e^S N_{\rm raw}e^{-S}P
=N+2\epsilon^2 B+O(\epsilon^3),\qquad
B=P A^\dagger D_1^{-2} A P.                         \tag{1}
\]

The displayed remainder denotes formal perturbative order in this scaling.
It is **not** an L-independent norm bound at epsilon=1/40. A global spectral
gap isolating the entire extensive P band has not been assumed or proved.
There may be higher-order charge-changing pieces in the effective observable.

## Exact local coefficient

For a no-adjacent configuration, creating a pair on the empty bond (j,j+1)
costs `V(1+n_(j-1)+n_(j+2))`. Thus the squared denominators in B are
1, 1/4 and 1/9, retaining *both* outer environments.

The diagonal contribution is the sum of

\[
(1-n_j)(1-n_{j+1})
\left[1-\frac34(n_{j-1}+n_{j+2})
       +\frac{11}{18}n_{j-1}n_{j+2}\right].         \tag{2}
\]

Using the hard constraint and summing translates gives
`L - (7/2)N + (3/2)D_2 + (11/18)D_3`.
Two different creation/annihilation bonds can share an intermediate state
only by moving an existing particle by two sites. This state has two
adjacent occupied bonds, giving the off-diagonal coefficient 1/4. Hence

\[
B=P\left[L-\frac72N+\frac32D_2+\frac{11}{18}D_3
             +\frac14(T_2+T_2^\dagger)\right]P.   \tag{3}
\]

Multiplication by 36 makes every entry an integer. `dressing.py` independently
constructs `36 A^dagger D_1^(-2) A` in the literal spin Hilbert space and the
local right side, then checks entrywise equality at L=8,9,10,12. This finite
check is corroboration; (2) and the shared-intermediate-state argument prove
the bulk local identity, without extrapolating those ring sizes.

Positivity follows from the first expression in (1). Also `D_r<=N` and
`T_2+T_2^dagger<=2N`, the latter from the sum of
`(b_j-b_(j+2))^dagger(b_j-b_(j+2))>=0`. Consequently

\[
0\le B\le L-\frac89N\le L.                       \tag{4}
\]

This is an operator bound **on the leading coefficient**. At density about
0.215 the leading excess raw density is between zero and 0.001011; the
omitted orders are not bounded by that interval.

## A background in number fluctuations

In a definite-N leading effective state, the virtual states all have N+2.
The same commutator expansion, or their probabilities, gives

\[
\langle N_{\rm raw}\rangle-N=2\epsilon^2\langle B\rangle+O(\epsilon^3),
\quad
\operatorname{Var}N_{\rm raw}=4\epsilon^2\langle B\rangle+O(\epsilon^3).
                                                               \tag{5}
\]

This contribution is short range in the effective description. It produces
an analytic q=0 background in the raw equal-time density structure factor;
the nonanalytic |q| coefficient of a matched smooth density is a different
quantity. There is no microscopic U(1) Ward identity that makes the raw
spin count an exactly normalized Gaussian charge. In particular a fit of
`S_raw(q)/|q|` through the origin can be biased even when epsilon is small.
Equation (5) does not determine a static compressibility: virtual spectral
weight is additionally divided by its excitation energy in a response.

## What can be evaluated from the saved medium

The C193 p153c records contain four central density-correlation anchors and
the short-distance complex one-body correlator. Inserting them in the bulk
representative (3) gives the following **plug-in estimates**, at L=256,
chi=768:

| Parity | Mean of four leading excess-density representatives | Anchor range |
|---|---:|---:|
| even | 0.00054724 | 0.00032392 to 0.00074045 |
| odd | 0.00055778 | 0.00033931 to 0.00074676 |

These numbers include the D2, D3 and real hopping terms; the vacuum
denominator alone would give 0.00125. They are useful scale estimates, not
measured dressed densities. The identity is a *summed bulk* identity;
different local representatives differ by telescoping terms. Four strongly
oscillating central anchors are not a full bulk average. Moreover these are
raw physical, not separately SW-matched, correlators, at one chi. There is
no controlled higher-order or thermodynamic error bar, and no split into
rho_A and rho_B follows from the total correction.

The open-chain endpoint convention matters. In the bulk,
`mu=4Jz+2hz=0.0482175`; on each end it is `mu-V/2` because one ZZ neighbour
is missing. Thus the end potential differs by 20 in physical energy units.
Equations (2)--(4) apply in the bulk; they must not be pasted onto the end
bonds without their different denominators.
