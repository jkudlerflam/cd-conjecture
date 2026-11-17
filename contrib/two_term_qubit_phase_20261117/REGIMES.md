# What weak and strong coupling do establish

## Flat point and the first generated velocity

**[Proved at g=0; derived first-order perturbation for g!=0.]**
In the Majorana normalization of EXACT.md, g=0 has a fast chain
`-2 sum e_A` with dispersion `8 |sin k|` and an entirely flat B chain.
The extensive zero-mode degeneracy prevents the complete low-energy
theory from being a finite-velocity compact CFT. Assigning only c=1/2
to the dispersive part does not describe the full chain.

Choose a finite-volume A ground state in a fixed boundary sector,
resolve its zero modes if necessary, and project the perturbation
onto it. In the thermodynamic limit its bond expectation is
`<e_A>=2/pi`, from the free-chain ground energy density `-4/pi`.
At first order the B Hamiltonian is therefore

```
H_B^(1)=-(2g/pi) sum e_B,
v_B^(1)=8 |g|/pi,        v_A^(0)=8.                       (1)
```

In C144's normalization the corresponding hopping is
`t=alpha+beta<a>`, with `<a> -> -2/pi`, alpha=0,beta=g/2.
Thus `t=-g/pi`, of order g, whereas the generated interactions start
at order g². The tuned cancellation curve instead has
`alpha approximately (2/pi) beta`; it is not this two-term line.

The C144/C157 second-order static kernel has a staggered 1/r tail and
a dimerizing sign. Its discussed wedge is controlled by competition
of a *cancelled* leading hopping with terms of order beta² log(1/beta).
The order-g hopping in (1) dominates that scale asymptotically.
This explains why their tuned dimer wedge cannot simply be assigned
to every perturbation of the flat point.

**[Heuristic continuum interpretation, not a uniform theorem.]** At
small fixed nonzero g the natural low-energy candidate is two moving
Majoranas, coupled through their energy densities, with total c=1
if this weak-coupling continuum description remains stable. Its bare
dimensionless interaction is of order `g/sqrt(v_A v_B)=O(sqrt(|g|))`.
This is useful power counting; it is not an all-order RG proof of a
gapless interval. In particular, finite-volume projection does not
give a remainder bound uniform in L: the A excitation gap vanishes
as 1/L, and static integration produces infrared logarithms. Nor do
we establish that two unequal Majorana velocities remain separately
observable after the interacting infrared flow.

The earlier [Nersesyan--Luther paper](https://arxiv.org/pdf/cond-mat/9404063),
section V, explicitly leaves the comparable-exchange/light-heavy
region outside its small-alternation analysis. Its asymmetric
Ashkin-Teller map is useful, but its name is not an exact solution
of the present line. We obtain no numerical g_* such that all
`0<|g|<g_*` are rigorously certified c=1.

The finite-size neutral estimator `v_N(L)=L Delta_N/(2pi)` is reported
separately in NUMERICS.md. At weak coupling its slow-bond spectral
weight and its value agree with (1). At g=1/2 or 1 it is a velocity
times an unmeasured scaling dimension, not a separately established
mode velocity. In particular the microscopic coefficient spreads
4 and 2 must not be confused with an inferred velocity ratio.

## A controlled large-|g| bulk exclusion

**[Proved duality; derived application of a read gap-stability theorem.]**
Use the real compass basis (1) of EXACT.md, and the ordinary Ising
duality

```
tau_j^x=Z_j Z_(j+1),     X_j=tau_(j-1)^z tau_j^z.
```

Since `Y_j Y_(j+1)=-X_j X_(j+1) Z_j Z_(j+1)`, the bulk dual is

```
H_dual = g sum tau_j^x
 -2 sum_(j even) tau_(j-1)^z tau_(j+1)^z
 -2 sum_(j odd)  tau_(j-1)^z tau_j^x tau_(j+1)^z.             (2)
```

On a ring the usual constraint `product tau^x=1` and the two boundary
twists corresponding to physical F must be retained. They change
global sector bookkeeping, not the bulk interaction. On an interval
there is an additional boundary spin; the displayed interior identity
does not omit a bulk degree of freedom.

For g>0, divide by g and add the irrelevant constant L. The unperturbed
part is `sum (1+tau^x)`, with a unique product vacuum, local gap 2,
and perturbations of norm 2/g and range three. Block pairs of dual
sites as `cell m=(2m-1,2m)` to make translation invariance literal.
The two perturbations centered at 2m and 2m+1 both lie within cells
m,m+1. Write the unperturbed density as half the sum of the onsite
Hamiltonians in these two cells: it has a unique product ground state
and local gap 1, and its translates sum to the original H0. Both pieces
now have the same bounded support. Each grouped perturbation has norm at
most 4/g. All assumptions of the small bounded-perturbation case of
[Yarotsky, Theorem 1](https://arxiv.org/pdf/math-ph/0412040) are then met
for sufficiently large g. There exists a finite, non-explicit G such
that the dual bulk remains gapped for g>G. Equivalently, choose the
theorem's allowed local norm epsilon for this cell convention and
take `g>4/epsilon`. No numerical value of epsilon is extracted here.

The dual constraint is compatible with the product vacuum for even
L. A boundary twist does not create an extensive gapless branch.
The physical Ising dual thus has a bulk gap above any broken-symmetry
ground-state partners. With the exact g-sign duality this supplies
an asymptotic exclusion of a positive-c critical bulk for |g|>G.
It does not certify g=4, locate a transition, or imply a positive gap
between the first two periodic physical eigenstates.

**[Numerical/perturbative ordering interpretation.]** The g->+infinity
states are the two Neel product states; g->-infinity gives the two
ferromagnets. The g=4 data exhibit their expected ordered continuation:
large staggered Z fluctuations, a shrinking symmetry partner splitting,
and a separate large excitation gap. The explicit rigorous result
claimed above is the asymptotic *bulk gap*. No quantitative lower
bound on the finite-g order parameter is derived here.

## Remaining interval

There is still no all-g c<=1 theorem or exact phase boundary. The
moderate rational points g=1/2 and 1 are numerically compatible with
c=1; g=2 is an especially simple equal-coefficient additional check.
Finite-size observations do not exclude a narrow additional phase
between the sampled points, a very small mass, or a long crossover.
