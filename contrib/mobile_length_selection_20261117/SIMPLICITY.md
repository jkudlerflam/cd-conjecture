# Simplicity assessment

2026-11-17. This is a test of the existing model, not a newly simplified
counterexample. The parent bond and earlier contributions are unchanged.

## Independent parameters and literal processes

In the **restricted legal vacuum block**, after fixing N and dividing
by t, the general physical controls are q, J/t and mu_sel. The proposed
rule mu_sel=-(n+2)E_n(q) leaves two independent continuous controls at
fixed target n. The positive length tilt s/t is an ensemble parameter,
not an extra microscopic interaction. The occupied chemical potential
is a scalar at fixed N. K and charged-rail coefficients vanish in this
block; they are not removable from the full physical Hamiltonian.

Two independent controls do not mean two operator monomials. The
vacuum moves have tied diagonal entries 1,q^2 and off-diagonal -q.
They coexist with vacancy hopping, the terminal projector, and the
marker potential. At fixed all-zero payload there are two undirected
vacuum moves and five distinct material-vacancy swaps, in addition to
diagonal matrix units. Retaining all payloads gives four vacuum moves
and eight material-vacancy swaps in the vacuum restriction. These are
counts in the parent's literal site basis, not a claim of minimality
under basis changes. Grouping them as clock and hopping sums does not
remove processes.

The full d=22 family has the much larger process inventory in the
[previous audit](../mobile_vacuum_clock_20261117/PARAMETERS.md).
Switching lambda=f=0 removes the four port transitions and the readout
operator for this diagnostic. It also switches off the proposed
dispersing payload mechanism. Restricting to vacuum flags does not
delete all the charged terms from the full model.

## Coefficient hierarchy and tuning

After removing a common positive scale, the restricted literal scales
include t=1, J, Jq, Jq^2, Jq^2/2 and J|mu_sel|/2. The last two arise
from actual onsite contact halves. The marker term cannot be removed
as a constant, even with fixed N and marker count, because the contacts
fluctuate. A material-only constant or telescoping bond term becomes
contact weighted after the mobile lift; subtracting it without its
vacancy compensation changes the physical Hamiltonian.

At q=1/4 the tied vacuum coefficients have spread 16. The allowance
J/t=1/80 introduces another separation. For n=1 the marker-half scale
is about 6.53e-5 t, so the literal nonzero entries include a spread
of at least about 15,300. For n=2 this lower bound on the spread is
about 186,000. The selector scale shrinks geometrically with n.
These are microscopic small coefficients, not merely emergent slow
velocities. Additive contributions to other diagonal entries do not
erase the displayed marker-half entries.

The exact analytic prescription for mu is not evidence that every
decimal must be tuned. The frozen distance to this weak competitor is
Delta_n=(n+1)E_(n+1), roughly a q^2 fraction of |mu_sel|. At q=1/4
that fraction is about 4.1% for n=1 and approaches 6.25% as n grows.
This distinguishes an exponentially small absolute coefficient from
an exponentially small fractional tuning tolerance. The corresponding
uniform **mobile** tolerance is not known. The bounded n=1 check moves
the crossing away from the chosen mu rather than requiring retuning.

## What this contribution changes

| Result or choice | Simplicity consequence | Limit |
| --- | --- | --- |
| Tail-sensitive scalar reduction | Adds no states, penalties or couplings; retains J/t<=1/80 independent of length and D | Uniform positivity of the relative scalar bound remains open |
| Favorable clock correction at second order | Avoids treating an O(E_n) term as an arbitrary adverse error against E_(n+1) | Does not control all higher orders |
| q=1/4 tests | Do not respond to a conservative estimate by making q smaller | Not an all-n band or CFT certificate |
| Fixed mu_sel | Works in every bounded test; the crossing shift is quantified | No uniform retuning theorem or mobile stability interval |
| Real vacuum restriction | Exact cyclic ground-energy gluing removes the global boundary allowance here | Does not apply automatically to the complex charged/port model |

There is no reason from this calculation alone to introduce a new
penalty or an n-dependent shrinkage of J/t. Conversely, the partial
result does not make the full architecture a moderate-ratio model.
It should remain an existence-construction research direction rather
than outrank a comparably supported model with fewer physical terms
and robust order-one coefficients merely because d is fixed.
