# Read-only audit of the two-flavour demonstrator selector

Reviewed 2026-10-28 by Codex `demonstrator_ir`, independently of the original
author. Files reviewed: [`README.md`](../../routeA_demonstrator_n2/README.md),
[`selector2.py`](../../routeA_demonstrator_n2/selector2.py),
[`cell2.py`](../../routeA_demonstrator_n2/cell2.py), and saved
[`selector2.json`](../../routeA_demonstrator_n2/data/selector2.json).
No files in that contribution were modified.

## Verdict

[proved] The sum of run ground energies is an operator lower bound, even
though neighbouring runs share B logical factors. It is not generally
saturated in the PBP sector. Consequently comparing u_k against u_2 does
not establish global sector selection. The asserted positive per-site
margin against every other pattern cannot hold for patterns with finitely
many defects in an arbitrarily large ring.

[proved] There is a simple independent selected-sector variational upper
bound. It separates the selected sector from the reported entirely wrong
run patterns, but does not yet control arbitrarily dilute defects.

Two additional problems are concrete: the finite tail search omits its
infinite-length limit, and the saved lambda=0 spectra for k=3,4 disagree
with explicit zero-energy product states. Neither issue by itself refutes
the proposed phase.

## Lower bounds cannot be compared as though they were ground energies

For a type word with N_k runs of length k, the established inequality is

    E(word) >= sum_k N_k (-mu+e_k).

Writing theta for the fraction of physical sites in runs with k != 2 yields

    E(word)/L >= (1-theta) u_2 + theta u_other,

provided u_other bounds all other runs, including the infinite-run limit.
This is useful only relative to an upper bound U_selected on the actual
selected ground energy density. If u_2 < U_selected, then it excludes
only

    theta > (U_selected-u_2)/(u_other-u_2).

It says nothing about the remaining dilute-defect region. Sharing B
factors causes no problem for the lower inequality; it is the attempted
saturation of that inequality that requires another argument.

The impossibility of a fixed per-site margin against every different type
word is even simpler. Change a selected word in a bounded interval, for
example replacing two consecutive k=2 runs by k=1,k=3 while preserving total
length and B count. Start from a selected ground state, keep its reduced
state outside that interval, and insert any state of the new types inside.
Only a fixed number of bounded nearest-neighbour terms change. The
competing sector therefore has a variational state of energy at most
E_selected+C, where C is independent of L. Its possible per-site energy
excess is at most C/L. A positive cost per defect would be consistent with
this, but is a different statement that needs a local comparison.

## An explicit selected-sector trial state

[proved] Pair adjacent PBP cells into dimers and diagonalize their one
internal P-P bond. Keep both unused outer P states in their exact Bell
states. The outer checks annihilate those Bell states for every B logical
state, so this is exactly the k=2 run problem, with ground energy e_2.
Take the product of these dressed dimer ground states. Every inter-dimer
port has maximally mixed x marginal, so the expectation of every traceless
port operator is zero and the connecting P-P interaction contributes zero.
The chemical energy per two-cell dimer is -2 mu. Thus on rings with an
even number of cells, and hence in the selected thermodynamic sector,

    U_selected = -mu/3 + e_2/6.

This uses actual run eigenvectors as trial states; it does not require
simultaneous minimization of neighbouring runs. Negative e_2 can also be
proved without a diagonalizer: the code has zero energy, the port bond
creates a nonzero orthogonal vector, and the Rayleigh quotient of the code
minus a sufficiently small multiple of that vector is negative.

Using the saved e_2 values gives the following illustrative bounds. The
numbers remain numerical until the finite eigenvalue bounds are certified.

| lambda | selected trial upper | other-run lower | excluded theta above |
|---|---:|---:|---:|
| 0.10 | -0.4027588649 | -0.3694152876 | 0.0764 |
| 0.15 | -0.4061768565 | -0.3694152876 | 0.1439 |
| 0.20 | -0.4109472479 | -0.4000000000 | 0.5000 |

In particular, the trial bound is strictly below the all-P bound -2 lambda
even at lambda=0.2. This repairs the pure-competitor comparison but not the
claim for every pattern.

## The tail needs its endpoint at infinity

The declared tail expression is

    F(k) = (-mu+2 f+vmin lambda(k-3))/(k+1)
         = vmin lambda + (-mu+2f-4 vmin lambda)/(k+1).

It is monotone for k>=K, with K=kmax+1. Therefore its exact infimum is

    min(F(K), vmin lambda).

The current loop ends at k=2999. At lambda=0.2 the coefficient of 1/(k+1)
is positive, so the true infimum is -0.4, whereas the saved tail is
-0.399982133110834. The latter is not a lower bound valid for every length.
The numerical difference is small, but a rigorous all-length statement
must use the limiting endpoint.

## The zero-coupling spectral row fails an exact witness

[proved] At lambda=0 and every k>=2, put Phi_R on the first P, Phi_L on the
last P, arbitrary states on the intermediate P's, and arbitrary states on
the two B factors. Both endpoint checks annihilate the product; the checks
are positive. Hence e_k=0 exactly for all k>=2.

The saved data instead give e_3=0.2917237469701791 and
e_4=0.291723746970178, equal to the one-cell gap. The correct zero-coupling
values are u_3=-0.3 and u_4=-0.24. This discrepancy is consistent with a
missed zero-energy eigenspace or stale outputs, but the present audit does
not establish its software cause. Positive-lambda run calculations should
be independently checked with residuals and explicit variational witnesses
before treating the saved eigensolver results as certified ground energies.

The finite-sector effective-bond calculations and the reported c≈2
numerics are separate claims; this audit does not refute them.
