# A one-vacancy mobile completion of the 23-state candidate

Historical construction note. The completed physics-level argument and
independent reviews are in [COUNTEREXAMPLE_D24.md](../COUNTEREXAMPLE_D24.md).
The conditional wording below records the initial stage of this work.

Initial status: explicit microscopic candidate; conditional infrared conclusion.
The proposed values are local dimension d=24 and total Virasoro central
charge c=24. The latter has not been measured or established for the
microscopic Hamiltonian. This note separates the exact construction from
the infrared transfer under investigation.

## 1. One repeated 576 by 576 bond

Let h be the literal 23-state bond of `../../explicit_candidate/bond23.py`,
with the n=23 finite-loop parameters of `../../../codex/work_C33/PROOF.md`
Section 9 and the sector bounds in `../route_a/`. Define

    H_site = C|vac> direct-sum C^23,
    n = I - |vac><vac|,
    h_occ = (iota tensor iota) h (iota tensor iota)^dagger,
    h_24 = -t sum_{a=1}^{23} (|vac,a><a,vac| + h.c.)
           + J h_occ - (mu/2)(n tensor I + I tensor n).

Here iota embeds occupied labels, t>0, J>0, and all coefficients are
fixed before the physical ring length tends to infinity. The sum of
h_24 over bonds is Hermitian, nearest-neighbour and one-site translation
invariant. No marker pattern or particle-number sector is imposed by
the Hamiltonian. Occupied number, original total charge and original
marker number are exactly conserved. The local dimension is exactly 24.

The entire original h, including the half-site terms, is lifted into the
occupied-pair block. Moving the original on-site terms outside that block
would change their first-order relative coefficient. The occupied chemical
potential mu is a separate physical parameter.

## 2. Exact squeezing and the first-order coefficient

For N occupied sites x_1<...<x_N, retain their labels in the same order.
A hop not crossing the physical ring cut changes only the positions.
A hop across the cut cyclically permutes the label list. Thus the zero-J
problem is the familiar hard-core charge problem with a boundary twist
given by the internal cyclic-translation eigenvalue. Mapping hard-core
particles to fermions also contributes (-1)^(N-1) to that boundary twist.
This is an exact coordinate change, not an approximation.

An occupied bond acts on consecutive squeezed labels only when their
physical positions are adjacent. At finite ring size the adjacency
probabilities of individual squeezed bonds need not be equal because
of the chosen cut. Within a definite internal cyclic-momentum sector,
all projected bond operators have the same matrix elements, so their
sum is the squeezed Hamiltonian times the mean number of adjacent
physical occupied pairs divided by N. Cross-twist sectors and finite-size
charge excitations must be retained for a complete low-energy theory.
See `../doping_check/` for the independent operator checks.

At thermodynamic charge density rho=N/L, the free spinless charge sea has

    <n_i n_(i+1)> = rho^2 - [sin(pi rho)/pi]^2,
    alpha(rho) = rho - sin(pi rho)^2/(pi^2 rho) > 0.

Consequently the leading internal Hamiltonian is J alpha(rho) sum_j h_j.
This preserves every relative internal coupling. The leading charge
sector has c=1 and velocity 2t sin(pi rho). At rho=1/2,
alpha=1/2-2/pi^2. Positivity follows from |sin x|<x for 0<x<pi.

This is a first-order and thermodynamic-reference statement. The charge
sector is gapless. Its finite-ring level spacing is not a denominator
that can be used uniformly when L grows.

## 3. Why this could complete the counterexample

The internal n=23 candidate has a period-b=24 marker pattern in occupied
label coordinates. Physical particle positions fluctuate, so this is
not a frozen physical crystal. Its b translates are tied to charge
boundary twists by cyclic label relabeling. Their finite-ring energies
belong to charge zero-mode sectors rather than an independent choice of
one frozen defect position on every incompatible physical ring.

For each physical length L, the occupied number can be chosen to be a
multiple of b within a bounded distance of rho L. The corresponding
charge density change is O(1/L). Internal wrong-length defects therefore
need not be present on prime physical rings. Showing that the complete
finite-J ground state actually selects those sectors remains part of
the transfer problem.

Below the internal defect scale the low-energy charge compactification
must respect this b-particle selection rule. Fractional charge harmonics
associated with the internal marker pattern cannot be silently discarded.
The independent infrared analysis in `../infrared/` tracks them.

## 4. Required finite-coupling argument

First fix the internal n=23 parameters inside the strict band and
interaction stability margins. Then take J/t sufficiently small but
strictly positive. The desired statement is that the fast charge sector
induces corrections to the internal Hamiltonian that are o(J) in a local
low-frequency norm, with any gapless mixed fields kept explicitly. The
remaining charge and internal sectors must have a positive joint Gaussian
form and no relevant non-Gaussian uniform vertex. All comparisons must
use the smallest internal velocity, however tiny, not only its bandwidth.

If these assertions and the internal c=23 inference hold, the complete
chain has c=23+1=24, exceeding d-1=23. A finite-size degenerate-perturbation
calculation alone does not establish them. We do not call this a completed
counterexample until the outstanding assertions have been justified at
the repository's theoretical-physics standard.

## Sources and provenance

[read, abstract] Ogata and Shiba, Phys. Rev. B 41, 2326 (1990),
https://doi.org/10.1103/PhysRevB.41.2326, derives the charge/spin
factorization for the large-U Hubbard chain. [read, abstract] Xiang and
d'Ambrumenil, Phys. Rev. B 45, 8150 (1992),
https://doi.org/10.1103/PhysRevB.45.8150, studies a separated squeezed
Heisenberg wavefunction for the t-J model at small J/t. These are
precedents, not proofs for the arbitrary 23-label h above.

The construction was independently suggested by the coordinator and the
route-A worker on 2026-10-28. All work is AI-assisted, GPT-6 family; no
human review claimed.
