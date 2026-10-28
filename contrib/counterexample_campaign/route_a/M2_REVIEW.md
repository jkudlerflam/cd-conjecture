# Local elimination of the finite-cell machinery: primary-source audit

[verdict] The finite-order local matching required by route A is covered
by established perturbation theory. The remaining passage to the exact
infrared phase is conventional stability of a strict Luttinger-liquid
phase coupled weakly to massive modes. An exact all-order separation of
the extensive code band is unnecessary and is not claimed.

## What the literature actually supplies

[read] Bravyi, DiVincenzo and Loss allow arbitrary finite-dimensional
onsite ground spaces whose tensor product is the unperturbed ground
space. Their local Schrieffer-Wolff construction, Section 4.4,
Eqs. (4.23)-(4.32), gives an exact identity

    exp(T) H exp(-T) = H_block^[m] + R_[m].

At any fixed truncation order m, the discarded interaction has strength
O(|lambda|^(m+1)) uniformly in volume. Their strength is the maximum over
sites of the sum of norms of interactions containing that site. The
factor proportional to volume appears only when converting that bound
to a global operator norm. Lemma 4.2 controls the coefficient strengths;
the linked-cluster property controls supports. This construction does
not require the perturbed extensive bands to remain disjoint. Their
ground-energy error theorem is weaker than a CFT theorem, and should not
be used as one. [Primary paper, Sections 4.1-4.4](https://arxiv.org/pdf/1105.0675).

[read] Datta, Fernandez and Frohlich formulate the conjugation method
directly for volume-independent interactions, including exponentially
decaying tails. They separate this perturbative construction from the
additional assumptions needed to determine phases, and explicitly
identify renormalization-group analysis as another possible use of the
resulting interactions. Their strong-coupling Hubbard examples provide
the familiar setting with highly degenerate local retained spaces.
Their Pirogov-Sinai phase theorems are not invoked for our gapless
payload. [Primary paper, Sections 1.1 and 2](https://arxiv.org/pdf/math-ph/9809007).

## Application to the actual selected cells

[proved applicability] Fix n=23 and all cell parameters first. Within
one selected marker sector, block each 24-site cell into one site.
Subtract its vacuum energy and call its Hamiltonian H_0,j. Its ground
projector P_j has rank 2^23 and its gap Delta_0 is strictly positive.
The inter-cell operator V_j acts on two neighboring blocked sites and
has norm at most sqrt(3). The graph degree is two. Every local dimension,
gap and norm ratio is finite and independent of the number of cells.
Diagonalizing each H_0,j is an onsite basis change, so an onsite
"classical" unperturbed Hamiltonian is available if that formulation
is used. The physical data field is a bounded local term at order
lambda^2 and can be included in the same polynomial expansion.

[consequence] At second order, P V P=0 and the projected block is exactly

    P H_block^[2] P = lambda^2 (H_2 + epsilon N_data).

The exact conjugated full Hamiltonian differs by a local interaction
R_[2] with norm per blocked site at most C_2 |lambda|^3. Thus

    ||R_[2]||_local / lambda^2 <= C_2 |lambda| -> 0.       (1)

Choosing a higher fixed truncation order improves this power. The
constants can depend on the tiny cell and band scales; no uniformity
in n is required. Connected clusters occupy finitely many cells at
each order. The convergent tail for a fixed truncated generator has
rapidly decreasing interaction strength with cluster size for small
enough lambda. Mapping back to physical coordinates multiplies lengths
by only 24. The local generator also dresses local observables locally.

This is a stronger starting point than an extensive eigenvalue-error
bound: (1) controls the actual local perturbations of an exactly unitarily
equivalent model. The retained block need not have a many-body gap.
Its arbitrary logical degeneracy is an allowed input, not an omitted
exception to the perturbation theory.

## A uniform machinery threshold above the ground

[proved consequence; coordinating worker's refinement] Let D denote the
finite-order locally block-diagonal Hamiltonian H_block^[m], and write
D=H_0+W. Let Q_j=1-P_j and take a=Delta_0/2. The modified unperturbed
onsite terms

    H'_0,j = H_0,j-a Q_j

still have ground projector P_j and gap at least Delta_0/2. For a fixed
truncation order and sufficiently small lambda, the local block-diagonal
interaction W satisfies the strength condition of BDL Lemma 4.1 (or its
stated extension to a sum of bounded-locality orders) with this smaller
gap. Consequently D-a sum_j Q_j has a ground state in global P.
Its restriction to P is exactly P D P. If E_P is the bottom of P D P,
the variational characterization of its ground energy gives the operator
inequality

    D >= E_P I + (Delta_0/2) sum_j Q_j.                   (2)

This step does not require W to commute with individual Q_j. Since D
does commute with global P, (2) implies

    D restricted to (1-P) >= E_P + Delta_0/2.             (3)

The machinery therefore has a uniform positive threshold **above the
ground** in the exactly block-diagonal reference problem. High-energy
states within P may still overlap the spectrum in (1-P), so (3) neither
assumes nor proves separation of the entire extensive bands.
Equations (1)-(3) reduce the final task to a uniformly small local
perturbation of a stable critical block with a massive complement.

## What is still physical input

[standard effective-field-theory inference] The reference payload has
strict positive velocities, positive stiffness and nonzero low-order
momentum margins at fixed n=23. Its first remaining neutral non-Gaussian
operators are irrelevant. The finite-g nonquadratic error can be made
arbitrarily small by reducing g before choosing lambda. The local
remainder (1) can then be made smaller than this fixed stability
neighborhood. Virtual coupling to machinery with gap Delta_0 supplies
local derivative corrections; it does not add a gapless branch.
Keeping the Gaussian forward terms and integrating out the massive
fields is the usual effective-field-theory treatment.

This last statement assumes infrared stability under sufficiently small
symmetry-preserving local perturbations, including weak coupling to
massive auxiliary fields. Neither cited paper proves this multicomponent
CFT conclusion. Conversely, the missing result is no longer a special
uncontrolled finite-cell compiler or a locality assumption about an
extensive spectral band. The candidate has no identified mechanism for
a nonperturbative instability inside its strict Gaussian phase.

The separate vacancy step has its own gapless charge mode and must use
the material-coordinate analysis in `DEFECT_TRANSFER.md`; it must not
be called another application of a gapped global SW theorem.

Source audit and application by Codex GPT-6-family route-A worker,
2026-10-28. No human review is claimed.
