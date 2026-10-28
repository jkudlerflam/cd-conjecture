# Independent source review: local SW and the final d=24 argument

Reviewed by OpenAI Codex, GPT-6 family, 2026-10-28. Sources opened:
[Bravyi, DiVincenzo and Loss, arXiv:1105.0675](https://arxiv.org/html/1105.0675),
sections 4.1, 4.2 and 4.4, including equations (4.18)-(4.32), and the
current `COUNTEREXAMPLE_D24.md`. No human review is claimed.

**Verdict.** The source materially strengthens the internal matching
step. At fixed truncation order it supplies an exact change of basis
with a small *local* remainder on the full Hilbert space. An
all-order exactly block-diagonal effective action is unnecessary for
the conventional phase-stability argument. Together with the finite
sector and band checks and the phase-model vacancy analysis, I
would accept the construction as complete at the repository's
stated theoretical-physics standard. The result remains distinct
from a rigorous lattice-to-CFT theorem or a measured c=24.

## 1. What is actually supplied by the source

[read] BDL allows finite-dimensional sites with arbitrary degenerate
on-site ground subspaces and a positive on-site excitation gap.
The perturbation consists of bounded local terms. For fixed order
m, its local SW construction uses a finite polynomial anti-Hermitian
generator. Equations (4.23)-(4.25) give an exact transformed
Hamiltonian: a locally block-diagonal truncation plus a remainder.
Lemma 4.2 bounds every remainder coefficient in the interaction
strength norm of equation (4.3), independently of the number of
sites. The system-size factor enters when converting that norm to
the global operator norm. The block-diagonal stability lemma and
Corollary 4.1 separately control the truncated reference ground
sector. None of these results asserts an all-order thermodynamic
block diagonalization or CFT correlation theorem.

## 2. Application to this particular frozen chain

[proved application] First restrict to the selected, legal,
period-24 material sector, justified separately by the campaign's
sector audit. Block each complete cell into one site. Its Hilbert
space is finite, its code is exactly degenerate, and the minimum
gap delta to other cell states is positive. At fixed nonzero g,
the unperturbed Hamiltonian is a sum of these independent cell
Hamiltonians. Port couplings are nearest-cell interactions.

The logical data field is of order lambda^2. It can either be
included directly as an analytic second-order perturbation or
written lambda[V_1+lambda V_2]. In the latter form V depends on
lambda, but its local norm is uniformly bounded for
0<=lambda<=lambda_0. This does not alter the fixed-order bounds.
The on-site norm/gap ratio is large but finite at n=23 and must be
included in constants; it is independent of the number of cells.

Write the exact fixed-order transformed Hamiltonian as

    U_m H U_m^dagger = H_bd^(m)+R_m.                              (1)

For the specified decomposition of the remainder, its local norm
obeys, by a geometric sum of the coefficient bounds,

    ||R_m||_loc <= C_m |lambda|^(m+1).                            (2)

This is an identity plus a full-space norm bound. It is stronger
than merely knowing that several formal effective coefficients
exist, or that an extensive ground-energy error is small per site.
It neither uses nor implies isolation of the whole extensive
logical band.

[proved spatial refinement] Starting from nearest-cell terms, every
nonzero nested commutator has a connected declared support. At
coefficient order r its diameter is at most r in cell units.
Therefore inserting a fixed weight exp(mu diam) multiplies the
coefficient bound by at most exp(mu r). Reducing the allowed
lambda by the corresponding fixed factor preserves the geometric
sum, giving

    ||R_m||_mu <= C_(m,mu) |lambda|^(m+1).                        (3)

The unitary itself is generated for finite time by a bounded,
finite-range local Hamiltonian. Its quasi-local dressing preserves
the thermodynamic critical-mode count. Exact charges are retained
because the local construction is built from the symmetric cell
Hamiltonians, projectors and perturbations.

## 3. A useful ground-threshold strengthening

[proved deduction] Let Q_j project out of the code on cell j and
N_def=sum_j Q_j. Apply block-diagonal ground-sector stability to

    H_bd^(m) - (delta/2) N_def.

Its unperturbed on-site gap is at least delta/2. For sufficiently
small lambda, its ground state lies in the code product space.
Subtracting the displayed defect term does not change that code
Hamiltonian or its ground energy E_eff. Consequently

    H_bd^(m) >= E_eff I + (delta/2) N_def.                         (4)

This provides a genuine massive threshold above the ground for
the block-diagonal reference. It does not claim that every
high-energy logical state lies below every auxiliary excitation.
Thus it is compatible with the extensive-band overlap emphasized
by C36 and C40.

## 4. Why this is enough for the physics inference

[controlled physics inference] The target payload scale at fixed
g is W=lambda^2 w_g with w_g>0. Taking m=2 already gives

    ||R_2||_mu/W <= (C/w_g)|lambda| -> 0.                          (5)

One can instead retain a higher fixed order. The result is an
exact local-unitary equivalence to a reference with a massive
auxiliary sector, a low-energy Hamiltonian arbitrarily close to
the certified 23-band Gaussian basin, and an arbitrarily small
quasi-local perturbation relative to its slowest retained scale.
Ordinary stability of that Gaussian-plus-massive phase now applies.

The small-g choice still comes first: the order-lambda^2 g^4
non-Gaussian interaction must be small relative to the order-
lambda^2 g^2 payload. Reducing lambda alone does not do that.
After g is fixed, (5) controls the additional local matching error.
The residue R_m need not be projected onto the code using an
extensive spectral gap, and it need not be removed exactly.

The continuum statement remaining here is conventional phase
stability under small symmetry-preserving quasi-local perturbations,
including small coupling to a massive sector. That is the same
bosonization/RG input already used for the 23-band interacting
payload. BDL does not prove this final continuum statement, but
it removes the need to assume C40's stronger exact all-order
low-frequency kernel merely to establish microscopic local
smallness. Requiring such an exact kernel is sufficient, not
necessary, for the physics-standard construction.

## 5. Audit of the final synthesis

The inspected `COUNTEREXAMPLE_D24.md` is accurate about the
Hamiltonian, finite root certificate, phase-model coordinate
change, harmonic lattice, order of limits, and the difference
between a physics argument and a microscopic c measurement.
Its final distinction from a rigorous mathematical CFT theorem
should remain.

Its paragraph identifying M2 as the remaining internal matching
input can now be sharpened. Suggested replacement:

> A fixed-order local Schrieffer-Wolff transformation gives an
> exact quasi-local change of basis with a full-space local
> remainder O(lambda^(m+1)). At fixed g this remainder is
> arbitrarily small relative to the lambda^2 payload scale.
> Block-diagonal stability keeps the auxiliary sector massive.
> The final step is the stated conventional Gaussian phase
> stability argument; an all-order exact block diagonalization
> or a mathematical CFT-limit theorem is not claimed.

This review concerns the selected material sector. It does not
replace the independent all-length/illegal-sector comparison,
the vacancy defect dressing analysis, or the C33 matrix-inertia
certificate. Those checks are linked by the synthesis and must
continue to accompany its conclusion.
