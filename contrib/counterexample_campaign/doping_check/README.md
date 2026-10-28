# Exact squeezing for the one-vacancy lift

**Status: [proved finite-ring identities; independently checked numerical
implementation].** These results establish the charge-sea projection of a
literal local Hamiltonian. They do not establish a thermodynamic CFT at
fixed nonzero coupling.

This audit is independent of the cell-model implementation. It applies to
an arbitrary Hermitian two-site bond `h` on `C^q tensor C^q`, including
complex bonds with neither reflection nor time-reversal symmetry.

## 1. Literal model and the exact ordered-coordinate map

Add a vacancy `|0>` to the `q` internal states. On a ring of `L` sites set

    H = -t sum_(i,a) (|0,a><a,0| + |a,0><0,a|)_(i,i+1)
        + J sum_i h_tilde_(i,i+1),

where `h_tilde` is `h` on the two-occupied-site subspace and zero otherwise.
Fix `2 <= M < L` occupied sites. Write their ordered positions as
`x_1 < ... < x_M`, with internal labels `sigma_1,...,sigma_M`. This is a
unitary basis identification

    H_(L,M) = C^(choose(L,M)) tensor (C^q)^tensor M.

Let `R` rotate internal labels right:
`R|sigma_1,...,sigma_M> = |sigma_M,sigma_1,...,sigma_(M-1)>`.
An interior vacancy hop leaves the internal word unchanged. A hop from
physical site `L-1` to site `0` applies `R`; the reverse applies `R^dagger`.
Thus, writing `T_bulk` for the ordered-position hopping away from the cut
and `B` for the hop `L-1 -> 0`,

    H_t = T_bulk tensor I + B tensor R + B^dagger tensor R^dagger.

The coefficients `-t` are included in `T_bulk,B`. In a conventional
spinless-fermion Fock basis the cut hopping instead has the prefactor
`(-1)^(M-1) R`, because `c_0^dagger c_(L-1)` itself has the fermion
reordering sign `(-1)^(M-1)` in the fixed-M basis. Omitting this sign shifts
the flux and gives the wrong finite-ring sectors.

For `j<M`, let `A_j(x)=[x_(j+1)=x_j+1]`; let
`A_M(x)=[x_1=0 and x_M=L-1]`. Then the interaction is exactly

    V = sum_(j=1)^M A_j tensor h_(j,j+1),

with `j+1` interpreted cyclically and with the bond's original orientation
preserved. This identity needs no approximation and no special property
of `h`.

## 2. Exact projection onto the charge ground seas

Resolve the internal cyclic operator into eigenvalues
`R=exp(i theta_r)`, `theta_r=2 pi r/M`, `0<=r<M`. Let `Pi_r` be its spectral
projector. In this sector the charge particles are free spinless fermions
with momenta

    k_a = [2 pi n_a - theta_r - pi(M-1)]/L.

The ground sea consists of `M` consecutive integers
`n_a=ell,...,ell+M-1`. Its minimizing `ell` is `0` for `r<M/2` and `1` for
`r>M/2`; at `r=M/2` both seas must be retained. Let `P_sea` include every
one of these charge ground seas, tensored with the corresponding complete
internal cyclic eigenspace.

This is a chosen low-energy subspace of `H_t`. Its different cyclic
sectors generally have different kinetic energies of order `1/L`; it is
**not** one exactly degenerate eigenspace.

The physical translation eigenphase of a sea is

    exp[2 pi i (r-M ell)/L].

The integers `r-M ell` lie between `-M/2` and `M/2`. They are distinct,
including the two endpoint values when `M` is even. Since `M<L`, they are
also distinct modulo `L`. Therefore physical translation invariance of
`V` kills every matrix element between different retained sea branches.
This removes both cross-twist mixing and mixing between the two
half-flux ground seas.

Within a fixed branch, cyclic translation gives

    Pi_r h_(j,j+1) Pi_r = (1/M) Pi_r H_s Pi_r,
    H_s = sum_j h_(j,j+1).

The charge factor is consequently the expectation of the **total** number
of adjacent occupied pairs divided by `M`. Wick's theorem for a consecutive
free sea gives

    rho = M/L,
    |<c_i^dagger c_(i+1)>| = sin(pi M/L)/(L sin(pi/L)),

and hence

    alpha_(L,M)
      = rho - sin^2(pi M/L)/(M L sin^2(pi/L)).

The full projected operator is the direct sum, over retained branches, of

    P_sea V P_sea = alpha_(L,M) Pi_r H_s Pi_r.

All charge-sea branches have the same factor: their position-space
probabilities differ only by a common plane-wave phase. At fixed interior
density,

    alpha(rho) = rho - sin^2(pi rho)/(pi^2 rho) > 0.

The effective internal chain is therefore the original `h` chain
multiplied by `J alpha`, with its cyclic-momentum sectors glued to charge
boundary conditions. The construction costs exactly one added local state.

## 3. Two pitfalls avoided by the proof

Individual `E[A_j]` are not equal in the ordered-coordinate basis: the cut
singles out `j=M`. For example the saved `L=4,M=2` charge ground state has
weights `3/8` and `1/8`. The uniform coefficient arises after internal
cyclic-sector projection; it is not an independent uniform-bond
probability assertion.

Also, the charge level spacing vanishes as `L` grows. Fixed-ring
perturbation theory is controlled only with a coupling restriction that
shrinks with that spacing. The exact projected identity does not justify
integrating out a gapless charge fluid at fixed `J`, and an error
`O(J^2 L)` cannot control conformal gaps `O(1/L)`. A local effective-action
or explicit bosonization argument is still required.

For the proposed `q=23` application, another global issue remains: the
frozen internal chain prefers occupied lengths divisible by its cell size.
Choosing the occupied particle number can satisfy this on arbitrary
physical rings, but the retained charge vertex lattice and the marker
harmonics must be audited. The present finite-ring algebra does not prove
that all unwanted cell-length sectors are massive or that every physical
ring has a unique ground state.

## 4. Independent numerical checks

`check.py` constructs the full physical-ring matrices directly, constructs
both squeezed formulas independently, and compares them. It then builds
all cyclic projectors, all twisted charge ground seas (including exact
degeneracies), and checks the entire projected interaction against the
formula above. Every internal bond is a seeded random complex Hermitian
matrix, so the check exercises oriented and off-diagonal interactions.

Six cases use `q=2,3`, `L=4,5,6`, and `M=2,3,4`, with full dimensions from
24 to 270. The two exact basis identities have zero numerical error.
The largest projected-operator discrepancy is `3.14e-15`.
See [results.json](results.json).

Reproduce with Python plus NumPy:

    python3 contrib/counterexample_campaign/doping_check/check.py

Author and checker: GPT-6-family research agent, assigned independently
from the literal d=24 bond author, 2026-10-28. No human review is claimed.
