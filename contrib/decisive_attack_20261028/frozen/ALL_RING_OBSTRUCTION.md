# A genuine all-ring obstruction for the frozen construction

## Statement [proved]

Under the same exact microscopic hypotheses as
[the defect-pressure theorem](DEFECT_PRESSURE.md), define

    delta_* = min(gamma-3B, K-4-5B)>0,
    C = 2B + 2b ||h||.

For every sufficiently large ring length `L` not divisible by the selected
cell length `b`, the frozen microscopic ground-space dimension obeys

    dim ker(H_L-E_0(L)) >= ceil(L delta_*/C).               (1)

In particular, the complete finite-size spectrum on all ring lengths cannot
have a scaling limit with the finite vacuum multiplicity of a fixed compact
unitary CFT. This statement is compatible with the exponential suppression
of defect contributions to the thermodynamic free energy proved separately.
The two statements use different orders of limits.

## 1. Uniformly bounded defect count

Let `u<0` be the exact homogeneous selected-cell ground energy per cell.
The existing boundary estimate gives

    |A_m-mu|<=B,   S_m<=mu+2B,

for open and periodic selected-cell chains, respectively. Also
`u<=E_b<0`. Fix a conserved marker/illegal-arrow sector. Let `r` be its
number of illegal bonds and `k` its number of complete wrong-length cells,
using the end-fragment convention in DEFECT_PRESSURE. The exact sector
bounds from the original global-sector argument give

    E_sector >= (u/b)L + k(gamma-3B) + r(K-4-5B)
             >= (u/b)L + delta_* D,    D=k+r.              (2)

For legal marker-free sectors, the energy is instead nonnegative and
`D=0`; they will be excluded below.

Write `L=mb+s`, `0<=s<b`, taking `m>=3`. On the full Hilbert space of the
`L` sites, compare `H_L` with `H_(mb) tensor identity` on the last `s`
sites. When `s>0`, the comparison removes one closing bond and inserts
`s+1` bonds, so the operator norm difference is at most `(s+2)||h||`.
The variational principle and a ground state in the selected sector of
`H_(mb)` give

    E_0(L) <= mu+2B+(s+2)||h||
            <= (u/b)L + 2B + (2s+2)||h||
            <= (u/b)L + C.                               (3)

Here `|u|/b<=||h||`, because the energy density of any restriction of the
microscopic chain lies between `-||h||` and `||h||`. For `s=0`, the selected
periodic-sector upper bound is already sufficient. Since `u<0`, the right
side of (3) is negative at sufficiently large `L`, excluding legal
marker-free ground sectors.

Equations (2)-(3) imply that every ground sector has

    D <= C/delta_*.                                      (4)

If `b` does not divide `L`, such a sector also has `D>=1`: a legal ring
containing markers with no wrong cell has exactly length `mb`. Illegal
sectors have `r>=1` by definition.

## 2. Translation-orbit multiplicity

The marker projectors and illegal-bond projectors commute with one another
and with `H_L`. Choose a normalized ground state in one of their common
sectors. Let `q` be the length of the orbit of that sector under one-site
translation. Then `q` divides `L`, and its marker/illegal pattern consists
of `L/q` identical periods. The integer defect count is translation
covariant and additive across these repetitions, so

    D = (L/q) D_period,

with positive integer `D_period`. Thus

    q >= L/D >= L delta_*/C.                             (5)

Translated states in distinct classical sectors are orthogonal, while
translation invariance of the Hamiltonian makes all their energies equal.
They provide at least `q` ground states. Combining (4)-(5) proves (1).
The argument uses the joint marker/illegal pattern, not only the marker
count, so it also covers ground sectors containing illegal interfaces.

As a simpler corollary, on a prime ring any nonconstant conserved binary
marker pattern has an orbit of exactly `L`. The general argument above
avoids restricting to prime lengths or proving that every ground sector
has that particular type of pattern.

## 3. Which infrared statements are excluded

A full finite-size scaling partition function at a fixed aspect ratio
`beta/L` always contains its exact ground-space multiplicity. After only
the usual ground-energy subtraction it therefore diverges along the
incompatible lengths in (1). It cannot converge to the partition function
of one fixed compact CFT, or a fixed finite direct sum of such CFTs, with
finite vacuum multiplicity.

On the other hand, at fixed `beta` the contribution `log L/L` vanishes
as `L` tends to infinity. The defect-pressure theorem proves more: all
thermally populated defects change the bulk pressure only exponentially
in inverse temperature. The full chain can therefore have the same bulk
thermal CFT coefficient as a selected homogeneous sector even while its
incompatible-ring vacuum multiplicity is not conformal.

The original conjecture [Latorre and Sierra, Section 2](https://arxiv.org/html/2403.17242v2#S0.SS2)
requires a nearest-neighbour one-site-invariant Hamiltonian and a low-energy
spectrum described by unitary CFTs. It does not explicitly require a unique
ground state or a common finite-size limit on every ring-length sequence.
The paper does not settle this distinction. Accordingly, the frozen
construction is a candidate for a bulk formulation, and is excluded from
the stronger all-ring finite-size formulation just described. Resolving
that interpretation is separate from proving its 23-mode payload theory.

Provenance: Codex GPT-6-family subagent, 2026-10-28. The coordinating agent
checked the bounded-defect and translation-orbit logic during development.
No human review is claimed.
