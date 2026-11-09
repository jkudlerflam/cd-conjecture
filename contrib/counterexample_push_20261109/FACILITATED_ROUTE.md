# A concrete modified qubit family to test next

OpenAI Codex, GPT-6 family, 2026-11-09. Proposal and calculations from
three parallel Codex agents and the coordinator. No human review or
counterexample claim.

## Why this route

The original C173 cusp has a bound pair. The exact binding criterion
suggests adding repulsion at effective distance two. A genuine
nearest-neighbour qubit term does generate it, but also generates
particle-number mixing. We retained both, computed their leading
effects, and tested the complete model on small rings.

The result is a concrete new search family, not just an adjustable
interaction in a longer-range microscopic Hamiltonian. It is a useful
next numerical attack alongside the more mature d=13 analytic route.

## Literal microscopic bond

Use Pauli matrices and repeat the same Hermitian bond on every link:

    h = Jx XX + Jy YY + Jz ZZ + Dz(XY-YX) + Dy(XZ-ZX)
        + (hz/2)(ZI+IZ).

Set Delta=1, lambda=r Delta, s=Delta^2/(2V), and choose

    w = -lambda^2/V + s(q r^2 + i),
    mu = (3Delta^2-2lambda^2)/V + s mu_eff,
    Jx=(Delta+Re w)/2, Jy=(Re w-Delta)/2,
    Jz=V/4, Dz=Im w/2, Dy=lambda/2, hz=(mu-V)/2.

The transverse DM term equals
lambda sum_j(n_j X_(j+1)-X_j n_(j+1)) on a periodic chain; the
single-site difference telescopes. Local dimension is exactly 2,
range is exactly nearest neighbour, and translation is by one site.

`microscopic_bond.py` writes all 4x4 matrices, in basis 00,01,10,11,
for V=40,80, r=0,.05,.10,.15 and two chemical-potential choices.
It checks Hermiticity and the occupation-basis identity up to the
explicit scalar and telescoping one-site term. The output
`microscopic_candidates.json` deliberately has no assigned central charge.

## Complete leading effective Hamiltonian

In units s, in the no-adjacent-occupation space P,

    H_r=P sum_j [ ((i+q r^2)b†_(j+1)b_j+h.c.)
                  -(b†_(j+2)b_j+h.c.)
                  +(-2+4r^2)n_j n_(j+2) -(2/3)n_j n_(j+3)
                  -mu_eff n_j
                  +2r(n_(j+2)-n_(j-2))X_j ] P.

**[proved local reduction; independent exact checks]** The last term
cannot be dropped. Its coefficient is order r, and its virtual energy
corrections are order r^2, the same order as the new repulsion. Three
independent enumerations agree on this operator. See
[the full reduction](../cusp_validation_20261109/FACILITATED_REPAIR.md)
and [the independent review](../route_audit_20261109/FACILITATED_REVIEW.md).

## Dressed branches and valley retuning

At r=0 the atomic minima are ka=-arcsin(1/4) and kb=pi+arcsin(1/4).
The density-assisted term couples an atom of momentum K to a pair at
separation two with amplitude 4ir sin K. The two valleys therefore
acquire unequal virtual shifts. This was checked directly on rings
with 9, 13 and 20 sites, then evaluated with two-body resolvents up to
1024 relative-coordinate sites.

**[numerical evaluation of second-order perturbation coefficients]**
At fixed mu_eff=-2.27, the shifts per r^2 are

    Sigma_a=-0.825224871396,
    Sigma_b=-2.663400053460,
    Sigma_pair=-1.163523967762.

The pair coefficient includes both the direct repulsion and its
virtual coupling to the interacting three-particle space. The latter
was computed through 57,121 gap configurations, with successive
cutoff differences below 4.1e-7. The pair-to-atom channel vanishes
exactly at its original momentum pi. The pair-to-trimer channel does not.

Choosing

    q=-(Sigma_a-Sigma_b)/(4 cos ka)=-0.474614791170

balances the two dressed valley energies to this order. The pair at
momentum pi has zero first-order response to this real hopping.
Consequently its binding relative to the balanced atomic threshold has
the formal expansion

    B(r)=0.0195265346664 - 2.32510095709 r^2 + O(r^4).

**[heuristic extrapolation]** Truncating this expansion predicts a zero
near r=0.0916. There is no bound on the r^4 remainder, so this is a
target for direct computation, not a proved unbinding point. At
mu_eff=-2.30 the corresponding retuning is q=-0.364132190385; the
retuning is source dependent, not a universal constant.

The fixed-source expansion treats isolated excitation branches. It
does not establish the thermodynamic vacuum, exclude larger clusters,
or determine finite-density scattering. The two valley momenta also
shift: their sum need not remain pi after dressing. The calculation
balances energies, not their momentum sum.

See [atom data](../cusp_binding_certificate_20261109/facilitated_atom_results.json)
and [pair derivation and data](../cusp_scattering_20261109/FACILITATED_EXTENSION.md).

## Direct finite-ring diagnostic

**[numerical; finite rings only]** The complete number-mixing H_r was
diagonalized in the exact constrained spin basis. The ground-state and
excitation calculations use periodic spin boundaries, with no
inappropriate fixed-fermion-parity twist.

At L=24, mu_eff=-2.30, q=-0.36413 (rounded for this run), the tracked pi-pair energy
minus the two valley-like excitation energies is

| r | E_pair-E_a-E_b |
|---:|---:|
| 0 | -0.041485743 |
| 0.10 | -0.019439073 |
| 0.15 | +0.002211339 |

This supports the sign of the perturbative trend. The positive margin
at r=.15 is small and strongly size dependent: a nearby retuning gave
about +.0285 on L=20 but only +.0026 on L=24. These levels are discrete
ring excitations, not an infinite-volume two-particle continuum. They
do not certify unbinding, a stable two-component liquid, or c>1.
Details and reproducible code are in the
[validation folder](../cusp_validation_20261109/README.md).

## Decisive next tests

1. Compute the dressed excitation and bound-state thresholds of the
   **complete** H_r at r around .08-.16, retuning the valley offset
   rather than assuming the old exchange symmetry survives.
2. Check atom-pair and pair-pair channels and larger-cluster onsets.
   Repulsion in one fixed-N matrix is insufficient when N mixes.
3. If a two-mode liquid survives at nonzero density, use converged
   size scaling and an independent estimator, then test the literal
   qubit bond at increasing V/Delta. Match bulk energies to catch
   metastable or symmetry-restored finite-size states.

The existing long DMRG runs for the original C173 candidate remain
useful controls. They do not test this modified family.
