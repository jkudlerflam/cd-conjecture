# C98. The bound from the Casimir energy

2026-10-29. Completed investigation. Proofs and qualifications are in [PROOF.md](PROOF.md); numerical checks are in [CHECKS.md](CHECKS.md) and [RESPONSE_CHECK.md](RESPONSE_CHECK.md).

## Verdict

[Open] I did not prove or falsify

    A <= pi M log_2 d,
    E_0(L)-e_0 L >= -pi^2 M log_2(d)/(6L)+o(1/L).

The logarithmic bound survives the specified solvable-chain tests. Its proof for the full symmetric SU(N) family and critical XXZ was already established conditionally in C35, so this work does not count that as a new result. The stronger statement with the same coefficient and no remainder at every finite L is false, including on the Heisenberg chain.

[Proved] The main new reduction is an exact representation of the doubling defect as an integrated response to two separated bond reconnections. It supplies a concrete sufficient extra hypothesis, stated below. That hypothesis has not been derived from locality or checked in the thermodynamic limit for the solvable families. No weakest noncircular physical hypothesis has been identified.

## What the elementary routes establish

[Proved] With D_L=E_0(2L)-2E_0(L), the vacuum Casimir expansion gives

    D_L=pi A/(4L)+o(1/L).

Thus the desired comparison is

    D_L <= pi^2 M log_2(d)/(4L)+o(1/L).                (1)

Product gluing instead gives D_L<=4M. More directly, tiling by L-site ring ground states gives the unconditional boundary estimate

    |E_0(L)-e_0 L|<=2M.

For XX, the product-gluing cost tends to 2/pi in the normalization h=SxSx+SySy, while D_L~pi/(4L). Moving or averaging the cut does not change this total cost.

[Proved] Open-chain closure also has a boundary term. On even XX chains,

    E_ring(L)=-csc(pi/L),
    E_open(L)=1/2-[2 sin(pi/(2(L+1)))]^-1,
    E_ring(L)-E_open(L) -> -(1/2-1/pi).

One must subtract the boundary energy before asking for an inverse-length estimate. Endpoint dimension does not control that subtracted energy.

[Proved] The Schmidt rank across a half-ring cut is 2^(L/2) already for XX. A d-dimensional endpoint is not a bound on the many-body Schmidt rank. A useful MPS or purification construction would need an average nearest-neighbor marginal error of order log(d)/L^2. Section 5 of PROOF.md gives the sufficient condition and all constants, but does not establish such a construction.

[Proved] A smooth product of local partial-copy swaps is not a low-cost twist. On the Heisenberg ground state, an intermediate swap angle theta adds

    0.704986412078... sin^2(2theta)

per doubled bond. The analogous XX coefficient is 0.447316023917.... A smooth angle ramp therefore has an extensive energy cost. The missing property is continuous symmetry under partial copy rotations; discrete exchange of the copies is insufficient.

[Proved logical obstruction] Reflection concavity and an O(M) boundary ceiling do not bound the Casimir coefficient. The abstract energy sequence

    E_L=e_0 L-B xi/(L+xi)

has a bounded positive deficit, increasing energy density, concavity, and positive doubling defects, while its 1/L coefficient is B xi. This is a counterexample to that inference, not a Hamiltonian counterexample. The actual [Lieb-Schupp comparison](https://arxiv.org/html/1410.5322) also requires special crossing interactions; arbitrary Hermitian bonds need not satisfy it.

## A precise remaining response estimate

[Proved] Start with two L-site rings. At each of two separated cuts, V_x or V_y replaces two straight bonds by two crossed bonds. Let H(s,t)=H_00+sV_x+tV_y. At the corners, switching one cut produces one 2L-site ring; switching both returns two L-site rings up to relabeling. For the Gibbs state of H(s,t), define the connected Kubo susceptibility

    chi_xy^(beta)(s,t)
      =Re integral_0^beta d tau <delta V_x(tau) delta V_y>.

Then, exactly at each finite L,

    E_0(2L)-2E_0(L)
      =(1/2) lim_(beta->infinity) integral_0^1 ds integral_0^1 dt chi_xy^(beta).    (2)

The finite-temperature derivation handles level crossings without differentiating a ground state. V_x and V_y themselves are invariant under scalar shifts and bond coboundaries, and their norms are at most 4M.

[Proved conditional implication] A sufficient additional hypothesis is

    limsup_(L->infinity) L limsup_(beta->infinity)
      integral_[0,1]^2 [chi_xy^(beta)]_+
        <= (pi^2/2) M log_2 d.                       (3)

Equations (2) and the vacuum expansion then prove A<=pi M log_2 d. More generally, replacing pi^2/2 by B proves A<=(2B/pi)M log_2 d. The positive-part condition permits negative response and does not impose pointwise decay throughout the interpolation.

[Open] The elementary estimate is only |chi_xy^(beta)|<=16 beta M^2, which is useless as beta tends to infinity. The zero-temperature response contains inverse excitation energies. Local dimension bounds entropy but does not control those denominators. The signed integral version of (3) would just restate the desired bound and should not be presented as a proof. An independent bound on the positive-part response, or a controlled variational extension, is the missing input.

## Checks

[Assumed continuum data; proved comparison inherited from C35] The asymptotic normalized ratios are

| Chain | A/(pi M log_2 d) |
|:--|--:|
| XX | 2/pi = 0.636619772 |
| XXZ, Delta=cos gamma | sin(gamma)/gamma |
| Heisenberg | 1 |
| Fundamental SU(3) | 0.841239671 |
| Fundamental SU(4) | 0.75 |
| Symmetric SU(2), k=2 | 0.630929754 |
| Symmetric SU(3), k=2 | 0.550190659 |
| Symmetric SU(7), k=2 | 0.422632833 |
| Symmetric SU(8), k=2 | 0.406195448 |

[Proved finite-size failure] At Heisenberg L=4,

    e_0 L-E_0(L)=3-4 log 2=0.227411277760...,
    pi^2 M log_2(d)/(6L)=0.205616758356....

This is a counterexample to the sharp finite-L statement without a remainder. It is compatible with saturation of the asymptotic bound.

[Numerical] New ED covers XXZ at Delta=-1/2,0,1/2,1 through L=12, fundamental SU(3) through L=9, and SU(4) through L=8. The XX ring/open formulas agree with ED to better than 1e-10. Eigenpair residuals are below 6e-12. Four-site matrix checks verify the partial-swap costs to about 1e-16. The finite-temperature two-seam identity agrees with quadrature to below 1e-15 on two three-site XX rings. These odd rings test the matrix identity and its sign, not a vacuum Casimir coefficient. Symmetric-representation numbers are evaluations of the inherited formulas, not new many-body simulations.

[Work record] Only work_C98 was modified. No parent STATE.md, shared report, repository metadata, or external application was changed. No cluster job was needed. Local STATE.md and PROGRESS.md provide the handoff. The next substantive step would be to control the response in (3), including interpolation regions with small gaps, rather than repeat a norm estimate for one seam.
