# C78 result

[assumed, unresolved target] Lemma TP remains unproved for weak XXZ off half filling. No actual mixed-chirality coefficient obstructing TP was found. This run controls the alternating bubble-ring vacuum subseries of both high chiralities, not the complete high-field integration or the actual BM infrared comparison.

[proved, mixed integration] PROOF Sections 2-4 strengthen the single-smeared high bubble to weighted L1 and L2, uniformly in the smooth sea cutoff. Its thermal replacement has arbitrary fixed power T^P, including two speed derivatives. The L2 bound controls the shared time in the order-g^2 mixed vacuum graph. All alternating bubble rings then sum analytically on a cutoff-independent coupling disk, with

    |partial_A^j(R_T-R_0)| <= C_P |g|^2 T^P,
    j<=2, P>1,

at the fixed preparation scale. These are infinite-spatial-volume statements. A uniform finite-circle construction is not supplied.

[proved, coefficient and normalization] Before a high/low split, the mixed density responses give

    Pi_+ Pi_-=p^2/[4 pi^2(Omega^2+A^2p^2)],
    f^(2)(T)-f^(2)(0)=-g^2 T^2/(48 pi A^3)+o(T^2).

Their ring determinant has sound dispersion sqrt(A^2-B^2 chi(p)^2)|p|, with g=2 pi B, agreeing with C62's already fixed trace. This nonzero T^2 term must be retained, but it is not an obstruction after matching the dressed reference.

[proved, first omitted mixed terms] PROOF Section 6 gives the complete non-ring order-g^4 vacuum coefficient. Its integrand is

    [C_(+,4) C_(-,4)+C_(+,4) P_-+P_+ C_(-,4)]/4!,

integrated against the four interaction kernels, with the common spacetime volume divided out. Here C_(omega,4) is the connected high-density four-cycle and P_omega is the sum of its three pair partitions. C73 shows that C_(omega,4) generally survives the high/low split. The expression is finite at zero temperature by C73's fixed-order theorem. Its presence invalidates a ring-only formula for the complete high action, but no failure of TP follows from it.

[proved, restricted infrared check] PROOF Section 7 compares two explicitly defined Gaussian density determinants with matching A,B and spatial versus spacetime smoothing. On scale q_h their scalar difference is O(q_h^4), with thermal replacement O(q_h^4(T/q_h)^P). This is not an identification with BM's full functional and does not control interacting histories.

[assumed, remaining work] Bound the complete mixed-cumulant and open-chain expansion exponentially in perturbative order in a fixed stronger norm, with its scalar and thermal derivatives. Then prove TP's complete lattice/reference defect, mixed thermal histories, common parameter domains, and vacuum Hessian bounds. Those tasks are unresolved, not reduced to a demonstrated bad coefficient.

[numerical] checks.py passes the final algebra and quadrature controls. An initial unsupported numerical threshold failed and was removed, with the failure and measured values recorded in PROOF Section 9. No interacting XXZ flow was computed. All writes are inside work_C78, and no independent review is claimed.
