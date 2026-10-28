# C25 local state

2026-10-28. Completed G3's specified C18 tagged-vacancy calculation.

[Proved] Four occupied sites, wrapped lengths 5 and 6. Exact code after
Laplacian replacement. All second-order source measures and their helix
excitations are retained. The limiting core is sum_j <D_j> h_j, so
m_p=1/r_p leaves an O((w/U)^2) mismatch. The additional finite-source
error is O(s/Delta_wrap) with compensator-dependent constants.
Delta_wrap includes the charge tunneling gap, which scales as w^8/U^7
for the two-vacancy cell. Mixed u lambda vanishes by Bell-syndrome
conservation. G1's displayed syndrome-mixing positive check is FF under
its strict local positivity condition; a compiled hopping core is untested.

[Numerical] 27 parameter points, 108 bonds. Maximum source sector 3150,
maximum pair product 9,922,500 per Pauli channel. Independent checks pass.
See RESULT.md Sections 10 and 11 for sizes, estimators and limitations.

Files: RESULT.md has the complete derivation; PROGRESS.md has the handoff;
composition.py and checks.py reproduce composition_results.json and
checks_results.json. Only this folder was written. No network or agents.

[Assumed next task] Apply the calculation to an actual hopping-producing
G1 cell once available, including its scalar length interaction and
potentially nonzero mixed term. No route-A exclusion or CFT result follows.
