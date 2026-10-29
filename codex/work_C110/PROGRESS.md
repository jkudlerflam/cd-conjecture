# C110 progress

2026-10-29.

[Work record] Read the project memory, the complete updated specs/_bound_background.md, REPORT.md Sections 9.68 through 9.77, and all six files in work_C107. All modifications are confined to work_C110. A local STATE.md will carry the memory handoff because the task forbids changes outside this folder.

[Proved, inherited] C107 establishes compactness of the gauge quotient, continuity of each Q_L, the sign-changing critical XY finite-size correction, and the distinction between an exact upper semicontinuous extension and a finite upper semicontinuous majorant. These results will be cited rather than re-derived.

[In progress] Testing the proposed infimum over twists and lengths, gauge covariance of boundary prescriptions, and gauge covariance of weighted smooth-cutoff Hamiltonians. The target remains a finite physical extension or an actual fixed-d divergent-A sequence.

[Open] No general semicontinuity proof or divergent-A bond family has been established.

## Completion

[Proved] Derived the XX flux formula and showed that twist infima select the wrong coefficient. Used the inherited critical XY finite-size expansion to rule out raw finite-tail infima even after a vacuum residue class is selected.

[Proved] Derived the seam coboundary defect and an explicit gauge orbit of a gapped product chain whose unrestricted twist supremum diverges. This is not physical A and does not violate the conjecture.

[Proved] Derived both an exactly solvable transverse-coboundary smooth response and a canonical-gauge interacting example with a unique product ground state and gap >=2-4|s|. Proved a strictly positive weighted deficit at order 1/L in the latter. Also identified its frustration-free gauge, which removes the response.

[Proved] Constructed a finite gauge-optimized open-cutoff deficit, proved upper semicontinuity and the bound 0<=D_w<=2WM, and stated the missing comparison with physical A. This is a concrete conditional route, not a completed identification.

[Numerical] checks.py passed. Spin matrices through dimension 256 verify the XX flux energies, the gapped-chain spectral bound, and the weighted variational bounds. Scalar sums verify the asymptotic trial coefficient.

[Open] The main problem remains unresolved. No actual divergent-A family and no general finite upper semicontinuous physical extension were found. All requested reports and the local memory handoff are written inside work_C110.
