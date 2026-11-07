# C134 progress

2026-11-07.

[Proved, inherited] C129 supplies flux minimization and the free point. Its ferromagnetic band is only (Delta+1)L^2 bounded. C125 positive spatial transfer applies to Delta>=0 on even rings, not to negative Delta.

[Assumed, work in progress] Testing a new real-axis quadrature argument for -1<Delta<=0. The negative scattering kernel gives a monotone counting function and an explicit cumulative root discrepancy. A dressed-energy cutoff may turn this into an O(1/L) energy bound. Ground-state identification and odd lengths need separate care.

[Work record] All changes stay in work_C134. Large scans will be prepared for SLURM, not run locally or submitted from this sandbox.


## Completed

[Proved] PROOF.md contains the all-length compact theorem C(delta)=3000/delta^3 on [-1+delta,0], smaller explicit constants, the odd-length ground-branch identification, and the exact dressed quadrature argument. It also derives C=2pi^2/9 at Delta=1/2 from the published supersymmetric twisted energy and C129's twist estimate.

[Proved, audit] LITERATURE.md records what was checked in Kozlowski, DDV, DKKMT, the later high-temperature theorem, and the supersymmetric-energy sources. Access failures for older full texts are explicit. Positive spatial transfer, even where established, does not by itself give a rate for D(L).

[Numerical] The 88 small and 64 moderate paired Bethe cases pass, with ED through L=12 and paired Bethe lengths through 258. Exact-identity and cutoff checks pass. CHECKS.md records the commands and numerical errors. No several-thousand-site scan was run.

[Assumed, unresolved] The full positive-anisotropy interval, Heisenberg endpoint, and a uniform ferromagnetic crossover constant remain open here. The new compact theorem has not had an independent referee review.

[Work record] RESULT.md, PROOF.md, CHECKS.md, LITERATURE.md, STATE.md and the numerical scripts are complete. run_typhon.slurm prepares the large scan but was not submitted. No outside files were modified.
