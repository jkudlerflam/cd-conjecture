# C122 handoff

2026-10-29.

[Open] The general scalar mu-ratio target and the bound on A remain open. No counterexample to either was found. All work stays in work_C122, so the parent project memory was not edited.

[Proved] Every full-order comparison Y<=c X has log c>=beta Tr(sigma V). On even XX rings, log c>=(2beta/L)csc(pi/L)-2log 2. At beta=tau L this grows as 2tau L/pi. Thus sharper full expansional norms cannot give the requested O(tau) comparison through operator order.

[Proved] Folding into d^2-dimensional rungs and counting support-growth events gives G(b)=exp[6b+2(exp(2b)-1)] for the needed half-time products, b=beta M. The resulting scalar bound is |log(muY/muX)|<=min(4beta M+2w log d,12beta M+4(exp(2beta M)-1)). No improved general R_L or A estimate follows. A factorial tail gives a sufficient locality width exponential in beta M. Local congruence approximation in mu costs q times the operator error.

[Proved] Exact low-energy decomposition is Q_C=beta Delta_L+log(muP0/muP1)+eta0-eta1. A gap delta~v/L alone does not control eta uniformly, since the bound contains d^(2L)exp(-beta delta). Even exact ground projection retains Delta_L. Remaining lemma is dual-witness transport with the partial-trace constraint and cost K'(d)+k M f(d) beta/L, PROOF Section 7. It is explicitly unproved.

[Numerical] 70 full finite-temperature and 14 ground-projector SDPs completed. Main cases XX/Heisenberg L=4 and SU(3) L=3, four physical seam sites; additional six-site qubit seams and L=2 calibrations. Ratio brackets have width <=3.11e-6. All sectors retained. Ground-state offsets are 0.11415, 0.01726, 0.04879 in the main cases. Widening changes the offset, not the beta slope. No interacting asymptotic scaling inference from these sizes.

[Work record] RESULT.md summarizes the outcome. PROOF.md has derivations and primary literature links. CHECKS.md records numerical accuracy and limitations. checks.py, sdp.py, ring_tools.py and analytic_checks.py reproduce the results. No large local work, external write, cluster run or independent review.
