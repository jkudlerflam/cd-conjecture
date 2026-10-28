# C14 progress

2026-10-27. Initial reading completed. No previous PROGRESS.md existed. Read project STATE.md, STATUS.md, refined_conjecture/README.md, both Astra STATE.md files, and REPORT sections 3.1 to 3.4 and 7 (plus relevant later context). All work stays here. No network, no agents, Python with numpy/scipy only. Parent STATE.md will not be edited because the task forbids outside writes.

## Active findings

1. General energy structure factor and thermal f-sum routes need an infrared coefficient bound stronger than a norm bound on a finite-temperature total weight. The exact momentum and temperature obstruction will be written out.
2. Candidate extension: double any even quadratic Majorana chain and combine corresponding Majoranas into complex fermions. The doubled Hamiltonian conserves number, its Casimir coefficient doubles, and the canonical bond spread doubles. Need handle arbitrary coboundary gauges explicitly.
3. Priority stress test: symmetric-rank-k SU(N) rational fusion chains. Recalled candidate formulas, NOT YET ESTABLISHED HERE: bond eigenvalues 2(H_{k-j}-H_k), j=0,...,k, norm H_k; c=k(N^2-1)/(N+k); v=2*pi/N in this normalization. These would give A/(2*pi*M)=k(N^2-1)/(N(N+k)H_k), exceeding one already at k=2,N=7 and unbounded for N much larger than k. Must derive the fused normalization and thermodynamic Bethe velocity, and calibrate at SU(2) spin 1 and SU(N) fundamental. Do not quote a counterexample until this audit closes.

## Next

Read existing fusion-chain notes for conventions. Derive the nested k-string density and dressed hole energy. Build independent small-chain checks where affordable, with at least two estimators for any numerical CFT claims. Record proofs and conditional assumptions distinctly in RESULT.md.

## Fusion audit, step 1

Derived the zero-temperature nested k-string density equation. The driving kernel B_kk has Fourier transform sum_{j=1}^k exp(-(2j-1)|omega|/2); the diagonal same-rank kernel is 2 cosh(omega/2) B_kk and the neighboring-rank kernel is B_kk. It cancels, leaving rhohat_r=sinh((N-r)omega/2)/sinh(N omega/2), independent of k. Dressed hole energy is 2*pi*rho_r and momentum derivative 2*pi*rho_r, giving v=2*pi/N. This supports the candidate violation. Ground energy for the ferromagnetic-reference-zero bond is e0=-(2/N) sum_{j=1}^k [psi(1+(j-1)/N)-psi(j/N)]. It reduces to ULS minus 1 at k=1 and TB/2 energy -2 at N=k=2.

Next decisive check: generalize the existing rank/string TBA by placing its source at string length k, and extract c separately from free energy and entropy. The high-temperature entropy must match log binomial(N+k-1,k), a normalization/completeness check.

## Fusion audit, step 2: counterexample supported

TBA generalized to source string k in fusion_tba.py. SU(2),k=2 yields c_free=1.50400 and c_entropy=1.50497 at T=.02, approaching 3/2; high-T entropy at T=100 is 1.098567, approaching log 3. SU(8),k=2 yields c_free=12.6517 and c_entropy=12.6689 at T=.005 (M=24,K=512,R=32), approaching 12.6. Both imply A/(2*pi*M) above 1, approaching the analytic 1.05. High-T entropy is 3.583458, approaching log 36. Residuals below 1e-9. These are two related but distinct thermodynamic estimators, not independent lattice ED. Need cutoff tests and small-chain bond/ED calibration.

Explicit k=2 bond: on Sym^2 C^N, E_ab are Schwinger-boson generators at occupancy 2, X=sum_ab E_ab tensor E_ba, and h=-X^2/8+3X/4-I. X eigenvalues are 4,0,-2, hence h eigenvalues 0,-1,-3. Gauge-minimized half-spread remains 3/2 by SU(N) twirling and irreducibility. At N=2 this is (S.S-(S.S)^2)/2. General eigenvalues are 2(H_{k-j}-H_k), and M=H_k. Standard SU(N)_k infrared identification gives ratio k(N^2-1)/(N(N+k)H_k). N=k tending to infinity makes it grow as k/(2 log k), so there can be no finite universal dimension-independent constant multiplying the bond norm, if the standard fusion-chain IR input holds.

## Fusion audit, step 3: algebra and two calibration routes pass

Independent fusion_checks.py builds Schwinger-boson generators and the two-site bond. At N=7,8,k=2 it confirms X levels 4,0,-2 and h levels 0,-1,-3 to <8e-15. A direct braided Yang-Baxter check of R_j(u)=product_{l=k-j+1}^k (l-u)/(l+u) passes at (N,k)=(2,1),(2,2),(3,2), maximum error 7.8e-16. Its derivative matches h. The A-type restricted-Y-system Rogers-dilogarithm sum reproduces c=k(N^2-1)/(N+k) to 2e-13 through N=k=16.

TBA cutoff increase N=8,k=2 to M=48,K=1024,R=48 yields c_free=12.614875,c_entropy=12.619620 at T=.002; ratios 1.051240,1.051635. N=k=5 (d=126) gives c_free=12.026884,c_entropy=12.030203 at T=.001 with M=40,K=1024,R=40, approaching c=12 and ratio144/137=1.051095.

Independent lattice ED normalization calibration implemented in ed_calibration.py. SU(2) fundamental rings L=8,10,12,14,16; spin-1 TB/2 rings L=6,8,10,12; SU(3) fundamental L=6,9. Both Casimir and spatial-entanglement estimators recorded. Small-size entropy estimates overshoot, as expected from existing notes; do not call them precise thermodynamic c. All sparse residuals <1.2e-10. Largest sector dimension 73789.

Main conclusion now: the proposed universal bound is false at the requested physics standard, using the standard thermodynamic Bethe/string and WZW/CFT identifications. The stronger assertion of any dimension-independent constant is also false in this family. Remaining work: write full derivation and precise assumptions, finish free-Majorana extension, general-route obstructions, LL and symmetry/coset scope, and reproduction instructions.

## Draft deliverable complete

RESULT.md now contains the explicit rank-two polynomial counterexample, full fusion normalization and Bethe derivation, the Majorana doubling theorem with arbitrary-gauge handling, general sum-rule/variational obstructions, conditional LL and current-algebra bounds, graded-swap and stack checks, and consequences for the entropy window. Scope exclusions are explicit.

Further cutoff check N=8,k=2,M=64,K=1024,R=48,T=.002 gives c_free=12.613358,c_entropy=12.618197, ratios1.051113,1.051516; residuals <1e-9. ED Casimir fits with free bulk energy yield c=1.005379 (SU2 fundamental) and1.517606 (TB/2); the spatial-entropy drifts are larger and are reported rather than treated as precision checks.

Running one stronger-excess TBA point family N=k=8 at T=.003,.001, cutoff48, grid1024, half-length48. It should approach c31.5 and ratio1.448752. Then final consistency/readability audit and file inventory.

## Completed, final resumable state

The last TBA run finished. N=k=8,T=.001 gives c_free31.634722,c_entropy31.650391, approaching31.5; the ratios are1.454948 and1.455668, approaching1.448752. Both temperatures and the exact reproduction command are now in RESULT.md. All TBA solver/tangent residuals across all saved runs are below1e-9. The m=1,2,3 Majorana nuclear-norm and arbitrary-gauge extraction tests pass, and direct doubling has zero displayed matrix and number-commutator residual. Every Python file compiles. RESULT.md and this file contain no em dashes. No calculation remains running.

Deliverables: RESULT.md (verdict, assumptions, derivations, class boundaries, numerical tables, reproduction commands and ranked next steps); fusion_tba.py with seven TBA JSON files; fusion_checks.py/.json; ed_calibration.py/.json; majorana_checks.py/.json. All are inside this directory. No outside writes, network, cluster jobs or delegated agents were used.

Final verdict: A<=2*pi*||h|| is false at the standard Bethe/CFT physics level. Symmetric SU(N) rank-k fusion chains give ratio k(N^2-1)/(N(N+k)H_k); N=7,k=2 already violates, and k=N disproves every dimension-independent constant. Positive extension: every even quadratic Majorana/BdG NN chain satisfies A<=pi*||h||, with the full allowed coboundary minimization. Finite-density graded-swap and symplectic IR surveys remain undone and are marked as such. No implication for unbounded c at fixed d is claimed.

Highest-value continuation: independently check primary-source fusion-chain normalization and IR identification, or solve finite-size Bethe equations for the SU7/SU8 rank-two model. Do not repeat the general proof search for the false universal statement. Parent STATE.md/STATUS.md were intentionally not edited because outside writes are forbidden.
