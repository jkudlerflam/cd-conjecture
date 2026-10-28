# C35 progress and handoff

2026-10-28. Finished the bounded investigation. All writes are confined to work_C35. No network, delegated agents, or project-code imports. RESULT.md is the final report. The project STATE.md was not edited because the task forbids writes outside this folder.

## Results

[Proved conditional on C14 Bethe/CFT inputs] A<=pi M log_2 d holds for all symmetric SU(N) rank-k chains. Equality occurs only at N=2,k=1. Every k>=2 member has Q=A/(pi M log_2 d)<3/4. The fundamental SU(3) chain is the runner-up, Q=0.8412396714. Proof is in RESULT Sections 3.1 to 3.3.

[Proved asymptotics from the same formulas] The complete symmetric-family envelope is Theta(log D/log log D) for dimension cutoff D. For fixed k/N, the asymptotic coefficient of [A/(2pi M)] log log d/log d is maximized at k/N=1 and equals 1/(4 log 2). The globally optimal leading coefficient, beyond fixed-ratio rays, was not proved.

[Proved conditional comparisons] SO(N) vector rational, zero-field critical XXZ, finite-density gl(1|1), and the C14 quadratic Majorana class obey the logarithmic bound. A general gauge-compatible stacking proof is included.

[Assumed additional finite fusion identities; proved conditional implication] The rectangular SU(N) extension follows if its rank-one fusion restriction and full SU(N)_k velocity identification hold in the specified normalization. The dimension comparison and norm implication are proved. The full rectangular R-matrix construction was not independently proved, so this extension has a stronger assumption label than the symmetric result. Antisymmetric spectra have the same explicit qualification.

[Located obstruction] Sp(2n) for n>=2 lacks a complete sector/velocity identification in the supplied material. Finite-density higher-rank graded chains lack the needed dressed velocities including chemical potentials in the norm. Zero-field mixed graded swaps have quadratic branches and do not test the CFT conjecture.

[Proved conditional consequence] If the proposed bound holds, c^2<=3 R(log d)^2/[(1-epsilon)log 2]. There is no c bound at fixed d without a lower velocity or thermal-window bound.

[Located obstruction] Entropy and energy capacities plus a velocity cap do not bound A. An explicit positive-heat-capacity test equation of state demonstrates this lack of implication. A thermal-current proof needs a bound of order M log d on the low-temperature IR current susceptibility divided by T^3. Unprojected current norms and sum rules have the wrong temperature powers and include UV weight.

[Heuristic, with proved finite ingredients] Route-A comparisons use A_phys=b A_log, the pair-resolvent bound M(K2)<=lambda^2 ||V||^2/(2Delta), and the independent mobile length velocity 2 J_len (b+nu) sin(pi nu). The recorded selector and weak-coupling hierarchies satisfy the resulting sufficient inequalities with large margins. This does not complete a missing microscopic all-n compiler or prove its infrared limit. C28's later restricted field-rank obstruction is included in the status audit.

## Computations

[Numerical, formula and local checks] checks.py scans N=2..2000,k=1..2000, totaling 3,998,000 pairs, checks rectangular dimensions for N=2..60,k=1..30 and all heights, and diagonalizes SO local bonds d=3..8 and Sp local bonds d=2,4,6,8. It also checks unrestricted Hermitian coboundary cancellation by swap averaging. Results are in checks.json.

[Numerical, ED and failed Casimir fit] sp_ed.py uses Sp(4) L=4,6,8,10; Sp(6) L=4,6,8; Sp(8) L=4,6. Maximum sector dimension 63,504, maximum eigenvector residual 8.04e-11. Full-Hilbert-space controls at (d,L)=(4,4),(4,6),(6,4),(6,6),(8,4) match sector ground energies. The runs take about 1.3 seconds total. Results are in sp_ed.json.

[Numerical failure, not a violation] Mixed-length fits give inconsistent A, including negative values and apparent Q>1. Every such fit is explicitly rejected in the JSON and RESULT. Sp(4) L=4,8 alone gives Q=0.780369, but two small sizes do not estimate the asymptotic value reliably. A vacuum-sector Bethe analysis or larger compatible rings is required.

## Reproduce

    OPENBLAS_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 sp_ed.py

## Remaining research

(i) Complete Sp(2n) and higher-rank finite-density graded Bethe comparisons. (ii) Independently prove the rectangular fusion restriction in the chosen normalization. (iii) Prove or disprove the stress/current susceptibility estimate in RESULT Section 7. (iv) Recheck any future Route-A implementation against its physical velocity sum and full optimized microscopic norm. No universal logarithmic theorem or counterexample is claimed.
