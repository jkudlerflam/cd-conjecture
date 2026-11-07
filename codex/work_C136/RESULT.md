# C136 result

2026-11-07. Task 1 is proved, with an explicit constant uniform in filling. The bound uniform through Delta=-1, and therefore the full requested theorem, remain unresolved. The new proof has not had an independent referee pass.

## All fillings on a compact attractive interval

[Proved] For every 0<delta<=1, every bond of the stated U(1) class with -1+delta<=J_z/J<=0, every field and DM angle, and every integer L>=2,

    L(E_(2L)-2E_L) <= [6 * 10^15 / delta^(25/3)] M.

Here J=sqrt(J_perp^2+D_DM^2), S=sigma/2, and M is the optimized half spectral spread of the physical bond including its field. The convention at L=2 is the sum of two oriented bonds. No parity or vacuum-compatibility restriction occurs.

[Proved] The underlying sector theorem is stronger. At every 0<=n<=L/2,

    |F_L(n)-J L e_bulk(Delta,n/L)|
                       <= [10^15 / delta^(25/3)] J/L.

[PROOF.md](PROOF.md), Sections 6 and 9, gives smaller explicit expressions K_all(delta) and C(delta)=5K_all(delta)+pi^2/4. Section 5 gives substantially smaller constants K_(delta,eta) when n/L<=1/2-eta. The large constants above are analytic majorants, not numerical fits. No claim of literature priority is made.

[Proved, mechanism] Negativity of the scattering kernel gives a globally increasing counting function at every filling. Projecting its cumulative function onto the occupied counting interval gives an O(1/L) cumulative-root estimate uniform in density. At the finite Fermi boundary, an exact dressed quadrature identity includes a chemical-potential subtraction. The bulk curvature is bounded by pi sin gamma, so the density mismatch contributes only a quadratic error. An exponentially decaying resolvent bound compares the finite-interval density to the full-line sech density and controls the boundary correction to the dressed velocity. This supplies a variation bound uniform in the Fermi boundary. C134's clipping argument then works at all fillings, including the half-filled limit and odd matched-sector doubling.

[Proved] This avoids a defect in the initially proposed dilute patch. C129's estimate n^3/L does not cover all sequences with n/L->0. The present theorem covers those sequences directly. PROOF Section 7 lists each half-filled ingredient of C134 and its replacement. [LITERATURE.md](LITERATURE.md) records what is and is not supplied by Kozlowski's finite-boundary expansion.

## What remains at the ferromagnetic crossover

[Proved, obstruction to the existing estimates] Write epsilon=1+Delta. C134's actual cutoff starts at L of order epsilon^(-1), whereas the C129 band reaches only L of order epsilon^(-1/2). Beyond that cutoff, C134's large-size coefficient still diverges as epsilon^(-3/2). Multiplication by L in the desired defect cancels its explicit 1/L suppression. Thus tracking its delta and L dependence does not close the crossover. The new all-filling constant also diverges.

[Assumed, exact missing input] It would suffice to prove explicit A,K_* independent of Delta and density such that

    |F_L(n)-J L e_bulk(Delta,n/L)| <= K_* J/L
                    whenever (1+Delta)L^2>=A.

Then the whole requested class would have C=max{A,5K_*}+pi^2/4, using C129 below that scale. A direct one-sided defect estimate beyond the same scale would also suffice. Neither has been proved here. At Delta=-1 the physical defect is zero, but that endpoint fact does not give uniformity in its neighborhood. No counterexample to the desired uniform bound was found.

## Numerical scope

[Numerical] All 1,032 paired Bethe cases passed on Delta in {-0.9999,-0.999,-0.99,-0.9,-0.5,0}. Every sector was checked at L=2,...,12 and L=16,17,32,33,64,65. Six sectors per length were checked at L=128,129,256,257, with paired rings through 514. ED covers every sector at L<=12. The largest energy discrepancy was 1.25e-14. The largest sampled sector ratio L Lambda/M_0 was 2.3431457505076203, at the two-site free point. All applicable variational and compact theorem ceilings passed.

[Numerical, incomplete] Of 54 independent finite-boundary continuum diagnostics, 43 passed after refinement. Eleven remain underresolved near Delta=-1. Their errors, convergence checks, and the quantities tested are recorded in [CHECKS.md](CHECKS.md). These are floating-point checks, not certified enclosures, and do not certify the new proof. The conservative interacting clipping threshold was not reached in the local Bethe grid. Four additional free-fermion checks reach that branch using constant-cost closed formulas.

[Assumed, pending computation] [run_typhon.slurm](run_typhon.slurm) prepares all 98,304 requested all-sector cases with L=2,...,256. [run_intermediate_typhon.slurm](run_intermediate_typhon.slurm) prepares higher-resolution continuum diagnostics. Neither was submitted. Typhon is unavailable in this sandbox. The full all-density grid through a few hundred sites and a universal final constant check remain unfinished; no universal final constant is claimed.

[Work record] PROGRESS.md, PROOF.md, RESULT.md, CHECKS.md, LITERATURE.md, STATE.md, scripts, and numerical outputs are confined to work_C136. Project-level STATE.md was not modified because this task prohibits outside writes. The optional positive-anisotropy task was not attempted because the crossover and full assembly remain open.
