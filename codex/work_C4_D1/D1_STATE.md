D1 closed-book working memory, 2026-10-21

Scope: work only in work_C4_D1, read no files predating this run, no browsing, no other agents. The project STATE.md was intentionally not accessed. This file records this attack only.

Result: no certified in-scope counterexample to c_T <= d-1 or to its effective-charge variant. This is not a nonexistence proof. D1_report.md contains the full derivations, assumptions, compact-family inventory, unresolved cases, and reproduction commands.

Main proved findings:
1. Reconstructed spin-1 twisted-projector Hamiltonian h=cot(2theta)P1-tan(3theta)P0 has a bulk anti-Hermitian entry i cos(2theta)/cos(3theta) between |+-> and |00>. It cannot telescope away. Its identification with the standard IK regime III and the k dictionary remain explicitly assumed.
2. TL staggered candidate H=sum[2 cos(gamma)e_j-{e_j,e_(j+1)}] has a surviving imaginary three-site interaction. Pair blocking gives d=4 NN but preserves non-Hermiticity. A small alternating six-vertex transfer-matrix logarithmic derivative matches this H up to scale and constant.
3. Complex eigenvalues at theta=.19 (IK L=3,4) and gamma=.4 (TL L=6) rule out any positive-metric quasi-Hermiticity for those finite-size Hamiltonians.
4. Bare positive-metric graded permutations with odd colors have polarized ground states and explicit quadratic impurity branches. Fixed density or chemical potentials are different problems.
5. The displayed mixed osp rational density is non-Hermitian. Pure orthogonal and symplectic rational densities are Hermitian. Their continuum classification is not fully resolved.

Numerics: d1.py is standalone and contains all routines. Commands: algebra, validate, symbolic, transfer, ed, momentum, audit, analyze. ED reached spin-1 L=12 (73789 states at zero magnetization), spin-1/2 L=20 (184756 states). All nonnegative charge sectors were audited at largest sizes, with spin-reversal conjugation relating negative charges. Momentum-projected residuals <=2.3e-11. Logarithmically shrinking neutral gaps are consistent with, but do not prove, a non-compact limit. Unverified velocity formulas are explicitly flagged in d1_diagnostics.json and the report. No complete conformal-tower identification was obtained.

Outputs: d1_algebra.json, d1_validation.json, d1_complex_spectrum.json, d1_spectra.json, d1_momenta.json, d1_ground_audit.json, d1_diagnostics.json. The d1_* companion Python scripts are earlier runnable versions also incorporated into d1.py. Progress text files contain run results. Do not read unrelated project files on continuation without relaxing the closed-book instruction.

Outstanding: literature/model dictionary verification; exact continuum velocity and k derivations; full classification of D_n^(2), A_(2n)^(2) compact branches, Sp and G2 continua, and finite-density graded chains; any alternative Hermitian non-compact realization; a successful range-shortening or extra-sector construction. Do not report these as settled.
