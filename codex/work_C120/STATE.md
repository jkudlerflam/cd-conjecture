# C120 handoff

2026-10-29.

[Open] No universal K(d), k and f(d) for the target ratio bound were obtained. No counterexample to the target or to the velocity-weighted bound was found.

[Proved] Arbitrary-reference Dmax, Petz and sandwiched divergences below one, and negative log fidelity have no dimension-only upper bounds. Faithful qubit counterexamples are explicit. Support D0 vanishes for finite-temperature Gibbs states; the zero-temperature support limit is different.

[Proved] Conditional min-entropy gives the exact identity R = Hmin(sigma)-Hmin(rho)+log[mu(Y)/mu(X)], with mu(W) = inf Tr Omega subject to W <= I_C tensor Omega. The entropy difference is bounded by 2 log dim C. The remaining ratio retains beta times the doubling-energy defect. Petz/sandwiched/fidelity identities likewise retain unnormalized mixed traces. Forward lower-order substitutions have the wrong direction. The reverse fixed-order substitution loses a factor L already on XX.

[Proved] Strictly positive completely positive qubit maps can violate any dimension-only ceiling on (Tr Phi^m)^2/Tr Phi^(2m). Explicit depolarized rotations have Birkhoff coefficient a and log ratio log[(1+a^m)/(1-a^m)]. A real transfer-matrix bound with a negative-pair penalty is proved. Zero penalty still leaves the growing history dimension.

[Numerical] Full ED through XX/Heisenberg length 8 and SU(3) length 6. Twenty-seven temperature/model pairs, six Petz orders, conditioned sandwich and max calculations, fidelity/trace distance, one-site conditional-min witnesses. Twenty-four finite Trotter networks checked by temporal contraction at both lengths. Full XX thermal traces through paired length 16384, ground-state Slater overlaps and four-site conditional min-entropy through paired length 1024. Last full run took about ten seconds. Exact Petz identity residual below 2e-13.

[Scope] Finite-temperature large-ring fidelity and full conditional-min-entropy optimizations were not run. No large local calculation, cluster submission or independent referee review. The project-level STATE.md remains untouched under the task's folder restriction.

[Work record] RESULT.md contains the summary, PROOF.md the derivations and primary references, CHECKS.md the candidate outcomes and limitations, PROGRESS.md the work record. checks.py and ring_tools.py reproduce checks.json and summary.json. A new proof needs independent control of a retained unnormalized quantity or of the effective spatial transfer dimension at separation L.
