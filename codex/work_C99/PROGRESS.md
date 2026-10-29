# C99 progress and handoff

2026-10-29. Completed the requested investigation. Read the parent project memory first, the full bound background and C99 specification, and every file in work_C69, including checks.py. All writes are confined to work_C99. The parent project STATE.md was not edited because the task forbids writes outside this folder.

[Proved] Extended C69's scaling obstruction to absolutely integrable signed combinations of calibrated rates, with explicit response-to-SIE-cost inequality D/|beta|>=ell_*/b. Derived the reduction of arbitrary single-cut unitary protocols with local controls to the unresolved capacity lower bound E(h)>=kappa A.

[Proved] Constructed an in-scope fixed d=4 Hamiltonian for which an arbitrary-state correlation length times entropy rate is unbounded. It is an XX critical sector with an onsite-gapped spectator. A long-correlated spectator state leaves the active entropy rate fixed. The Hamiltonian's A remains 2J and obeys the proposed bound.

[Proved] Constructed a one-qubit Euclidean filter with entropy rate M log N and arbitrary spectator Schmidt multiplicity N. This refutes an ancillary Euclidean SIE estimate and any direct analytic continuation of the real-time bound.

[Proved] Constructed a distant-MI example where preexisting Bell entanglement is transferred locally, giving I(A:B)=H_2(sin^2 Jt) independently of separation. A uniform nearest-neighbor chain realizes the example. It refutes a state-independent separation-suppressed MI rate, not LR.

[Proved continuum identity] Derived A=(6epsilon/pi)int dot S(t)dt/t for a common-spatial-smoothing local join with arbitrary sector speeds. SIE's integrated envelope diverges; truncation retains the regulator length and a logarithmic window cost.

[Assessed] Inhomogeneous/Floquet controls, modular generators, and Euclidean flows require respectively an actual crossing-action bound, a microscopic modular-generator estimate, and a state-restricted covariance estimate. None is supplied by the physical M and infrared CFT identification alone.

[Numerical] Ran python3 checks.py successfully. It uses only small NumPy/SciPy matrices and quadrature. CHECKS.md records the exact output and its limits. No large computations, cluster work, subagents, commits, pushes, or shell network requests. The web tool checked three primary references for SIE, local joins and the CFT modular generator.

[Open] No uniform A bound and no counterexample to it. No universal theorem excluding all entanglement-rate methods. The three explicit auxiliary counterexamples and the weighted-rate obstruction have the narrower scopes stated in PROOF.md. The main remaining route is E(h)>=kappa A or a comparably strong restricted response theorem.
