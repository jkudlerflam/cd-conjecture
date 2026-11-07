# C132 progress

[Proved] Read the project memory, target, full cd_free_route sources, updated bound background, required REPORT sections and reconnection background. All writes were confined to this folder.

[Proved] Audited the marginal equivalence, periodization, gauge cancellation and arbitrary CPTP extension. Identified equality of the dual-M and trace norms on the marginal-difference space, and the missing log-L step in the source trace-norm obstruction.

[Numerical] Copied and ran the unmodified quick checks. All pass. Original --long is prepared for Typhon rather than reported as reproduced.

[Proved] Derived the support-reduced recovery SDP and its dual. Read the explicit rotated recovery formula in arXiv:1509.07127v3 through the browsing tool and checked its admissibility directly.

[Numerical] Solved small SDPs for every requested model. Refined qubit L = 6 values with Clarabel. Retained residuals, dual diagnostics and inaccurate statuses. An early larger SCS run was interrupted for slow convergence, and larger tasks were written as SLURM scripts. A process-inventory command was unavailable in the sandbox; the tool session interrupt stopped that run instead.

[Proved] Solved the ferromagnetic Heisenberg mixture SDP analytically by angular momentum addition. Petz is optimal and its positive defect is an ansatz artifact.

[Numerical] Scanned near-ferromagnetic XXZ, near-critical TFIM and 16 generic gauge-quotient bonds. Corrected incomplete sparse FM mixtures using exact Dicke vectors. Followed the worst random cases through L = 12 and refined their full small-window SDPs. Nonmonotonicity found, no divergence established.

[Proved] Extracted the explicit quadratic-class finite-size constant K = pi^2/6 from momentum smoothing and a parity projection of the optimized gauge.

[Proved analysis] Recorded the exact stationarity inequalities and why they do not give the needed upper bound, the singular/noncommuting modular issue, and the prescribed-marginal obstruction to arbitrary longer-ring preparation.

[Numerical] verify_c132.py passes the FM formula, arbitrary non-Petz stabilization and an independently contracted complex Choi objective.

[Assumed, open] General P prime, even at d = 2, remains unproved. The asymptotic optimal factor and uniform large-L boundedness cannot be decided from completed finite-size SDPs. RESULT.md distinguishes this from the proved class results and the unsubmitted cluster continuations.
