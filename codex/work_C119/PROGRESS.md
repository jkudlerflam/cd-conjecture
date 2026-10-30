# C119 progress

2026-10-29.

[Work record] Read the project memory, the complete bound background with updates, REPORT sections 9.68 through 9.79, and the C116 derivations, code and checks. Read C106's model and normalization. All writes and Python execution are confined to work_C119. No cluster connection or network command will be attempted.

[Proved] The relevant parameter for a sector is y = beta v/L. Suppression as v tends to zero holds at fixed beta/L and fixed conformal spectrum, not at beta proportional to L/v. Shrinking conformal gaps require a separate check.

[In progress] Exact parity-resolved fermion traces, subleading scaling terms, a uniform quadrature bound for the requested free models, and attempted counterexamples using C106 crossover and compactification-radius growth.

[Open] No general interacting estimate or fixed-d counterexample to the velocity-weighted bound has been obtained.

[Proved] Derived a uniform free-family estimate R<=3log d+(10pi/3)M beta/L on even rings using trace convexity and Fourier quadrature. This bypasses divergent band curvature at the C106 avoided crossing.

[Proved] Derived the full parity products and the first lattice corrections for XX and aligned C106 rings. C106's shared parity projection changes its order-one torus ratio relative to independent spin XX chains.

[Proved within the specified CFT description] Compactification growth can make the ratio large at fixed central charge, so a large ratio does not alone count fast sectors. The XXZ attempt near the ferromagnetic endpoint yields a finite fixed-tau theta function because Kv stays finite.

[Numerical] Small checks completed in less than one second with matrices no larger than 256. Exact C106 trace error is 4.98e-14. Subleading formulas, crossover checks, slow and massive limits, and a 360-point uniform-bound scan pass. CHECKS.md records the numerical scope and roundoff limits.

[Work record] RESULT.md, PROOF.md, CHECKS.md and local STATE.md written. The extended scan is prepared for SLURM but not submitted. No network command or Typhon connection was used. General interacting existence/impossibility and the universal seam estimate remain open.
