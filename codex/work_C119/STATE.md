# C119 handoff

2026-10-29. All work is confined to this folder. Project-level STATE.md was read but not modified under the task's write restriction.

[Proved] New uniform bound for even rings in XX, independent spin XX stacks, C106, and its a=0 massive comparison:

    R_L(beta)<=3log d+(10pi/3)M beta/L.

For one shared parity sum the intercept is only 3log 2. The proof controls the negative part of the second k derivative of tr log[2cosh(beta H(k)/2)], then its Fourier quadrature error. It uses a nearest-neighbor Bloch matrix and an even band multiset for the spin parity step. It does not prove the claim for arbitrary quadratic spin chains or interacting chains.

[Proved] Exact finite products include both fermion grids and parity projections. Spin XX stacks have independent projections; C106 has a shared one and a different torus partition function. First L^-2 lattice corrections are explicit for XX and aligned C106 rings. C106 coefficients diverge as its avoided-crossing scale closes, while the full ratio stays uniformly bounded.

[Proved] Slow XX sectors vanish exponentially at fixed beta/L as v->0, but contribute order one at beta proportional to L/v. For changing CFT spectra, shrinking dimensions can prevent suppression. A compact boson has r~(1/2)log K at fixed y and c=1. The XXZ microscopic attempt has Kv->pi/4 and a finite fixed-tau ratio, so does not violate the target.

[Proved] Fixed-gap free modes disappear from the beta=tau L limit. Near-critical massive modes contribute at crossover scales. Large R at one scale is not alone evidence for many fast infrared sectors.

[Open] No fixed-d,M interacting counterexample or universal interacting proof. The sharp free ordered slope is k=pi/2,f=1 from the inherited A<=2M. The strongest stated interacting data envelope corresponds to k=pi^2/4,f=H_(d-1), conjectural and only an ordered-slope statement. A Euclidean seam response/effective transfer-rank lemma with the required constants is stated in PROOF.md sections 7-8.

[Numerical] All small checks pass, largest matrix dimension 256, complete C106 trace error 4.98e-14, 360 scan points, run under one second. CHECKS.md gives the measured limits and tolerances. checks_small.json is reproducible output, not orientation material.

[Work record] RESULT.md is the summary, PROOF.md contains derivations, CHECKS.md records verification, and checks.py reproduces it. run_typhon.slurm prepares an extended scan. The sandbox cannot reach Typhon and no job was submitted. No independent review was performed.
