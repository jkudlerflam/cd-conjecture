# C104 local state

[Proved, scope] Completed the requested investigation and wrote RESULT.md, PROOF.md, PROGRESS.md, check.py, and CHECKS.md. The parent project STATE.md was read but was not edited because the task explicitly forbids modifying files outside work_C104.

[Proved, conditional theorem] Full simple WZW infrared factors seen by on-site semisimple symmetry obey A_WZW <= (4 pi/3) d(d^2-1) M, conditional on the WZW flux-response matching. The simply laced coefficient is 2 pi/3. Multiple levels and velocities coupled to the same microscopic factor are allowed. See PROOF.md Sections 5 and 6. No upper bound on levels is used.

[Proved, conventions] Canonical generators T have usual spin normalization on a long-root SU(2); tr_R(T^A T^B)=I_R delta_AB. The normalized generators are T/sqrt(I_R), and their stiffness is sum_a k_a v_a/(2 pi I_R). I_R <= d(d^2-1)/12. All normalized Cartans share the supplied budget 4M.

[Proved, limitations] Compact-boson integer labels and mutual locality allow every K>0. They do not prove a current-weight lower bound. No fixed-d nearest-neighbor lattice example with K tending to zero is constructed here. The general charged Gaussian bound remains open, as does a bound for independent neutral sectors. Parafermions contained in a full WZW theory with its common velocity are already included in the WZW theorem.

[Numerical] check.py passed all finite-dimensional checks. It performs no many-body simulation or lattice-to-CFT test. No cluster work was needed.

[Proved, next mathematical tasks] Establish a Gaussian aggregate estimate sum v_a <= G(d) sum v_a |Z_a|^2, or a stronger mode-weight lower bound, using a microscopic input beyond integrality. Separately control stress tensors not fixed by the measured current algebras. Do not treat the continuum radius example as a finite-d lattice counterexample.
