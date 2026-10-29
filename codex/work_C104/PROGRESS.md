# C104 progress

[Assumed] The setting and infrared interpretation are those of ../../specs/_bound_background.md, including both updates. Read the pairing README and REPORT sections 9.68 through 9.75. The project STATE.md was read for orientation.

[Proved] The supplied stiffness budget bounds velocity times current weight. Integer vertex labels alone do not bound the Gaussian metric. A compact-boson calculation will make the missing implication explicit.

[Proved, conditional on the WZW response dictionary] An upper bound on WZW levels is unnecessary for the requested relaxed bound. The ratio c/k decreases with level. A representation-dependent bound follows from normalized Cartan stiffnesses. A uniform polynomial bound in d is being checked using restriction to a long-root SU(2).

[Assumed] The WZW stress tensor must exhaust the sectors being bounded and the lattice symmetry must couple to their affine currents. Neutral cosets require separate treatment.

[Proved, scope of work] No files outside work_C104 have been modified. No numerical simulations or network shell commands have been run. Literature is being checked with the browsing tool. Project memory will be recorded locally because the task forbids modifying the parent STATE.md.

[Proved, conditional on the stated WZW dictionary] The representation estimate now gives A_WZW <= (4 pi/3) d(d^2-1) M for any finite collection of fully symmetry-visible simple WZW factors, including several levels and velocities coupled to the same microscopic simple factor. This uses the Cartan torus budget and no upper bound on k. The sharper constant retains the on-site Dynkin indices.

[Proved] The pairing README mixes fundamental-trace normalization with trace in the actual on-site representation. For generators normalized by tr_R(t^a t^b)=delta_ab, the correct WZW response is D_t = k v/(2 pi I_R), with I_R defined in the long-root SU(2) normalization. For spin S, I_R=S(S+1)(2S+1)/3. The missing factor matters beyond fundamental representations.

[Proved, limitation] No fixed-d nearest-neighbor lattice family with K tending to zero has been established in this run. Arbitrary compact-boson radii disprove the proposed continuum inference, but do not settle the stronger finite-d lattice question. The general U(1) conclusion remains open.

[Proved, completed write-up] RESULT.md and PROOF.md now contain the polynomial WZW theorem, the Gaussian quantization obstruction with its precise scope, the neutral-sector distinction, and normalization corrections. STATE.md records local memory. The general charged Gaussian goal is still open, and no universal A/M counterexample is claimed.

[Numerical] check.py passed small deterministic checks: spin matrices through d=8; all 914 SU(2) decompositions through d=16; 401 root-system ratio checks; Heisenberg and spin-1 stiffness conventions. Maximum spin trace residual was 4.441e-16. No heavy numerics were performed.
