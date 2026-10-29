# C111 session memory

2026-10-29

[Proved] The requested smooth-modulation lemma is false as stated, including for bonds that attain the gauge-optimized M. PROOF.md contains a fixed-dimension counterexample; RESULT.md is the summary. This does not disprove the velocity-weighted bound A <= C M f(d).

[Proved] Four-state spectator V = A direct-sum B, each a qubit. Parameters a = 3/4, lambda = 1/10, 0 < delta <= 1/20. AA bond is -a + delta(XI+IX)/2 + lambda(ZI-IZ), BB is swap, AB and BA are +I. Bond spectrum lies in [-1,1]. B singlet and B product witnesses with equal one-site marginals force optimized M = 1. The minimizing representative is unique up to scalar shift for this spectator. Uniform ground state is all A with X = -1, energy density -a-delta and gap 2 delta. Other label sectors are excluded by open-segment bounds, using P12+P23 >= -I on three B qubits.

[Proved] Stack an ordinary XX qubit at each site. Now d = 8, optimized M = 3/2 attained by the displayed bond, and uniform infrared A = 1, using the supplied XX identification and the proved spectator gap. No uniqueness of the optimal gauge is asserted for the stacked bond.

[Proved] In the weighted all-A sector, define u_i = (w_i+w_(i-1))/2 and v_i = w_i-w_(i-1). Its deficit D is sum[sqrt(delta^2 u_i^2 + lambda^2 v_i^2)-delta u_i]. The unrestricted spectator deficit is at least D. The XX deficit is nonnegative on open intervals and even rings. For fixed smooth positive periodic w_i = W(i/L), D/Q -> lambda^2/(2 delta). Therefore liminf Def/(M Q) >= 1/(300 delta), unbounded at fixed d and M. For quadratic modular weights, D = (lambda^2/delta) log N + O_delta(1), Q = 2 log N + O(1). Take length to infinity first at each fixed delta.

[Proved] This also invalidates the general identification of the modular log coefficient with A alone. A gapped local gradient response contributes at the same order, even in an optimizing energy-density gauge. A corrected representative/subtraction prescription is not supplied here.

[Proved] C105 already excluded positive decompositions into arbitrary-length overlapping uniform intervals, not merely equal-length ones. Varying the lengths does not evade it. Trial-state dilation has the wrong direction for an upper deficit bound; inverse-gap curvature is explicit in this counterexample.

[Numerical] checks.py and checks_results.json contain passing small checks of the counterexample, XX, Heisenberg, XXZ(.5), SU(3), and four seeded random bonds. Full Hilbert spaces, maximum dimension 729, runtime .103 seconds, maximum residual 3.08e-12. Smooth scalar evaluations use at most 65536 sites, without a large diagonalization. NUMERICS.md gives the values and limitations. Random bonds use raw half-spread, not optimized M.

[Proved, infrastructure record] Typhon connection timed out during banner exchange. No SLURM job was submitted. large_gaussian.py and cluster_job.sbatch are prepared but unrun. No large computation ran locally or on a login node.

[Proved, file boundary] All edits are inside work_C111. Parent STATE.md and REPORT.md were not changed under the task-specific restriction. The next integration step is to report the counterexample and correct the unrestricted smooth-modulation and modular-matching claims, while leaving the original A bound open.
