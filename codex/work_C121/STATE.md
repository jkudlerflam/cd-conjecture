# C121 handoff

2026-10-29.

[Proved] For the positive midpoint measurement-channel trace Q, N d-outcome seam pinchings at any spacing satisfy d^(-N) Z <= Q <= Z. Lower bound: D(X) >= X/d propagated through completely positive maps. Upper bound: expand pinching into d unitary conjugations and use Schatten Holder on each imaginary-time trajectory. This is independent of locality, gaps and Markov properties.

[Proved] At N = max(1, ceil(beta M/(a L))), the pinching loss is at most (1+beta M/(a L)) log d per seam. Both coherent and positive measurement histories admit the d^N Cauchy-Schwarz step. Pinching is a reflection-paired positive channel, distinct from inserting resolutions of identity into the cyclic Euclidean trace.

[Disproved] No-cost spatial reconnection of these coarse histories fails, including stoquastic XX and Heisenberg. For positive histories define G = log(sum q_s^2/Q_D), with both doubled-ring seams pinched. Exact reduction R <= 3N log d+G. At tau = L/M, N = 1, G is 0.762551 (XX L4), 1.888216 (Heisenberg L4), 2.292663 (SU3 L3). Large XX L8192 beta4L violates even sum q_s^2 <= Z_(2L). Analytic XX witness at L4 beta64 is also proved.

[Proved] Continuous pinching tends to the Gibbs trace for D(H). XX loses both seam bonds and has E_Zeno-E_L = 1/2+(1/2)tan(pi/(2L)) -> M. The cheap coarse measurement leaves quantum transmission between samples. Nonnegative microscopic weights do not eliminate those unobserved histories.

[Proved] R-3N log d <= G <= R+2N log d. The remaining lemma is comparable in strength to the target, so the pinching theorem alone does not bound the infrared slope.

[Open] Need G <= b_0(d)+b_1(d) N with explicit universal constants and a. Then A <= 4[3 log d+b_1(d)]M/(pi a). No general ratio or velocity-weighted bound proved. Zero-cost failure at a = 1 does not exclude smaller a or a finite cost per slice.

[Numerical] Full small ED through qubit L8 and SU3 L6. Forty-eight channel cases, history enumeration for N<=4, lower gluing witnesses for N8, eighteen complex/commuting algebraic controls, fixed-beta Zeno tests. Exact parity-retaining one-pinch XX formula agrees with ED to 3.22e-15 and reaches L8192, paired thermal length16384. No large interacting run, large multiple-pinch determinant, network shell command or independent review.

[Work record] RESULT.md summarizes results, PROOF.md contains derivations and the missing lemma, CHECKS.md records numerical scope and failed claims, checks.py reproduces checks.json and summary.json. PROGRESS.md records the work. The parent STATE.md was read but remains unmodified under the task-specific folder restriction.
