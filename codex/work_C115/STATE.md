# C115 local state

2026-10-29.

[Open] General A<=C M f(d) is not proved. No counterexample to the target. Work completed as a negative result for unfiltered finite-size counting, with explicit sufficient velocity-window lemmas.

[Proved] The proposed tower-window direction is reversed. Delta E=2 pi v Delta/L means an energy window 2 pi u x/L reaches Delta<=u x/v. Slow sectors dominate unfiltered low-energy counts.

[Proved] Fast-plus-slow XX stack: d=4, t_f=1/(1+epsilon), t_s=epsilon/(1+epsilon), M=1 by physical swap convexity, v_f=2t_f, v_s=2t_s, N(1)=1, A=2. After the infrared limit, excitation trace at beta=tau L is at least [1-exp(-2 pi tau v_s)]^-1 and is unbounded across epsilon. Finite O(1/L) window counts also diverge. This uses only finite-level convergence before summing, not uniform trace convergence.

[Proved] A spatial quantum transfer matrix carries the imaginary-time column, not a fixed d^2 state space. E/|P| does not isolate fast sectors because counterpropagating slow excitations cancel momentum.

[Assumed missing lemma; proved implication] With V=4eM, I_j=(V2^(-j-1),V2^(-j)], C_j=sum_(I_j)c_a, C_j<=K(d)2^(j/2) implies A<=4e(2+sqrt(2))M K(d). More generally C_j<=K(d)2^(gamma j), gamma<1, gives coefficient 4e/[1-2^(gamma-1)]. C_j<=K(d)(j+1) gives A<=16eM K(d).

[Proved conditional spectral extraction] For full modular window state content, Q_j=Tr exp[-pi sum(h+hbar)] has log Q_j>=pi C_j/4. A new lattice upper bound log Q_j<=B(d)2^(j/2) would give A<=16e(2+sqrt(2))M B(d)/pi. General spectral selection, modular completion for nonfactorized factors and the upper bound remain open. Full traces also depend on radii and are stronger objects than c. Unequal velocities require sector-resolved energy and conformal momentum, not an implicit reflection symmetry.

[Proved dynamical limitation] SIE plus even an ideal width-ell packet protocol gives only C_w<=const M ell log d/w. Independent window estimates are borderline nonintegrable. C106 excludes a universal IR wavepacket width at fixed v/M for its shrinking pockets. A new simultaneous multi-scale capacity input is still needed.

[Proved tests] Free N(u)<=min(log_2 d,2M/u). XX stacks saturate A=2M and satisfy the spectral estimate with B(d)=b_XX log_2 d, b_XX=0.868313632810 from the full boson dimension trace in the stated radius convention. Exponential suppression of the maximum payload velocity yields a logarithmic fast count. Fixed-alphabet comparable slow velocities can concentrate many sectors in a few far-infrared windows, so uniform C_j is too strong for that scenario.

[Proved gauge witness; inherited numerical input] The C22/spatial-cell vacuum hop X_(v,a) R_b -> L_a X_(v,b) has amplitude -t and changes both site letters. All coboundaries have zero such matrix element, proving optimized M>=t. The n=23 documented maximum cell velocity bounds v_max,physical/M by b v_max,cell/t, reported as approximately 24 times 10^(-205.2) lambda^2 g^2 alpha_n^4 after b=n+1 conversion. The printed values are not certified enclosures. This is an effective-spectrum test with the campaign's existing infrared qualifications, not a new all-n interacting proof.

[Numerical] checks.py passes, maximum matrix 16 by 16. Its first launch had a stray unary plus on a string in the report initializer; this was corrected before any mathematical check ran. No failed assertion remains. CHECKS.md records the final run.

[Proved file map] RESULT.md is the handoff; PROOF.md contains all derivations; PROGRESS.md is the work record. No writes outside work_C115, no heavy computations, no cluster use, no shell networking, no commits or pushes. The parent STATE.md and REPORT.md were not edited.
