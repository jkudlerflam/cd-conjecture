# C94 handoff

[Proved] C92's literal pointwise response bound (29) is false before the thermodynamic limit. At u=0 and a crossing between two conserved vacancy-sector ground states, D_sigma-D_good has coefficient (Delta a)^2[tr_normalized(Hmat_sigma^2)-tr_normalized(Hmat_good^2)]/(4T). The material subtraction does not cancel it.

[Proved, full C88] N=15, good=(15), wrong=(3,3,3,3,3), W=15. At s/t=0.610116148316937... the free sectors v=15,16 are the only minima. The coefficient T(D_wrong-D_good)>534.73. crossing_certificate.py certifies the crossing with exact rational interval arithmetic. The internal moment bound uses all legal auxiliary states and includes the scheduled port. The good one-marker norm is below 1016; the wrong normalized mean is at least 30000/13.

[Open] This does not disprove C92's integrated estimate (26), with its global o(N) term. Even its finite-N version without that term fails at a free crossing because the integrated residual can be linear in J. The thermodynamic order of limits matters.

[Proved] At fixed vacancy number, an all-order Schur-complement bound gives Delta R >= -2J^2 N^2 M^2/Delta for 4JNM<=Delta. It has no internal-gap assumption but is not volume-uniform.

[Proved conditional] For zero-port product clock ground states and chi=k-c/N, |Delta B|<=4M^2 b W(K_1+|c|), K_1=sum(1+dist)|k|. Arbitrarily long wrong clocks are allowed. A uniform free-kernel bound and the uniform relative higher-order remainder remain missing. This does not include the coupled auxiliary thermal trace.

[Numerical] 50 canonical comparisons, 10 thermal comparisons, maximum joint dimensions 13365 and 880. Largest eigenpair residual 1.29e-11. All canonical defect costs positive. Unfavorable residuals included, observed epsilon decreases with J in both a b15 high-density quantum example and a small-clock thermal example. Separate half-density second-order calculations reach N=v=9, without forming the joint interacting matrix. No half-density b15 many-cell simulation.

[Files] Read RESULT.md and PROOF.md. CHECKS.md gives tables. Scripts and JSON outputs are local. Project-level STATE.md was not edited because the user prohibited external writes. No network, agents or cluster jobs.
