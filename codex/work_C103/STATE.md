# C103 local state

2026-10-29. Completed investigation. The general velocity-weighted bound remains open.

[Proved] New obstruction 1. A nearest-neighbor d=8 quasi-free chain obtained by blocking hopping at distances one and three has A_ring(L)=A+7pi^2(t1-27t3)/(810 L^2)+O(L^-4), A=2(t1-3t3)/3. For t1/27<t3<t1/3 it approaches A from below. This refutes raw Casimir monotonicity toward the infrared, not the conjecture or every endpoint c-theorem.

[Proved] New obstruction 2. For the XX chain with hopping t, the normalized thermal current variance is R_J(T)=2t-7pi^2 T^2/(20t)+O(T^4/t^3). It increases toward the infrared and tends to zero at high temperature. The current is conserved. This refutes raw current-coefficient monotonicity, not A<=2M.

[Proved conditional implication] For a complete set of conserved-current sectors with D=sum_a v_a G_a/pi and W>=0 satisfying tr(WG_a)>=c_a, A<=pi tr[W(K-R)]. A microscopic upper bound on this weighted optical budget would prove the target. Retaining regular weight gives A<=pi[tr(WK)-2s0^2/s1] through Cauchy-Schwarz. No general d,M estimate of that expression is established.

[Proved from inherited exact data] At Heisenberg, the spin twist has K=(2/3)(ln2-1/4), D=1/4, W=2, so WK=0.5908629074 exceeds M=1/2. Discarding regular weight cannot give the sharp coefficient. Neutral sectors are the other obstruction to extending the conserved-current route.

[Numerical] python3 checks.py passes. CHECKS.md records blocked dispersion and ring formulas, independent XX thermal integrals, inherited C101 ratios, and finite-ring optical moment and stiffness identities. Computations are small, with no cluster work. The initial command name python was unavailable; python3 ran successfully.

[Work record] Deliverables are RESULT.md, PROOF.md, PROGRESS.md, checks.py and CHECKS.md. All writes stayed within work_C103. The parent project memory is intentionally unchanged under the task's folder restriction. Next work would need an independently bounded response including neutral sectors, or a genuine weighted c-function with a lattice ultraviolet budget. Do not repeat the raw monotonicity proposals or remove a matching scale by setting it to one.
