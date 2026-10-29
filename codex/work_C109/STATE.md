# C109 local memory

[Proved, record] Completed 2026-10-29. Task was a crossover-free bound on A_neutral=sum c_a v_a for sectors invisible to conserved on-site charges. Read the full background and REPORT 9.68 through 9.77 first. The parent project STATE.md is outside the allowed write folder and was left unchanged.

[Proved, outcome] No general C M f(d) bound or counterexample. Single bounded-c factors are controlled by v<=4 e M, conditional on a nontrivial local dynamical scaling limit. Ising cv<=2 e M; unitary minimal models cv<4 e M; standard Z_N parafermions cv<8 e M. Uncontrolled multiplicity prevents summing this into a dimension-only bound.

[Proved, discrete route] A symmetry seam changes the energy by at most 2M, whereas its CFT shift is 2 pi sum v_a x_a(g)/L. Missing: an upper bound on limsup L delta E and sector coverage. Standard parafermion disorder dimensions cover c with factor 8 for the best twist or factor below 12 for the group average. Replica cyclic twist is exactly E0(nL)-nE0(L), so its universal version reduces to C98.

[Proved, current route] Equal-space energy-density tail measures B2=sum c v^2; energy-current tail measures B4=sum c v^4. Local norms do not bound those asymptotic coefficients. A positive covariance countermodel establishes the logical gap. Even a moment bound needs integrable counting at small velocities.

[Proved, new exact obstruction] Critical anisotropic XY at field one has h=-(1+gamma)XX/2-(1-gamma)YY/2-(ZI+IZ)/2. Bond eigenvalues are +/-sqrt(1+gamma^2), +/-1. Swap symmetry proves M=sqrt(1+gamma^2) is gauge optimal. For every gamma>0, c=1/2 and v=2gamma. Hence no positive v/M floor exists even at d=2 with one neutral Ising sector.

[Assumed, sufficient lift] Need d'<=g(d), M'<=K(d)M, and A'_WZW>=eta(d)A_neutral in the same units, plus C104's microscopic-current matching. Then A_neutral<=[4 pi K/(3 eta)]g(g^2-1)M. No general construction is supplied. A U(1)-only lift would still need the missing stiffness-weight comparison.

[Numerical] checks.py completed successfully. Largest diagonalization 64 by 64. No heavy local numerics, cluster jobs, network shell commands or writes outside this directory. Results and limitations are in RESULT.md; complete derivations and sources in PROOF.md; check output in CHECKS.md.
