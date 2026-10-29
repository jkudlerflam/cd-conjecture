# C95 handoff

[Status] C92 equation (26) remains open, including at infinite grammar penalty. Read RESULT.md and PROOF.md. Do not describe the bounds below as a proof of the full d=13 mobile defect comparison.

[Proved] For the free canonical gap kernel and 1/3<=rho=N/(N+v)<=2/3, 0<=sum_r chi_N(r)<=48/t uniformly in N. Proof uses an exact isometry on the cyclic-invariant gap subspace, total contact C=sum_x n_x n_(x+1), and free-fermion two-particle-hole matrix elements (2/L)[cos(p1-k1)-cos(p2-k1)]. This controls only the signed spatial sum, not K_N.

[Proved] Free e(q)=-2t sin(pi rho)/(pi rho), a(q)=rho-sin^2(pi rho)/(pi^2 rho), q=(1-rho)/rho. The new derivative identities are e''=2pi t rho^3 sin(pi rho), a_q=-rho^2+(2rho/pi)sin(pi rho)cos(pi rho)-sin^2(pi rho)/pi^2. The thermodynamic one-sided density-mode coefficient is c_dens=a_q^2/(2e''). At half filling it is 0.07857579749228.../t.

[Proved] At positive T, c_dens=rho D a_rho^2/2, where D=d rho/d mu for the physical free fermion gas and its pressure equals the length tilt s. Uniformly for t/2<=s<=3t/4, t c_dens<=4 at sufficiently small T. This does not identify a finite-N kernel decomposition with bounded K_N.

[Proved] Adjacent interior vacancy-sector contact jumps are O(1/N), so the first-order crossing residual identified by C94 is O(JM) for every material pattern. Only that term is bounded. The O_(N,sigma)(J^2) remainder is still not uniform.

[Proved diagnostic] Gapped bath H0=t tau_z, A=(1+tau_x)/2, h=diag(0,1), frozen contact 1/2. Its exact residual has expansion -tanh(t/T)J^2/(16t)+tanh(t/T)J^3/(64tT)+O_T(J^4), but -J^2/(8t)<=R<=0 uniformly in T. Replicated systems retain the nonuniform Taylor coefficient. This is not a C88 counterexample.

[Open] The exact tagged-contact spatial first absolute moment and a resummed W-relative thermal bound for all legal auxiliary states, nonzero port, and unbounded wrong cells. A bare J=0 thermal Taylor coefficient uses normalized material trace and cannot be substituted for the frozen ground-state coefficient. PROOF equation (34) gives the exact interaction-picture formulation retaining Ja_0 Hmat in the reference.

[Checks] python3 -B checks.py passes. It performs small quadratures and finite form-factor sums, not large diagonalizations or a full d=13 simulation. All files and execution stayed here. No network, agents, cluster jobs, or external edits. The project STATE.md was not changed under the explicit folder restriction.
