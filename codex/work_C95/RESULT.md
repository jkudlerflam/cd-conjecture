# C95 result

[Status] C92 equation (26) remains unproved. The full mobile defect comparison for the d=13 chain at infinite grammar penalty cannot be concluded from this work. The exact results below control the uniform contact response and the conserved-density contribution, but not the spatial first absolute moment or the full quantum thermal remainder.

[Proved, exact free-fermion calculation] At fixed vacancy number and material density 1/3<=N/(N+v)<=2/3,

    0 <= sum_r chi_N(r) <= 48/t.

The cyclic-invariant gap subspace maps isometrically to the translation-invariant physical hard-core chain. The total contact becomes sum_x n_x n_(x+1). Its exact two-particle-hole form factors give the bound without a volume-independent bath gap. This is a signed spatial sum. It does not bound C94's K_N=sum_r(1+dist)|k_N(r)|.

[Proved, thermodynamic density contribution] The one-sided response gained by allowing vacancy number to relax is

    c_dens(T,s) = a_q(T,q)^2/[2 f_qq(T,q)],
    q=v/N.

For t/2<=s<=3t/4, t c_dens<=4 at sufficiently small positive T. At half filling and T=0,

    c_dens = (2/(pi t))(1/4+1/pi^2)^2
           = 0.07857579749228.../t.

This does not establish that subtracting c_dens/N leaves a kernel with bounded first absolute moment.

[Proved, finite-volume crossing term] At adjacent free vacancy crossings in an interior-density window,

    |a_N(v)-a_N(v+1)| <= C/N.

Since every material ground energy has magnitude at most MN, C94's explicitly identified linear crossing residual is O(JM), uniformly in the pattern. That term is global and subextensive. The higher-order finite-volume remainder is not bounded by this argument.

[Proved limitation of the proposed expansion] Frozen subtraction does not automatically produce a Taylor series with a temperature-uniform remainder. An exact two-level bath with a bounded material source gives

    R(J,T) = -tanh(t/T) J^2/(16t)
             +tanh(t/T) J^3/(64tT)+O_T(J^4).

The cubic coefficient diverges as T->0, even with a gapped bath and no vacancy crossings. Independent copies retain the obstruction in the thermodynamic limit. Nevertheless -J^2/(8t)<=R<=0 uniformly in T. This diagnostic does not refute equation (26). It shows why a resummed estimate can be needed even when the integrated correction is small.

[Open] Two estimates are still missing: (i) the first absolute spatial moment of the tagged-contact resolvent, after the canonical constant term is removed; (ii) a comparison of the full resummed thermal residual proportional to W, including auxiliary states, the port, and unbounded wrong-cell lengths. The exact interaction-picture target is PROOF.md equation (34). A bound on separate extensive residuals cannot replace their relative comparison.

[Numerical] All checks in python3 -B checks.py pass. They include the free-fermion form-factor sums, agreement with three independent C94 gap-resolvent rows, finite-size crossing scaling, thermal quadratures, and the exact two-level diagnostic. No full interacting d=13 calculation was performed.

[Files and scope] PROOF.md contains the derivations. CHECKS.md records the checks. PROGRESS.md and STATE.md provide the handoff. All new files and Python execution stayed in work_C95. No network commands, delegated agents, cluster jobs, or external edits were used.
