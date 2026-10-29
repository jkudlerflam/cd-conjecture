# C90 progress

[Scope] All writes and any Python execution stay in work_C90. No network or delegated agents. The parent STATE.md is read but will not be edited under the task's explicit restriction.

[Read] Project memory, the two S background files, REPORT sections 9.25 through 9.51, C48/C58/C51 proofs and handoffs, both requested referee reports, and C88's full proof and handoff. Additional dependency checks cover C81/C83/C85, the d=24 frozen matching and vacancy audits, and the corrected N and TP statements, including C79's progress.

[Proved, scope finding] C48's abstract theorem and C51's Gaussian identification are independent of the alphabet size. C58's pinning, bounded seam, pressure-zero and canonical-duality arguments are also reusable. C58's invariant illegal-cell flags and full microscopic verification are specific to the wrapper. They are not verified for the colored mobile bond by changing d0 to 13.

[Open, identified] The vacancy fluid must remain among the gapless fields. Its phase-model stiffness is K_gap=4 at half occupancy, whereas the physical particle normalization has K_particle=1. The finite-J matching into this reference is explicitly open in the supplied vacancy notes. The small alphabets also have odd-charge auxiliary states for which exp(i pi Q/2) is not a fermion parity. The existing separate-Jordan-Wigner lemma does not by itself fix this.

[Completed] PROOF.md supplies the reusable algebra and thermodynamic implications, finite-penalty budget, coordinate factors, and a conditional theorem with initialization, defect-pressure, RG/source, scalar and completeness obligations stated separately. RESULT.md summarizes what transfers and what does not. Local STATE.md records the handoff.

[Numerical, passed] python3 -B checks.py checks 54 fixed-density pinning inequalities and the compression identity, gap-energy curvature at five densities, independent free-fermion thermal quadrature, the coordinate conversions and harmonic margin, and synthetic finite-penalty budgets. Largest matrix dimension 64. Compression residual 1.11e-16. At half density the curvature calculation gives K_gap=4.000000001001962. The thermal coefficient approaches pi/(12t), with relative error 1.70e-5 at T=0.01 and t=1.3. All assertions pass. No C88 certificate rerun or microscopic infrared computation.

[Remaining research inputs] I, D_E, D_T, N_mob, TP_mob and G_mob are named in the final proof. They can be grouped as expanded N, TP and dressed-defect comparison only with the extra initialization, thermal and completeness clauses stated explicitly. No claim that these inputs have been established for one numerical finite-J member.
