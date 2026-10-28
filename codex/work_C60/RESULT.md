# C60 result

[proved, counterexample to a specified source majorant] C55 (30a) fails with its overall interaction factor epsilon when the comparison is to the local continuum pair primary and the finite-crossover term comes from finitely many smooth integrations at fixed scale. At zero interaction, for P_x=f_x f_(x+1),

    <P_x P_0^dagger>
      = sin^2(p_F)/(pi^2 x^2)
        + [sin^2(p_F)-sin^2(p_F x)]/[pi^2 x^2(x^2-1)],
        x>=2.

The first term has the matched free continuum amplitude and X_pair=1. At p_F=pi/3 and x=3m, the second term is

    3/[4pi^2 x^2(x^2-1)],

which is nonzero although epsilon=0. A finite smooth crossover changes the correlation only by rapidly decreasing terms, so it cannot remove this x^(-4) tail. The exact calculation and cutoff argument are in [PROOF.md, Sections 5-6](PROOF.md).

[scope of the counterexample] This refutes the interaction-prefactored continuum-source remainder in that convention. It does not refute the physical Luttinger limit, the source-free normalized endpoint majorant, or a comparison retaining the full free lattice source problem. If “finite-crossover term” means an arbitrary algebraic remainder with only one exponent M<4, the displayed example can be hidden in that term. Such a term is not a finite smooth crossover. This distinction is explicit in the proof.

[necessary correction] Bare pair-source Taylor terms and fixed free curvature have nonzero size at the Gaussian point. They need their own remainder budget, or an explicit free subtraction, before an interaction-proportional source estimate is possible. The collaborator's IR_SOURCE_REVIEW already distinguishes this budget from interaction-induced source corrections. No corrected nonlinear theorem is claimed here.

[proved, finite-initialization result] The matched action satisfies the degree hierarchy

    ||K'_(2m)|| <= C^(2m) epsilon^max(1,m-1)

after the fixed finite integrations, with fixed norm reserves. For a connected term with n vertices and c contractions, c>=n-1 and m=sum_i m_i-c; hence its artificial analytic order is at least max(1,m-1). Cauchy estimates give the displayed bound uniformly in degree. This is proved in PROOF Section 3 using the imported finite-integration theorem.

[read, BM route] [BM 0907.2837, Section 3](https://arxiv.org/html/0907.2837) permits time-nonlocal, all-degree kernels at scale zero. Its displayed degree bound is therefore compatible with the matched action. It is not a stand-alone replacement-kernel hypothesis. The referenced tree proof also uses localization, position moments, and simultaneous derivative/zero scale allocation. Matching the degree powers does not make the induction applicable verbatim.

[unresolved, action majorant] I did not establish or refute the newest source-free endpoint lemma. PROOF (17) records a concrete scale-allocation inequality for its proposed routing, to be proved in the same expansion that controls edge loads and anchor integration. C55's all-scale linear theorem covers the one-endpoint sector, not this multi-vertex inequality. Full Lemma N also requires the source extension and a corrected remainder statement.

[numerical] python3 -B checks.py passes. The free Wick formula agrees with its closed form to 3.47e-18; the fixed-cutoff absorption identity agrees to 4.97e-16; 35,014 finite degree-count cases pass. These are controls of exact formulas, not an interacting RG calculation.

[record] PROOF.md, PROGRESS.md, STATE.md, checks.py and checks.json are local to work_C60. No subagents or large computations were used. Parent project memory was not modified because the task restricts writes to this folder. No independent review or central-charge conclusion is claimed.
