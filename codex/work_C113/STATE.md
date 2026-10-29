# C113 session memory

2026-10-29

[Proved] The smooth-modulation lemma is false even for the unique representative minimizing the bond spread. Section 9.1 of PROOF.md is the strongest result, a d=4 embedding with M=8. The desired bound on A remains open. A different infrared-adapted representative with controlled but nonminimal norm is not excluded.

[Proved within CFT] With A = (1/2) sum_chiral c_alpha v_alpha, the stress coefficient is kappa = 1/(24 pi). Circle formula relative to the uniform finite-circle ground state includes pi A/6 [1/ell_w-W/L^2], ell_w = integral dx/w. Take infinitely many modulation periods before wavelength R -> infinity to eliminate it. Smooth weights are positive and bounded away from zero.

[Assumed, needed for equivalence] Pure source matching P means no extra two-gradient identity coefficient and no gradient-induced operator contributing at order R^-1 or larger. Under P, Def/G -> A/(24 pi), giving both directions with factor 24 pi. Uniform CFT spectral or Casimir matching alone is insufficient. Renormalized irrelevant infrared contributions are suppressed, but short-distance analytic response need not be.

[Proved] Gapped example h = -g(ZI+IZ)/2 + lambda(XI-IX), d=2, M=g, gap=2g. Def = sum [sqrt(g^2 a_i^2+lambda^2 b_i^2)-g a_i], a_i=(w_i+w_(i-1))/2, b_i=w_i-w_(i-1). Def/G -> lambda^2/(2g), arbitrary at fixed M and gap. Tensor with XX: d=4, A=1, M=g+1/2, arbitrary extra coefficient.

[Proved] Fixed-bond z=2 violation: h_i=-J S_i.S_(i+1)+lambda(S_i^z-S_(i+1)^z), M=J/2. For w_i=1+epsilon cos[2pi(i+1/2)/R], coherent product spins n_i=(cos(2pi i/R),0,sin(2pi i/R)) give Def_R >= lambda epsilon R sin(pi/R)/2 - JR[1-cos(2pi/R)]/4 -> lambda epsilon pi/2. G_R=O(1/R). Canonical ferromagnet without the divergence instead has Def=0. Alternating z rotations map to XXZ Delta=-1.

[Proved within XXZ Luttinger matching] Longitudinal divergence adds chi lambda^2/2 to Def/G, chi=K/(pi v), diverging as lambda^2/[4(1+Delta)] near -1. First take R -> infinity at fixed Delta. This does not make A diverge.

[Proved, conditional] Finite-amplitude gap bound: if H_s=H+s(H_w-H) has a unique ground state with gap Delta_* and current covariance Schur bound K_*, then Def <= K_*/Delta_*^3 sum(delta w)^2. Exact identity i[H_s,V]=sum(delta w) i[h_i,h_(i+1)] is independent of s.

[Proved, conditional sufficient target] A bound on the small-q symmetrized inverse-cube current moment implies A <= 24 pi C M f(d) given spectral stress matching S. It allows extra positive spectral response and is weaker than assuming pure weighted-energy matching. No universal upper bound was proved; it fails for the unique minimum-spread prescription in the constructed families.

[Numerical] python3 checks.py passed. Arrays <=4096, matrices at most 16 by 16. No heavy computations or remote jobs. checks.json records gapped ratio approaching 0.245 and the ferromagnetic ratio increasing past 271 at R=2048. Initial python command unavailable, then python3 succeeded.

[Proved, file map] RESULT.md is the summary; PROOF.md gives definitions, derivations, exact witnesses, assumptions and primary links; PROGRESS.md records work; checks.py and checks.json contain small verification. All files are in work_C113. The parent project memory and shared reports were left untouched as required.

[Proved, strongest counterexample] Site=C^2_A direct-sum C^2_B. h_AA=-S.S+SzI-ISz; h_BB=16 S.S+8I; h_AB=h_BA=12I. BB extrema -4,12 witness M=8. Every optimal gauge must annihilate the BB singlet, forcing off-diagonal gauge blocks zero and o_B scalar. The two mixed upper-endpoint blocks then force o_A=o_B, so the gauge divergence is zero. Uniform full-chain ground space is all-A ferromagnetic: AA domain energy >=-n_AA/4-1; BB domain energy >=-4 since each two adjacent BB bonds sum to a positive operator; mixed bonds cost 12. With m interfaces, E>=-n_AA/4+19m/2. The all-A coherent-state witness proves Def/G diverges in this unique minimum-spread representative.

[Proved, strongest gapped variant] Replace AA by -g(ZI+IZ)/2+XI-IX, 0<g<=1. Same unique minimizing gauge, M=8, all-A all-up unique uniform ground state with gap2g. Weighted all-A trial energy gives liminf Def/G>=1/(2g), disproving uniform constants at fixed d=4,M=8.

[Proved within XXZ infrared dictionary] AA can instead be critical XXZ plus SzI-ISz. Active open energy Eopen(s)>=s eA-2 follows from product tiling, so interfaces remain costly and the uniform infrared is XXZ alone. Its susceptibility response diverges near Delta=-1 at fixed M=8. This disproves the minimum-spread converse within the CFT class at the stated physics standard.
