# C113 result

2026-10-29

[Proved] The proposed smooth-modulation lemma is false, even for the unique bond representative attaining the minimum spectral spread M. An explicit d = 4 chain with M = 8 has a z = 2 ferromagnetic infrared and a positive limiting deficit under increasingly slow modulation, while the proposed quadratic-gradient allowance tends to zero. Gapped and conformal variants also obstruct a universal coefficient. These counterexamples do not violate the proposed bound on A.

[Proved within unitary diffeomorphism-covariant CFT] The universal stress contribution is

    Def_CFT(w) = A/(24 pi) integral w'^2/w,
    A = (1/2) sum_chiral sectors c_alpha v_alpha.              (1)

Thus kappa = 1/(24 pi). This includes arbitrary positive velocities, charged and neutral sectors, and unequal left and right data. No internal symmetry is needed. The formula follows by conjugating the weighted Hamiltonian with f'(x) proportional to 1/w(x); the Schwarzian is the curved-background anomaly. PROOF.md Sections 2 and 3 give every step.

[Proved within CFT] On a finite circle, the deficit relative to its uniform ground state also contains

    (pi A/6)[1/(integral dx/w) - (integral w dx)/L^2].         (2)

This term cannot be dropped at order 1/L. For example, it exactly cancels (1) for w = 1+epsilon cos(2pi x/L). To isolate (1), take infinitely many modulation periods before taking their wavelength to infinity.

[Proved, conditional equivalence] Fix a representative and assume its weighted lattice energy has no extra response at the two-gradient order, including no gradient-induced relevant source. Then

    lim Def(w_R)/G_R = A/(24 pi),
    G_R = sum_i (w_i-w_(i-1))^2/w_i,

for any fixed nonconstant positive smooth profile of wavelength R. Under this matching assumption,

    Lemma(C,f) implies A <= 24 pi C M f(d),
    A <= C' M f(d) implies
        Def(w_R) <= [C'/(24 pi)] M f(d) G_R + o(G_R).          (3)

The remainder and the wavelength threshold can depend on the bond. Uniform CFT finite-size scaling does not establish this stronger source-matching assumption.

[Proved] Infrared irrelevant operators of dimension Delta > 2 have renormalized long-distance corrections of order R^(1-Delta), under regular conformal perturbation theory. However, microscopic modes can supply an analytic term B_h integral w'^2/w at the same order as (1). A derivative of a dimension-one current also contributes at that order. Therefore the requested assertion that all non-CFT corrections are higher order is not valid without extra hypotheses.

[Proved, exact gapped counterexample] For Pauli matrices X,Z, take

    h = -(g/2)(Z tensor I + I tensor Z)
        + lambda(X tensor I - I tensor X),    g > 0.

The uniform chain consists of independent spins, with gap 2g, e0 = -g and gauge-optimized M = g. For a_i = (w_i+w_(i-1))/2 and b_i = w_i-w_(i-1),

    Def(w) = sum_i [sqrt(g^2 a_i^2+lambda^2 b_i^2)-g a_i],
    lim_(R -> infinity) Def(w_R)/G_R = lambda^2/(2g).         (4)

The coefficient is arbitrary at fixed d, M and gap. Tensoring this chain with XX gives a d = 4 chain with unchanged conformal coefficient A = 1 but the same arbitrary extra response. This disproves the unconditional converse even in the presence of a CFT sector.

[Proved, exact z = 2 counterexample] Consider

    h_i = -J S_i dot S_(i+1) + lambda(S_i^z-S_(i+1)^z),
    w_i = 1+epsilon cos[2pi(i+1/2)/R],

with J,lambda > 0 and 0 < epsilon < 1. The uniform chain is the ferromagnet, with M = J/2 and quadratic magnon dispersion. A product of coherent spins with n_i = (cos(2pi i/R),0,sin(2pi i/R)) proves

    Def_R >= lambda epsilon R sin(pi/R)/2
                - J R[1-cos(2pi/R)]/4.

The right side tends to lambda epsilon pi/2 > 0, whereas G_R tends to zero as 1/R. Thus the lemma fails for one fixed bond, even with its raw bond norm in place of M. Alternating pi rotations map this example to the ferromagnetic XXZ point Delta = -1. Removing the divergence restores Def = 0, so this does not refute a theorem that first prescribes a suitable representative.

[Proved within the XXZ Luttinger description] Approaching Delta = -1 from the critical side, the same longitudinal divergence contributes chi lambda^2/2, where chi = K/(pi v). This coefficient diverges as lambda^2/[4(1+Delta)] while A = v tends to zero. The limits must be ordered with wavelength tending to infinity first. The exact coherent-state witness handles the endpoint without Luttinger theory.

[Proved, failure survives unique gauge minimization] Give each site two spin-1/2 label sectors, A and B, and use the block-diagonal bond

    h_AA = -S_1 dot S_2 + S_1^z-S_2^z,
    h_BB = 16 S_1 dot S_2 + 8 I,
    h_AB = h_BA = 12 I.                                     (5)

Its spectrum lies in [-4,12]. The BB singlet and triplet force M >= 8 in every gauge, so M = 8. Saturation of the spectral endpoints forces every gauge operator preserving this optimum to be scalar. Thus the minimizing representative is unique. Positive interface costs force the uniform ground space to lie entirely in A, where it is the ferromagnet. The coherent-state witness above then proves failure for this minimizing representative. PROOF.md Section 9.1 gives the block spectra, gauge-uniqueness proof and domain-energy bounds explicitly.

[Proved, gapped minimum-spread variant] Replace h_AA in (5) by the independent-spin bond in (4), with lambda = 1 and 0 < g <= 1. The full chain has d = 4, M = 8, unique minimum-spread representative, unique uniform ground state and gap 2g. Its deficit satisfies liminf Def/G >= 1/(2g). Sending g to zero disproves a universal coefficient even within this gapped class.

[Proved within standard XXZ matching] Replacing AA by critical XXZ plus S_1^z-S_2^z similarly preserves the unique minimizing gauge and M = 8. The label sector remains massive. The lower bound Def/G >= chi/2 diverges near Delta = -1 although A = v tends to zero. Thus minimum-spread matching does not repair the converse for conformal infrared sectors either. The z = 2 and gapped variants above are exact lattice results and do not require this infrared calculation.

[Proved, gapped control with stated inputs] Along H_s = H+s(H_w-H), assume a unique ground state, gap at least Delta_*, and a connected bond-current covariance row sum at most K_*. Then

    Def(w) <= K_*/Delta_*^3 sum_i (w_(i+1)-w_i)^2
            <= w_max K_*/Delta_*^3 G(w).                    (6)

PROOF.md Section 8 derives this from the spectral formula for E'' and an exact commutator identity. It is a gap-dependent estimate, not a crossover-free bound in terms of M and d alone.

[Proved and heuristic, distinguished] The canonical ferromagnet and the hard-core-boson quadratic-band vacuum edge have Def = 0 for all positive weights, by local positivity. More generally, hyperscaling would give an infrared contribution of order R^-z, which is subleading to G_R for z > 1. That scaling argument requires regular response and excludes singular ground-space effects. The exponent z alone cannot control the deficit, as the explicit counterexample shows.

[Proved, conditional sufficient replacement] A weaker lattice target avoids asserting equality with the full deficit. If Sbar(q) is the symmetrized inverse-energy spectral moment of the bond density, continuity gives

    Sbar(q)/[4 sin^2(q/2)]
       = symmetrized sum_(m>0) |<m|j_q|0>|^2/(E_m-E0)^3.

A universal upper bound C M f(d) on its small-q limit implies A <= 24 pi C M f(d), provided its low-energy stress spectral weight matches the CFT. Positivity retains all additional microscopic contributions. This target needs a specified representative and is sufficient, not equivalent to a bound on A.

[Numerical] All small algebra and quadrature checks passed with python3 checks.py. At R = 2048 the exact gapped ratio is 0.2449998826 versus its limit 0.245, and the ferromagnetic witness gives Def/G >= 271.8032521. No heavy numerics, cluster submission, or many-body diagonalization was needed.

[Proved, completion status] The exact continuum reduction, its conditional converse, a gap-dependent bound, and explicit infrared edge cases are written in PROOF.md. The all-representative lemma is disproved. The minimum-spread repair is also disproved. A different, infrared-adapted representative with norm controlled by M is not excluded, and the general bound on A itself remains open. All edits are confined to work_C113; the project-level STATE.md was not modified under the task-specific restriction.
