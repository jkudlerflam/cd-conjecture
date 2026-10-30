# C119 result

2026-10-29.

[Proved] The requested free families admit a uniform estimate, including their full parity-resolved spin traces. For even L,

    R_L(beta) <= 3 log d + (10 pi/3) M beta/L.                 (1)

This covers XX, arbitrary independent stacks of XX with different velocities, the C106 two-band family throughout its stated parameter range, and its massive a=0 comparison. Thus k=10pi/3 and f(d)=1 work in these families, with no crossover scale. At beta=sL/M their ratios cannot diverge, even when the Hamiltonian parameters depend on L. The coefficient is not sharp. [PROOF.md](PROOF.md), section 6, gives the proof from a matrix trace-convexity and quadrature estimate.

[Proved] Equations (1)-(5) of PROOF.md compute the complete finite-size XX ratio and its first L^-2 correction. For the scaling variable y=beta v/L,

    r_XX(y) = 2 log Z_XX(y)-log Z_XX(y/2),
    r_XX(y) = 4 exp[-pi/(2y)]-6 exp[-pi/y]+...       at small y,
    r_XX(y) = pi y/4-2 exp[-pi y/4]+6 exp[-pi y/2]+... at large y.

Independent spin stacks add these expressions. A slow XX sector contributes little at fixed beta/L, but contributes r_XX(s) at beta=sL/v, independent of v. At its own inverse-velocity scale it is not negligible.

[Proved] C106 has one shared fermion-parity projection, so its torus trace is not the product of two spin XX traces. Its exact finite products, compatible-length CFT function, and curvature correction are in PROOF.md section 4. The correction coefficients diverge as the avoided-crossing scale closes, while the exact ratio remains bounded by (1). Generic Fermi phases must be retained; arbitrary even lengths need not converge to a vacuum torus.

[Proved for free massive modes] At fixed positive gap and beta=tau L, the massive ratio tends to zero. The proof and the Bessel-function vacuum correction are in section 5. Modes with mL/v of order one still contribute before the infrared limit. Finite vacuum multiplicity instead leaves a constant log g.

[Proved within the stated CFT description] A large ratio is not by itself a count of fast sectors. A compact boson at fixed y and large Luttinger parameter has

    r_K(y) = (1/2)log K + O_y(1),

with c=1. Light primary states cause this growth. The microscopic XXZ attempt near the ferromagnetic endpoint does not violate the target: K grows while v falls, Kv tends to pi/4, and R at fixed beta/L has the finite theta-function limit in PROOF.md equation (21). It also shows that slow-sector suppression is not uniform over changing conformal spectra.

[Open] I did not construct a fixed-d, fixed-M interacting chain with R_L(sL/M) unbounded, and did not prove that such a construction is impossible. The free-family obstruction is rigorous. The interacting construction problem and the general bound A<=C M f(d) remain open.

[Proved and heuristic, distinguished] The sharp free asymptotic slope uses k=pi/2, f(d)=1, by the inherited A<=2M theorem and saturation by XX stacks. For general interacting chains, the strongest specified solvable-data envelope is the conjectural ordered-slope estimate with k=pi^2/4, f(d)=H_(d-1). Heisenberg fixes its d=2 coefficient. No claim that H_(d-1) is optimal at larger d is justified. This coefficient cannot be used in an all-length estimate without slope corrections, as C116 already showed.

[Assumed sufficient missing lemma] The remaining input is a lower bound on the normalized Euclidean reconnection insertion, or equivalently an upper bound on its integrated seam response, of order M beta f(d)/L. PROOF.md equations (23)-(25) state it with the constants and an allowed ordered o(beta/L) remainder. A positive-transfer formulation is an effective Renyi-2 rank bound at separation L. Gate rank and the norm of separate seams do not supply it.

[Numerical] The complete C106 trace agrees with direct spin diagonalization to 4.98e-14. Subleading formulas, slow and massive limits, compactification growth and 360 free-family scan points pass their checks. See [CHECKS.md](CHECKS.md). The small run took less than one second; its largest matrix has dimension 256.

[Work record] All requested files and a local STATE.md are in work_C119. No files outside this folder were changed. The extended scan is SLURM-ready in [run_typhon.slurm](run_typhon.slurm), but was not submitted because the sandbox cannot reach Typhon. No general interacting numerical search or large computation was performed.
