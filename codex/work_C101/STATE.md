# C101 state

2026-10-29. Main deliverables RESULT.md and PROOF.md. The full requested all-integrable theorem is not established. There is no violation in completed checks.

New results, conditional on standard complete continuum identifications:
- Balanced repulsive Hubbard at all metallic fillings: optical charge/spin stiffness and K_c>=1/2 give A<=6t sin(pi n/2). Exact bond M>=2t proves Q<=3/(2pi). Attractive balanced branch with spin gap and K_c>=1 gives Q<=1/(2pi), with that Bethe input explicit.
- SU(q) symmetric supersymmetric t-J: exact bond M is the half-range of {0,-mu/2+/-t,-mu,-mu-2t}, q>=2. Weighted stiffness gives A<=2pi t for q>=3. For q=2, SU(2) averaging and the exact four-site certificate give A<=(3/2)pi t. K_c>=1/q follows from Kawakami's exponent range and the SU(q) LL vertex formula. General grading or flavor-split phases are not covered.
- Exact rational certificate: tj_certificate.py proves 27I-6 sum_{j=1}^3 Bbar_j/t positive definite in 15 charge blocks totaling dimension 81, Bbar/t=1.5T-(8/3)S.S. Fractions only for certification. Output tj_certificate.json.
- Sp(2r) M=1 for r>=2, M=1/2 for r=1. Martins's sector assignment yields A/pi=2(r-1)/r+1/(r+1), so the bound follows conditionally. Published completeness evidence is finite-size, not an all-rank proof.
- General graded permutation with fields has eigenvalues J(-1)^p_a+f_a and (f_a+f_b)/2+/-J. Ordinary reflection proves exact gauge optimality. General A unresolved.
- XYZ width is (largest two |J|)/4 in spin normalization. Critical zero-field loci reduce to XXZ.

Inherited C35 results are cited rather than re-derived: SU(N) symmetric all rank/level, SO(N), XXZ, product closure, and the symmetric-family Theta(log d/log log d) envelope. Heisenberg saturates, XXZ approaches it; SU(3) is only the discrete symmetric-family runner-up.

Numerics: bond spectra errors <8e-15. Hubbard 66-point scan maximum Q=0.3168709574, final quadrature change <=1.01e-6, half-filled Bessel spin-velocity checks <=4.2e-15. Charge velocity at the maximum is 2.0712t, so a purported v_c<=2t proof is invalid. Ladder screen 14 couplings, L=4,6,8, sectors <=4900; maximum two-size Q=.8745965622 at SU(4) (exact .75), residual <=1.15e-11. Small-ring numbers do not establish conformality or unfamiliar full ground states.

SSH to Typhon timed out during banner exchange. No job submitted and no login-node computation. Larger ED is prepared in typhon_ed.slurm, execution guarded by SLURM_JOB_ID. Work was confined to this folder, so project STATE.md was not edited.
