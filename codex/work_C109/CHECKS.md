# C109 checks

[Numerical] All calculations are small deterministic normalization checks. The largest diagonalization is 64 by 64. They do not establish infrared matching.

[Numerical] PASS. XY two-site spectrum and swap-gauge identity, maximum eigenvalue error 0.000e+00; no tested gauge lowered the spread. Small-k slope relative error 1.250e-09.

[Numerical] PASS. Periodic spin-ring ground energies versus the independently evaluated Jordan-Wigner vacuum sum at L=4,6 and three anisotropies, maximum error 9.770e-15.

[Numerical] PASS. Ising Z2 seam energy delta E=2 tan(pi/(4L)), maximum error 9.326e-15. This gives lim L delta E=pi/2, consistent with v=2 and x=1/8.

[Numerical] PASS. Cyclic replica seam versus a single long ring, entire spectrum maximum error 5.773e-15; untwisted product ground-energy error 2.665e-15. A generic complex, oriented bond was used.

[Numerical] PASS. Parafermion twist identities verified using rational arithmetic for N=2 through 200, including the maximal-twist factor 8 and averaged factor below 12.

[Numerical] PASS. Positive covariance countermodel: finite Gram matrices positive and Fourier-transform integration agrees with the analytic formula to 2.243e-14.

[Proved, analytic conclusions checked above] The tests neither prove a uniform seam budget nor a general bound on A_neutral. The velocity-floor counterexample and replica identity are established analytically in PROOF.md.
