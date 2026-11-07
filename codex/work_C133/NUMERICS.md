# C133 full-ring numerical results

[numerical] Ground energies come from complete d^L Hilbert spaces, without charge-sector restriction. Reflection symmetry proves the M used here is gauge optimal. Each Q uses the same bond on L and 2L. A_eff=4Q/pi is a finite-size estimator, not an established infrared A. The independent fit uses E_L/L=e_inf-pi A/(6 L^2)+b/L^4 where at least four compatible sizes are available. Three sizes give no fit diagnostic, so no extrapolation is reported for them.

| Family | Parameter | d | Doubling L | Q_L | A_eff/M | Fit A/M | pi H_(d-1) |
| :-- | --: | --: | :-- | :-- | :-- | :-- | --: |
| xxz | -0.9 | 2 | 4, 6, 8 | 0.279660, 0.321753, 0.350455 | 0.356074, 0.409669, 0.446213 | 0.458718 | 3.141593 |
| xxz | 0 | 2 | 4, 6, 8 | 1.722410, 1.635560, 1.606735 | 2.193040, 2.082460, 2.045759 | 1.994650 | 3.141593 |
| xxz | 0.5 | 2 | 4, 6, 8 | 2.282296, 2.145546, 2.099606 | 2.905909, 2.731795, 2.673302 | 2.593278 | 3.141593 |
| xxz | 0.9 | 2 | 4, 6, 8 | 2.695999, 2.528522, 2.470850 | 3.432652, 3.219414, 3.145983 | 3.049106 | 3.141593 |
| xxz | 1 | 2 | 4, 6, 8 | 2.791253, 2.617924, 2.558247 | 3.553933, 3.333245, 3.257262 | 3.156991 | 3.141593 |
| xxz | 1.05 | 2 | 4, 6, 8 | 2.767975, 2.596033, 2.536809 | 3.524296, 3.305371, 3.229965 | 3.130515 | 3.141593 |
| ashkin_teller | -0.5 | 4 | 2, 3, 4 | 0.371003, 0.362732, 0.359539 | 0.472376, 0.461845, 0.457779 | 0.453295 | 5.759587 |
| ashkin_teller | 0 | 4 | 2, 3, 4 | 0.608964, 0.578258, 0.568067 | 0.775357, 0.736261, 0.723285 | 0.705215 | 5.759587 |
| ashkin_teller | 0.5 | 4 | 2, 3, 4 | 0.794593, 0.746983, 0.730989 | 1.011708, 0.951089, 0.930724 | 0.902863 | 5.759587 |
| ashkin_teller | 1 | 4 | 2, 3, 4 | 0.930418, 0.872641, 0.852749 | 1.184644, 1.111082, 1.085754 | 1.052330 | 5.759587 |
| spin_orbital | 0 | 4 | 4 | 2.791253 | 3.553933 | unresolved | 5.759587 |
| spin_orbital | 1 | 4 | 4 | 1.283119 | 1.633718 | unresolved | 5.759587 |
| spin_orbital | 2 | 4 | 4 | 3.402852 | 4.332646 | unresolved | 5.759587 |
| spin_orbital | 4 | 4 | 4 | 4.315961 | 5.495252 | unresolved | 5.759587 |
| spin1_blbq | -0.785398 | 3 | 4 | 2.884678 | 3.672887 | 3.176231 | 4.712389 |
| spin1_blbq | 0.785398 | 3 | 3 | 4.183346 | 5.326402 | unresolved | 4.712389 |
| fz_parafermion | 2 | 2 | 4, 6, 8 | 0.568067, 0.560954, 0.558496 | 0.723285, 0.714229, 0.711099 | 0.706998 | 3.141593 |
| fz_parafermion | 3 | 3 | 3, 4 | 0.723011, 0.708556 | 0.920567, 0.902162 | 0.878520 | 4.712389 |
| fz_parafermion | 4 | 4 | 2, 3, 4 | 0.858095, 0.805112, 0.786991 | 1.092561, 1.025101, 1.002029 | 0.971293 | 5.759587 |
| fz_parafermion | 5 | 5 | 2, 3 | 0.911315, 0.850511 | 1.160322, 1.082904 | 1.020955 | 6.544985 |
| fz_parafermion | 6 | 6 | 2, 3 | 0.947833, 0.880407 | 1.206818, 1.120968 | 1.052290 | 7.173303 |

[numerical] Largest eigenvector residual: 4.062e-10. 98 ring ground states. Maximum Hilbert dimension 65536.

[assumed, infrared interpretation] The self-dual ferromagnetic FZ and Ashkin-Teller rings use the untwisted vacuum, XXZ and the spin-1 TB point use even rings, ULS uses multiples of three, and spin-orbital comparisons use multiples of four. These are the known vacuum-compatible sequences at the named critical points. Away from known critical points Q is still a finite-ring observable, and is not labeled an infrared central-charge measurement.

[numerical, limitation] The fit at XXZ Delta=1.05 remains close to critical on these rings although its asymptotic phase is gapped. The Heisenberg finite-size estimator exceeds pi slightly. These are controls against identifying small-ring drift with a counterexample. The q=5,6 results only reach L=6 and do not establish asymptotic A. No DMRG was run.
