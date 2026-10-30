# C122 result

2026-10-29.

[Open] No general bound log[mu(Y)/mu(X)] <= K'(d)+k M f(d) beta/L was obtained. The operator-order version of the proposed locality route is ruled out already on XX. This does not disprove the scalar target or the bound on A.

[Proved] The correct comparison is

    c_ord = min{c : Y <= c X}
          = ||exp(beta H1/2) exp(-beta H0/2)||^2,
    log[mu(Y)/mu(X)] <= log c_ord.

But log c_ord >= beta Tr(sigma V). For XX on every even L>=4,

    log c_ord >= (2 beta/L)csc(pi/L)-2 log 2.

At beta=tau L this is (2/pi)tau L+O(1), while the requested ceiling is O(tau). Thus even an optimal estimate of this full operator norm loses a factor L. Taking L to infinity first does not repair it.

[Proved] An explicit nearest-neighbor counting argument gives a volume-independent half-time expansional constant

    G(b)=exp[6b+2(exp(2b)-1)],    b=beta M.

This is an upper bound, not a claim about optimal growth. Combining it with the settled trace bound gives

    |log[mu(Y)/mu(X)]|
      <= min{4 beta M+2w log d, 12 beta M+4(exp(2 beta M)-1)}.

The general partition ratio still has only the previously proved |R_L|<=4 beta M. No new universal coefficient for A follows.

[Proved] The spatial tail from the same counting has factorial decay at fixed beta, but the width it guarantees is O(exp(2 beta M)) for a fixed approximation error. A local approximation in the mu comparison also pays a factor q in its norm error. Enlarging the seam therefore does not balance the entropy cost into a finite bound at beta=tau L. Even an assumed cost proportional to beta M/w would balance against 2w log d at order sqrt(beta M log d), still divergent in that limit.

[Proved] Projecting to low energies leaves the exact vacuum normalization

    log[mu(Y)/mu(X)]
      = beta(E_(2L)-2E_L)+log[mu(P0)/mu(P1)]+eta_0-eta_1.

A gap delta~v/L alone gives eta_i<=log[1+d^(2L)exp(-beta delta_i)/g_i], which is not uniform at beta=tau L. Even perfect ground-state projection leaves the unknown Casimir slope. The remaining lemma must compare the optimized SDP witnesses, allowing them to change between the two Gibbs operators, with a cost K'(d)+k M f(d) beta/L. PROOF Section 7 states it explicitly. A full positive-operator comparison is too strong.

[Numerical] Solved 70 full finite-temperature SDPs and 14 ground-projector SDPs using NumPy and SciPy, including XX and Heisenberg on 8 total sites and SU(3) on 6 total sites. The minimal seam contains all four reconnection sites. Additional checks use a six-site seam and the length-two calibration. All physical sectors are retained. Numerical primal and dual bounds give log-ratio intervals of width at most 3.11e-6. These are full SDPs solved to numerical accuracy, not symbolic exact values or directed-rounding certificates.

[Numerical] At beta/L=16 with the four-site seam,

| Model | L | log[mu(Y)/mu(X)] | beta(E_(2L)-2E_L) | Difference |
| :--- | ---: | ---: | ---: | ---: |
| XX | 4 | 13.89343 | 13.77928 | 0.11415 |
| Heisenberg | 4 | 22.34728 | 22.33002 | 0.01726 |
| SU(3) | 3 | 66.98232 | 66.93354 | 0.04879 |

[Numerical] The differences agree with the respective ground-state offsets. Widening the qubit seam to six sites changes those offsets to 0.08386 and 0.04043, leaving the slope unchanged. The finite sizes do not establish an interacting large-L estimate.

[Proved, inherited] For XX, C118 already implies log[mu(Y)/mu(X)] <= 2(w+1)log 2+(5pi/3)beta/L. For diagonal classical bonds it implies the ceiling (2w+1)log d. These remain the available target-scale subclass bounds.

[Work record] Full derivations and checked primary references are in [PROOF.md](PROOF.md). Numerical details are in [CHECKS.md](CHECKS.md), with reproducible implementations in checks.py, sdp.py and analytic_checks.py. No large interacting calculation, cluster submission, independent agent review, or write outside work_C122 was performed. The local [STATE.md](STATE.md) carries the handoff.
