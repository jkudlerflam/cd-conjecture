# C94 result

[Proved] C92 equation (29), with its stated uniformity in finite N, T, s and 0<=u<=J, is false. At a free vacancy-number crossing, the frozen-material subtraction leaves

    D_sigma(0)-D_good(0)
      = (a_v-a_(v+1))^2/(4T)
          [tau_sigma(Hmat_sigma^2)-tau_good(Hmat_good^2)]
        + O_N(1).

Here tau is normalized trace on the full legal internal space. This result retains the quantum sources and the auxiliary states.

[Proved, full C88] Take N=15, selected lengths (15), wrong lengths (3,3,3,3,3), W=15 and t=1. At s=0.610116148316937... the free vacancy sectors 15 and 16 cross. Rational interval arithmetic certifies that these are the only free ground sectors. The coefficient of 1/T in the response difference is greater than 534.73. The full scheduled auxiliary and port terms are included in the analytic bounds. Continuity in u also rules out a bound uniform over 0<u<=J.

[Proved limitation] This does not refute C92 equation (26), which allows a global subextensive term and takes the thermodynamic limit before T->0. At the finite-volume crossing, even the integrated residual can have an unfavorable term of order J. The global finite-volume allowance cannot be dropped. The thermodynamic integrated comparison remains open.

[Proved, finite volume] A Schur-complement bound gives R_sigma-R_good >= -2J^2 N^2 M^2/Delta when 4JNM<=Delta, where Delta is the fixed-vacancy free bath gap. It applies to arbitrary quantum material sources and needs no internal excitation gap. Its volume-dependent constant does not prove the thermodynamic estimate.

[Proved, conditional second-order theorem] In the zero-port vacuum-clock model, write the canonical contact resolvent kernel as chi_N(r)=k_N(r)-c_N/N and set K_N=sum_r(1+dist(r,0))|k_N(r)|. Then

    |B_sigma-B_good| <= 4 M^2 b W (K_N+|c_N|).

This matches translated selected runs and retains the full quantum source Gram matrix. It is uniform over wrong-cell lengths at a fixed bath kernel. A uniform bound on t(K_N+|c_N|), and an all-order relative remainder bound, are still unproved.

[Numerical] All 50 canonical logical-clock defect costs are positive. Ten thermal comparisons subtract the full frozen material free energy, with controlled vacancy tails. Both signs of the residual occur. In an unfavorable b=15 canonical comparison, the observed epsilon decreases from 3.16e-5 to 8.22e-7 as J/t decreases from 0.04 to 0.0025. In an unfavorable small-clock thermal comparison at T=J, it decreases from 4.69e-4 to 9.47e-6. These are finite-system results, not a uniform thermodynamic bound. Separate half-density resolvent checks reach N=9.

[Files] PROOF.md contains the derivations and their scope. CHECKS.md summarizes the computations. crossing_certificate.py reproduces the rational certificate. checks.py and coefficient_scaling.py reproduce the numerical data, using local copies of C89's quantum-clock implementation. PROGRESS.md and STATE.md record the handoff.

[Scope] No finite-K theorem, full mobile selection theorem or central-charge conclusion is claimed. All writes and Python execution stayed in work_C94. No network commands, external edits or delegated agents were used.
