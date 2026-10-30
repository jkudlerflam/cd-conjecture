# C122 numerical checks

2026-10-29.

[Assumed normalizations] XX has h=Sx Sx+Sy Sy and M=1/2. Heisenberg has h=S dot S and M=1/2. SU(3) has h=P and M=1. These are the campaign conventions, not fitted values.

## Full SDP scope

[Numerical] checks.py solves the conditional-min-entropy optimization for both the joined and split Gibbs operators at beta/L=0,0.1,1,4,16. The main cases use C={0,L-1,L,2L-1}. They are XX L=4, Heisenberg L=4 and SU(3) L=3. Additional cases are XX/Heisenberg L=2 and XX/Heisenberg L=4 with C={0,1,3,4,5,7}. Thus there are 35 Gibbs pairs, 70 finite-temperature optimizations and 14 ground-projector optimizations. The main suite took 91.06 seconds with one BLAS thread.

[Numerical scope] The length-two rings retain the definition H_2=h_01+h_10. Their seam is the whole four-site Hilbert space, and mu reduces analytically to the largest eigenvalue. They are useful calibrations, not evidence of asymptotic scaling. SU(3) has one vacuum-compatible comparison, L=3 versus 6. No interacting large-L extrapolation is justified by these sizes.

[Proved reduction] Realness and conservation of each total color justify restricting Omega to real symmetric complement-charge blocks. This is symmetry averaging of the full SDP, not restriction to a selected thermal sector. The largest physical matrix has dimension 729. For the four-site seam the complement dimensions are 16 for the qubit L=4 cases and 9 for SU(3) L=3. The implementation retains all physical charge blocks in its positive-semidefinite constraint.

[Numerical] Each Gibbs operator is shifted by its own ground energy before exponentiation. The stored log_mu_Y and log_mu_X restore those scalar factors. checks.json contains these individual unnormalized values, their ratio, the direct partition ratio, both SDP bounds and residuals for every case. No low-energy truncation is made. Terms below floating-point range can underflow, as in any double-precision exponentiation, but no spectral threshold or rank cutoff is imposed by the code.

## Accuracy and independent checks

[Numerical] The maximum relative primal-dual gap among finite-temperature solves is 1.62e-6. The widest resulting interval for log[mu(Y)/mu(X)] is 3.11e-6. The smallest repaired primal slack eigenvalue is 1.10e-9. The largest eigenvalue of the rescaled dual partial trace is 0.999999999990001. These are floating-point checks of feasibility and weak duality, not interval-arithmetic certificates.

[Numerical] At beta=0 every bracket contains mu(I)=dim D. All 14 ground-projector brackets contain the independent pure-state formula (sum of Schmidt singular values)^2. The largest difference between a ground SDP midpoint and that formula is 3.60e-6. All C120 identities and entropy ceilings tested in the code pass.

[Numerical] analytic_checks.py adds two independent diagonal SDPs with the analytic optimum sum_D max_C W(C,D). It also checks the new expansional bound on twelve small Gibbs pairs, the explicit majorant recurrence and generating sum, and four width estimates including the q factor. All pass. The order norms are evaluated from Hamiltonian eigensystems only at conditioned inverse temperatures. checks.py likewise omits direct inverse-Gibbs order norms once beta times the larger spectral width exceeds 25. It does not replace tiny eigenvalues by a floor.

## Values of the normalization ratio

[Numerical] Four-site seam values, rounded consistently with the SDP accuracy, are

| Model | L | beta/L | R_L | log[mu(Y)/mu(X)] |
| :--- | ---: | ---: | ---: | ---: |
| XX | 4 | 0.1 | 0.00039 | 0.03398 |
| XX | 4 | 1 | 0.61939 | 0.94895 |
| XX | 4 | 4 | 3.37038 | 3.55897 |
| XX | 4 | 16 | 13.77927 | 13.89343 |
| Heisenberg | 4 | 0.1 | 0.00066 | 0.05236 |
| Heisenberg | 4 | 1 | 1.12241 | 1.40193 |
| Heisenberg | 4 | 4 | 5.58181 | 5.59977 |
| Heisenberg | 4 | 16 | 22.33002 | 22.34728 |
| SU(3) | 3 | 0.1 | -0.01153 | 0.45290 |
| SU(3) | 3 | 1 | 3.92328 | 4.23212 |
| SU(3) | 3 | 4 | 16.73338 | 16.78217 |
| SU(3) | 3 | 16 | 66.93354 | 66.98232 |

[Numerical] The ground-state calculations independently give

| Model and seam | mu(P0) | mu(P1) | log[mu(P0)/mu(P1)] |
| :--- | ---: | ---: | ---: |
| XX L=4, four sites | 8.49264069 | 7.57648311 | 0.11415 |
| Heisenberg L=4, four sites | 9.00000000 | 8.84598716 | 0.01726 |
| SU(3) L=3, four sites | 9.00000000 | 8.57146336 | 0.04879 |
| XX L=4, six sites | 4.00000000 | 3.67824450 | 0.08386 |
| Heisenberg L=4, six sites | 4.00000000 | 3.84151440 | 0.04043 |

[Numerical] The slopes E_(2L)-2E_L are 0.215301195 for XX L=4, 0.348906591 for Heisenberg L=4 and 1.394448725 for SU(3) L=3. Their scaled values L Delta_L/M are 1.72240956, 2.79125273 and 4.18334617. The length-two qubit calibrations give 2.34314575 and 4. These are finite-size quantities. They are not claimed to equal the infrared limits.

## The order comparison is much more expensive

[Numerical] At beta=4 and L=4, the direct full-order log constants are 10.30670627 for XX and 13.01627429 for Heisenberg. The corresponding optimized log mu ratios are 0.94895 and 1.40193. At beta/L=16, the rigorous lower bounds beta Tr(sigma V), evaluated from the Gibbs matrices, are 45.25483400, 64.00000000 and 128.00000000 for the three main cases. They already exceed the respective optimized ratios 13.89343, 22.34728 and 66.98232. The large-temperature-parameter order norms themselves are not evaluated through ill-conditioned inverse matrices.

[Proved interpretation] The full-order comparison controls positive test operators without the SDP partial-trace restriction and is therefore stronger than needed. PROOF Section 2 proves its failure at the target scale on XX for arbitrarily large even L. The numerical table is a small-size illustration of that analytic obstruction.

## Reproduction

[Work record] Run inside work_C122 using Python, NumPy and SciPy only.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 analytic_checks.py

[Work record] The outputs are checks.json, analytic_checks.json and the optional captured console file run_output.txt. sdp.py implements the optimization. ring_tools.py is copied from the inspected C120 source. No network command, large local interacting ED, Typhon connection or job submission was used. The conclusions distinguish a full finite-dimensional optimization from an unproved large-L scaling inference.
