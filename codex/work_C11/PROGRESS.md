# C11 progress

Closed-book calculation begun 2026-10-27. Only files created in this session will be read. The supplied context is treated as input, with any additional modeling assumptions recorded explicitly.

The first issue is whether the one-port Hamiltonian plus the stated mirror compression fixes a two-port cell Hamiltonian and its excited-state resolvent. First-order port compressions do not in general determine the second-order interaction. I will separate conclusions requiring only the given data from explicit one-port calculations and any proposed two-port completion.

The exact one-port resolvent has two excited gaps, D = sqrt(Delta^2 + 4 r g^2) and d = (Delta + D)/2. Its spectral port tensors were derived analytically and checked by dense diagonalization for q=2 and q=3. For q=2 the isolated-bond four-Majorana coefficient is positive at every nonzero finite g/Delta and starts at (83/2) lambda^2 g^4 / Delta^5. This statement is conditional on each endpoint having the supplied one-port resolvent. It is not yet the full two-port-chain answer.

An explicit additional model, Hcell = Hleft + Hright on two independent port-partner pairs sharing the logical register, has been diagonalized. At g/Delta=0.1 its ground degeneracy and Clifford port compressions hold numerically for r=2,4, but its compression coefficient differs from the supplied single-port alpha. Its second-order bond coefficients also differ. Thus importing the single-port resolvent into this completion is unjustified. The q=3, r=4 bond correction is anisotropic, with nonzero mixing between distinct flavor pairs. Orthogonality of port labels does not establish O(r) symmetry of higher-order processes.

The leading one-port bond tensor has now been obtained for arbitrary normalized odd F_a. In a basis Q_A=gamma_c gamma_d, A=(c<d), define J_ab,A=delta_ac delta_bd-delta_ad delta_bc and L_ab,A=tr(F_c F_a F_b F_d-F_d F_a F_b F_c)/q. The coefficient of Q_A tensor Q_B is lambda^2 g^4/Delta^5 times [119 J^T J/4 - 7(J^T L+L^T J)/4 - L^T L/2]. No complex conjugation is taken in these contractions. q=2 gives 83/2. q=3 gives an anisotropic 6 by 6 tensor, recorded by leading_tensor.py.

For q=3 there is also an earlier one-cell charge term. Each incident bond contributes (3/2) lambda^2 g^2/Delta^3 times i(gamma_1 gamma_2+gamma_3 gamma_4). This is a dimension-one chiral-current perturbation; translation forbids R-L masses but not uniform flavor currents.

The natural two-port completion was evaluated including adjacent-bond SW processes. Their quartic contribution starts at g^6 in the tested r=2,4 cases, while the same-bond g^4 tensor agrees with the one-port result. At r=4 the commensurate continuum theory separates into an isotropic SU(2) sector with coefficient 712 and an axial SU(2) sector with longitudinal/transverse coefficients 232 and 1024, in units a lambda^2 g^4/Delta^5. The axial sector is on the marginally relevant BKT side despite all three coefficients being positive. However, the generated charge term detunes its umklapp by a polynomially small wavevector, parametrically larger than the exponentially small commensurate gap at weak lambda. The distinction between the commensurate and detuned theories is essential.

The tensor calculation extends analytically to q_+=1, q_-=m, r=2m. Writing q=m+1, the continuum current matrix has three invariant eigenvalues, in units a lambda^2 g^4/Delta^5:

- U(1) charge: 952 - 8 q^2(q+7).
- SU(m) currents: 8(119 - 7q - q^2).
- Pair currents (commensurate umklapp): 8(119 + q^2).

These follow from L_ab,A = tr([F_d,F_c] F_a F_b)/q, J^T L = -tr(C_A C_B)/q, and L^T L = -2 t_A t_B, where t_A is the trace of C_A over the plus block. Numerical evaluations for m=1 through 6 confirm the resulting full current spectra. The SU(m) coefficient becomes negative already at m=7, q=8, r=14. Detuning the common filling does not remove this flavor backscattering. Thus the small-g, flavor-symmetric family does not supply an automatically stable unbounded-c route.

A different conditional stabilization is available. Preserve exact U(1)^m flavor charges and add distinct small chemical potentials for the m complex fermions. Distinct fillings make flavor-exchange backscattering oscillatory; a common nonzero shift removes two-particle umklapp. The remaining weak interactions are forward scattering and give a multicomponent Luttinger liquid. This requires an explicit microscopic implementation of the charges and the cell-local chemical-potential terms. It is an open region of couplings, not tuning a marginal coupling to zero. It changes the first-order target from coincident Majorana cones to m Dirac bands.

Numerical validation of the full second-order SW formula, including adjacent bonds, has passed for the natural r=2 two-port completion. Exact sparse spectra of two and three cells at lambda=0.004, 0.002, 0.001 differ from the SW predictions by O(lambda^3), with error/lambda^3 approaching 0.1716 and 0.3179 respectively. See validation_output.txt.


RESULT.md is complete. Final algebra checks extended the U(m) current-spectrum calculation through m=8, confirming the sign change at m=7 and the analytic multiplicities. All created Python sources compile. The result distinguishes exact finite-cell algebra and SW coefficients, finite-chain numerical validation, and conditional weak-coupling continuum conclusions. No claim of unconditional stability of the supplied unspecified mirror implementation is made. The proposed distinct-filling repair requires exact U(1)^m symmetry in the full history-cell Hamiltonian. No existing files, external sources, or project memory files were read.
