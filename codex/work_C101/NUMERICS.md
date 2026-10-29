# C101 computations

All completed computations are small local numpy/scipy calculations. No Typhon job ran. SSH with `-o LogLevel=ERROR -o BatchMode=yes -o ConnectTimeout=10` failed with `Connection timed out during banner exchange`. The prepared `typhon_ed.slurm` runs the larger L=12 spin-orbital screen only inside SLURM. Its partition and resource availability have not been verified because the connection failed. The Python large mode refuses to run without SLURM_JOB_ID.

Reproduce the completed checks with one BLAS thread:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 bond_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 tj_certificate.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 hubbard_quadrature.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 adversarial_ed.py

`stiffness_check.py` is the small exploratory diagonalization that located the four-site t-J certificate. The exact certificate, not its floating eigenvalues, is used in the proof.

Bond checks independently construct two-site Hubbard Fock matrices, generalized t-J, graded permutations with fields, and XYZ matrices. They compare the complete spectra or widths against the formulas in PROOF.md. The maximum spectral discrepancy is below 8e-15.

The Hubbard quadrature solves finite Fredholm systems at 128 and 192 Gauss-Legendre points, increasing to 288 when the change exceeds 2e-6. Its 66 points have U/t=0.5,1,2,4,8,16 and Q/pi=0.1,...,0.9,0.95,0.99. It records density, chemical potential, both velocities, exact M at that chemical potential, and A/(2pi M). Four half-filled spin velocities are compared against the independent Bessel formula. This is a finite numerical scan, not the all-density proof.

The spin-orbital search uses a rung as the d=4 site and

    h=a P_s+b P_t+c P_s P_t.

Both P_s and P_t swap spin-1/2 factors across adjacent rungs. The exact optimized half-spread is half the range of a s+b t+c s t, s,t=+/-1. This covers products, unequal leg speeds, interacting ladders, SU(4), and competing signs in a single one-site nearest-neighbor family.

The completed search has 14 coupling points at L=4,6,8. Fixed zero magnetization on each leg gives dimensions 36,400,4900. This symmetry sector is known to contain the ground state for the two reference models. For unfamiliar couplings it is a restricted-sector screen, not a certified search over the full ground state. The maximum residual is 1.15e-11. No sampled two-size Q exceeds one. The largest is 0.8745965622 at SU(4), where the exact Q is 0.75. At two decoupled Heisenberg legs the fit is 0.5656260606, while the exact Q is 0.5. Thus the finite-size bias is substantial even on calibrated conformal examples.

The L=4,8 sequence is used for the displayed two-point estimates. L=6 is retained in the raw output but not mixed into the SU(4) vacuum sequence. A fitted A at other couplings is not an established CFT coefficient. Neither a small same-sector gap nor a stable-looking two-size estimate proves conformality. The saved gaps provide a diagnostic only.
