# C127 progress

2026-11-06.

[Assumed] The task rules restrict every write and Python execution to this folder. The project STATE.md was read; a local STATE.md will carry the handoff because outside writes are forbidden. Read the shared target, full background and reconnection note, REPORT 9.76 through 9.88, and the C125 proof, result, checks, progress, handoff and scripts.

[Proved] A fixed-frequency insertion map has rank at most d^2-1. Its range need not be invariant under crossed evolution. The spectral cyclic span can have arbitrarily larger dimension, so that rank is not a count of excitation species or of crossed eigenstates.

[Assumed or unresolved] The universal scalar boson-product bound has neither been proved nor refuted. Failure of a momentum-preserving spectral comparison would not itself refute the scalar trace bound with an adjustable constant.

[Heuristic, under audit] Small checkerboard transfer matrices will be diagonalized in discrete time-rotation sectors. XX will also be tested against the exact infinite-column product, keeping finite-Trotter error separate from numerical diagonalization error.

[Proved] XX supplies a lattice rank obstruction. At fixed external Matsubara frequency, the reflected density two-point function has infinitely many positive spectral masses. The fixed-frequency insertion map still has rank at most three. This refutes the identification of insertion rank with the number of detected energies, not the scalar boson-product bound.

[Proved] Using convex midpoint quadrature in the exact XX crossed product gives log Ztilde(L) <= log 3 + beta[csc(pi/L)-L/pi] for every even L>=2. This implies an explicit two-species bosonic scalar majorant with V=4eM and K=log 3-2 log(1-e^(-1)).

[Proved] The diagonal ferromagnetic d-state chain forces K(d)>=log d in any oscillator-only bound of the proposed type. A Gaussian charge factor, if its spectral comparison and stiffness identification are assumed, obeys log Q <= 4(beta/L) sum rho <=16 M beta/L even for arbitrarily many charge directions.

[Numerical] The first local run covers 48 momentum-resolved positive transfer matrices, m=2,3,4, beta=0.25,1,2,8, Delta=0,0.5,1,2. Scalar S=6,V=4eM,K=log 3 passes this sample. Momentum-sector comparisons fail. Coarse m=2 at beta=8 violates the continuum LR cone, which is not a finite-step theorem. The exact XX product cutoff at L=2,beta=64 leaves a 0.00519 log-trace tail, bounded explicitly rather than treated as roundoff.

[Assumed or unresolved] Allowing arbitrary dimension-dependent S and K makes existence of the scalar bosonic bound equivalent, up to explicit constants, to the original scalar trace target. No general proof or lattice counterexample to that target has been found. The final proof will state this limitation explicitly.

[Proved] Completed PROOF.md and RESULT.md. Besides the rank obstruction, the notes prove the XX scalar bound, K(d)>=log d from a diagonal lattice chain, the conditional charge bound, equivalence of the unrestricted scalar product formulation to the scalar target, and a conditional reference bound using logarithmic ultraviolet dispersion. Independent XX stacks inherit the proved scalar comparison.

[Numerical] Replaced cancellation-prone Fourier sector extraction by positive cyclic convolution. Independent audit checks the exact XX resolvent and reflected bubble, including a frequency-tail bound; the XX integral and quadrature inequalities; both elementary boson-product bounds; the Gaussian charge bound; the square-root Fock count; both Gibbs comparison directions; and the nonlinear reference product with its tail. All stated assertions pass. The momentum-sector failure and the coarse-step cone failure remain explicit expected failures.

[Work record] CHECKS.md and local STATE.md record numerical scope, errors, and the unresolved interacting question. Larger m=5,6 work is SLURM-ready but unrun. No outside files, network shell commands, Typhon connection, or sub-agents were used. The general scalar decision requested in tasks 2 and 4 remains unresolved; no claim of a completed C(d) proof is made.
