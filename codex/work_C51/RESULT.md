# C51 result

[proved conditional on a complete Gaussian scaling limit] For N=n+1 compact bosons with any positive quadratic Hamiltonian, including density-current mixing and arbitrarily separated velocities,

    c_R=c_L=n+1,        c_cd=(c_R+c_L)/2=n+1.

This is the sum of the Virasoro central charges in the convention assigning c=1 to a nonchiral compact boson. It does not require one common Lorentz velocity. [PROOF.md](PROOF.md) Sections 2 and 3 give the canonical diagonalization and current-algebra proof. The central-term calculation works in nonzero charge sectors, so irrational-density rounding does not alter c.

[assumed, precise boundary] A complete Gaussian scaling limit means normalized chiral current algebras, a positive quadratic Hamiltonian, the compact charged sectors and exhaustion of the infrared theory by these fields. Power-law exponents and a vanishing beta function for selected correlations do not establish all of this. C50's PROOF.md appeared during this run and was read. It explicitly leaves its all-scale RG comparison and source theorem open. It does not yet prove the hypotheses needed for a microscopic c=n+1 conclusion.

## Thermal statement and several velocities

[proved for the Gaussian limit; proved for the lattice under scalar estimate P] In physical length units the safe formula is

    f(T)=e_0-(pi/12) sum_(a=1)^(n+1)
               [1/u_(a,R)+1/u_(a,L)] T^2+o(T^2).

The requested -(pi/6) sum_a 1/v_a coefficient follows if u_(a,R)=u_(a,L)=v_a. The complex two-bus construction does not impose that equality. Velocities must be the final dressed speeds, not bare band slopes or eigenvalues of a noncanonical stiffness matrix. The thermal coefficient alone is weighted by inverse speeds and does not determine an unweighted c.

[proved, new sufficient scalar criterion] Section 7 gives a usable quantitative implication. Let q_h=E_*2^(-h). After subtracting the dressed Gaussian model and both ground energies, suppose the scalar RG contributions obey

    |delta_h(T)| <= C q_h^2 eta_h min{1,(T/q_h)^p},
    p>2,        eta_h bounded and tending to zero,

with a UV remainder O(T^p) and an exact summed representation. Then their sum is o(T^2). If eta_h=O((q_h/E_*)^omega), the error is

    O(T^2 (T/E_*)^min(omega,p-2)),

with one extra log(E_*/T) at p=2+omega. The proof sums the high- and low-energy scales separately. It also covers logarithmically slow approach to the fixed point.

[assumed, remaining thermal input] Neither C50's finite-scale scalar norm nor its conditional beta-function estimate establishes this temperature-sensitive bound. The scalar normalization and physical ensemble transfer must be controlled through the full RG. This estimate is an explicit sufficient extension to prove, not a result claimed for the present microscopic action. C46 already established the diagonal Gaussian coefficient and common-bond reference quotient pressure; those are cited as prior results.

## Irrational fillings and finite rings

[proved in the Gaussian theory] Write the chiral charge lattice as x=B(m+alpha_L), with m integral and alpha_L a bounded, possibly moving offset. Then

    E=e_0 L+(2pi/L){(1/2)x^T U x
                    +sum_(i,r)u_i r n_(i,r)-(1/24)sum_i u_i}.

At fixed charge and chiral levels k_i, the oscillator multiplicity is product_i p(k_i). These are current Fock towers carrying Virasoro representations. They must not be confused with one irreducible Virasoro module in each charge sector; the zero-charge Fock space already has a level-one current state absent from the Virasoro vacuum descendants.

[proved] The ground-state zero-mode energy is O(1/L), the same size as the Casimir term. It can oscillate with L at irrational density. A raw Casimir fit need not converge. Subtract the explicit minimizing quadratic charge energy, or retain the moving sectors and take convergent subsequences of alpha_L. Uniform shifted-lattice Gaussian bounds show that these sectors contribute only O(log L) to log Z at fixed temperature. They therefore do not change the thermodynamic coefficient. Thermodynamics takes L->infinity first, then T->0, at each fixed n and fixed bond.

[assumed, remaining spectral input] Passing these towers to the microscopic ring requires o(1/L) spectral matching with exhaustion and the actual physical gluing rules. A fixed local interaction error or an O(1) seam bound does not suffice. The interacting wrapper's complete compactification lattice and tower multiplicities were not established here. Gaussian charge labels may be emergent and need not be exact microscopic conserved charges.

## Scope and checks

[proved, definition delivered] The continuum theorem delivers the Virasoro central charge of the complete multichannel Gaussian infrared theory. The thermal theorem delivers its operational low-temperature signature under scalar estimate P. No lattice interval-entanglement asymptotic is claimed. If only the Gaussian current subtheory is identified, it gives c=n+1 for that subtheory; with the usual unitary CFT stress assumptions it gives the lower bound c_full>=n+1, not equality without completeness.

[numerical] `python3 -B checks.py` passes small controls of canonical diagonalization, nonzero-charge Virasoro commutators, oscillator multiplicities, moving-sector subtraction, theta bounds and scalar summation rates. Maximum canonical and Virasoro residuals are 3.56e-14 and 6.22e-15. These are algebra and asymptotic checks, not measurements on the 4464-state chain.

[proved, verdict] The implication from a complete compact Luttinger limit to c_cd=n+1 is proved with the velocity and zero-mode qualifications above. The stronger claim that the presently available constructive RG results already give the microscopic thermal coefficient or full conformal towers remains unproved. The fixed-d counterexample remains conditional. All files and the local [STATE.md](STATE.md) handoff are inside work_C51; no network or delegated agents were used.
