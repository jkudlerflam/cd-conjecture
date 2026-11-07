# C125 result

2026-11-06. Full derivations are in [PROOF.md](PROOF.md). Reproducible numerical scope and failures are in [CHECKS.md](CHECKS.md).

[Proved] The proposed crossed-channel inputs do not imply the Casimir bound or eventual convexity of Theta in L^2. There are two separate obstructions. The global crossed spectrum condition and finite real momentum tilt fail already for XX. Even if both are granted as extra assumptions, the finite-dimensional positive spectral family below satisfies them and the trivial bounds while its infrared coefficient A is unbounded at fixed labels d and M.

## Corrected lattice reduction

[Proved] A positive T with finite traces gives

    R_L=S_2(T^L/Tr T^L),
    dR_L/dL=2(<epsilon>_(2L)-<epsilon>_L)<=0.

Thus D(L)=E_(2L)-2E_L is nonnegative and nonincreasing on the allowed length subsequence. Neither L D(L) monotonicity nor Theta convexity follows. Positivity supplies E''<=0 when the real-length beta limit exists, which is weaker than L^2 E''+L E'-E<=0.

[Proved] The transfer construction must use a compatible reflection and may require two spatial sites per step. For a checkerboard Trotterization, positive realigned gates give

    Z_m(L,beta)=Tr U_m^(L/2), U_m>=0.

Taking T_m=sqrt(U_m) gives the claimed powers on even rings. An eigenvalue compactness argument constructs an exact common positive spectrum reproducing all even physical traces L>=4. At every finite Trotter number the construction also covers L=2; the simple compactness argument does not establish equality of its limiting trace at that endpoint. XX has an explicit exact spectrum covering L=2 as well. These statements do not identify the conventional one-site QTM with a positive Hermitian matrix.

| Model, zero field | Established transfer statement |
| :--- | :--- |
| XX | [Proved] Positive crossed spectrum for even rings, explicitly known. |
| XXZ, h=SxSx+SySy+Delta SzSz, Delta>=0 | [Proved] Positive even-ring construction, including Delta>1. |
| XXZ, Delta<0 | [Proved] No on-site DLS reflection certificate. [Numerical] Negative or complex small-Trotter spectra at Delta=-0.9,-0.5. At Delta=-0.9 exact ring traces also violate R monotonicity. |
| Antiferromagnetic bilinear Heisenberg, any spin | [Proved] Positive even-ring construction. |
| SU(2) swap | [Proved] Covered by the Heisenberg construction. |
| SU(N) fundamental swap, N>2 | [Proved] The on-site reflection criterion fails. [Numerical] SU(3) has D(4)<0, excluding a positive spectrum reproducing all even ring traces. |
| Integrable SU(2) Sym^k, every k | [Proved] Positive even-ring construction. Its non-scalar reflected coupling eigenvalues are (k+1)/[2l(l+1)], l=1,...,k. |
| Integrable SU(N) Sym^k, N>2 | [Proved] The same-representation on-site reflection criterion fails, including the supplied Sym^2 families. No positive representation on their selected vacuum subsequences is established. |

[Proved] Odd spin rings cannot generally be included. For example XX has D(3)=-1. SU(N) vacuum compatibility, such as L a multiple of N in the supplied data, does not supply reflection positivity. Alternating a representation with its conjugate changes the model.

[Proved] Real powers require a summability domain. For the exact XX crossed spectrum Tr T^L is finite for real L>1 and divergent for L<=1. Existence of the real-length zero-temperature limit is another hypothesis in general.

## The Matsubara condition is false globally

[Proved] At fixed beta, XX has crossed states reached by a local Sz insertion with

    p_m=2pi m/beta,
    epsilon_m-epsilon_0
       =asinh(pi/beta)+asinh((2m-1)pi/beta).

The right side grows only logarithmically with m. No finite V makes epsilon-epsilon_0>=|p|/V valid for all such states. The exact crossed fermion product is given in PROOF section 3.

[Proved] Consequently the XX trace Tr[T^L exp(-a P_tilde)] diverges at every real a!=0. Commuting unitary circle rotations do not imply a finite real tilt or its analyticity.

[Proved from C109] The KMS commutator Laplace formula and V=4eM give, for 0<omega<V,

    |G_AB(i omega,x)|
      <=2||A||||B|| V/[omega(V-omega)] exp(-omega x/V).

Optimizing the exponential Lieb-Robinson bound gives the asymptotic rate omega/(4eM) below 4eM and log(omega/(4M)) above 4eM, with a polynomial prefactor at the matching boundary. A positive reflected spectral measure converts this into a bound only for detected transfer sectors. This is a valid low-frequency replacement for the false global claim.

## A countermodel even under the stronger assumptions

[Proved] Fix d=2, M=1 and v=1. For each n take n compact c=1 bosons in the crossed channel. At each beta retain complete conformal-energy shells up to total rank D_beta<=floor(exp(beta/2)). Keep their usual integer multiplicities, circle momenta and energy gaps. The resulting T_beta is a finite-dimensional positive operator, its real tilt is entire, and every state satisfies epsilon-epsilon_0>=|p|/v. Its full infrared torus limit is the n-boson CFT.

[Proved] Nevertheless

    0<=R_L<=log D_beta<=beta/2

for every L. This is stronger than the stated trivial bound and the inherited |R_L|<=4 beta M ceiling. A scalar vacuum shift, explicitly given in the proof, also makes the elementary dimension and bond-norm bounds on Z hold and gives Z_L(0)=2^L.

[Proved] Set B=pi n v/6 and b=M/2. The real-length ground energy is

    E(L)=-b+L b^2/(4B), L<=ell=2B/b,
         -B/L,           L>=ell.

For L<ell,

    d^2[-L E(L)]/d(L^2)^2=-b/(4L^3)<0.

For L>=ell,

    lim_(beta->infinity) L R_L/beta=pi n v/4,
    A=n v.

Increasing n makes both the ratio coefficient and the onset of constant c_A arbitrarily large while d,M and the entropy ceiling stay fixed. This is an ordinary positive spectral model, not a fixed-d local-chain counterexample. Its spectra at different beta are not asserted to come from a single microscopic Hamiltonian.

[Proved] For a species with v_s<V, the crossed null edge is |a|=L/v_s, outside the cone |a|<L/V guaranteed by a global V-bound. Convexity inside the common cone does not access all sector null limits. Even with equal velocities and an entire tilt, the countermodel still fails Theta convexity. The additional missing ingredient is a relation fixing the tilted free energy from its untilted slice, such as the continuum rotation identity, or a direct scale comparison of the finite-size anomaly.

[Proved conditional closure] One sufficient scalar addition is L D(L)<=L_0 D(L_0). Another is Theta convexity for L>=L_0. The latter gives

    A<=c_A(L_0)<=(12/pi)L_0^2 M.

These are extra hypotheses, not consequences of entropy monotonicity. Ground-energy convexity alone does not bound the finite-temperature remainder in R. Also multiplying the C116 estimate by L/beta leaves K L/beta, not K, at arbitrary aspect ratio.

## Numerical findings

[Numerical] XX full traces agree with ED to 2.28e-13. The explicit infinite-column spectrum agrees with free-fermion traces through L=1024, beta=32768, to 5.59e-9 in log Z. Free-fermion ring checks reach paired length 16384. All sampled C118 ratio bounds pass.

[Numerical] At XXZ Delta=-0.9 and beta=21.928721377,

    R_4=1.073477335, R_6=1.137473641.

Thus even unscaled R monotonicity fails outside the positive-transfer class. Its c_A also increases from 0.190911175 at L=4 to 0.243662250 at L=12, and Theta fails discrete convexity. [Proved] Perturbation about Delta=-1 pushes this failure beyond any fixed dimension-only length cutoff. This is not a violation of the proposed A bound.

[Numerical] Heisenberg and tested XXZ Delta=-0.5,0,0.5,2 pass discrete Theta convexity on lengths 4 through 12. SU(3) passes on vacuum-compatible lengths 3,6,9, despite E_8-2E_4=-1.406345439 on the nonvacuum even sequence.

[Numerical] All supplied Sym^2 Bethe families pass both c_A monotonicity and discrete Theta convexity. The JSON starts at 2N, so the exact E_N=0 point was added separately. For N=8,

    c_A(8)/M=7.125295634,
    A/M=6.597344573.

[Numerical] Momentum-resolved countermodel checks verify the cone, positive covariance Hessian and Renyi derivative in 81 samples. The SU(2) recoupling recurrence and harmonic-sum coefficients were checked with exact rational arithmetic through k=12; direct two-site matrix spectra were checked through k=8.

[Assumed or unresolved] No general theorem bounding A for physical fixed-d chains was obtained. A positive representation restricted to selected N>2 vacuum subsequences is not classified, and general reconstruction of continuous time rotations and local insertions from the limiting transfer spectrum is not supplied. These points are not needed for the negative implication result. No independent referee review was performed.

[Work record] All writes and computations stayed in work_C125. The project-level STATE.md was not modified because the task explicitly forbids outside writes; STATE.md here is the handoff record. run_typhon.slurm is ready for larger interacting checks but was not submitted. No Typhon connection was attempted.
