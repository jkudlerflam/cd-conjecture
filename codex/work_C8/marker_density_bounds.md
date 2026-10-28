# Legal-sector marker selection at fixed nonzero coupling

Status: proved energy and density bounds for the explicit C8 atomic
formulation, using its inherited exact first-order projection. No
conformal limit, effective action or illegal-sector comparison is
claimed. In particular, "ground state" below means the ground state
among legal sectors with at least one marker.

## 1. Ingredients and normalization

Subtract the common base energy e_base L. The legal microscopic
Hamiltonian in cell coordinates is

    H=sum_j A_j+epsilon V0, epsilon>0,
    J=sup_j sum_(X contains j)||V0_X||.

Use the same A, n=ell-b, Delta, gamma and
C=max(2/Delta,1/gamma) as canonical_length_bounds.md. In particular

    |n-P1| <= C A.                                        (1)

The inherited code compression is the sum of a length hopping chain
of amplitude J0>0, the specified translation-invariant payload chain,
and kappa per marker. Write e_pay for the infinite-chain payload
ground energy per cell, and B_pay for the norm of its intercell bond,
with all onsite fields assigned onsite. These are finite and explicit
for the existing free-fermion target. For example run9_audit.md gives

    e_pay=sum_f[-2 sin(k_f)/pi+mu_f(k_f/pi-1/2)]

at unit payload hopping. The implemented marker coefficient is

    kappa=2J0/pi-e_pay.                                   (2)

Define beta=J0+B_pay. This is a projected boundary bound, whereas J
is the full unprojected perturbation norm. They must not be conflated.
All constants remain fixed as L increases. The direct-sum on-cell
Hilbert space remains untruncated.

## 2. A chemical tilt without an unbounded perturbation norm

[proved] Set mu=epsilon alpha with |alpha|<=2J0. On one cell define

    A_mu=A-mu(n-P1).

The operators commute, and (1) gives

    A_mu >= (1-C|mu|)A,
    gap(A_mu) >= Delta(1-C|mu|).                          (3)

Its kernel is still P0+P1. For N=sum n_j,

    H-mu N=sum_j A_mu,j+epsilon V0-mu sum_j P1_j.         (4)

Only the last two terms in (4) are treated as bounded interactions.
This is why no norm of the unbounded length operator enters.

Apply the arbitrary-retained-theory energy bound of
defect_energy_bounds.md with W=-mu sum P1 and F=epsilon V0. Its
diagonal norm is at most |mu|, leakage bound is epsilon J, and

    g_mu=Delta(1-C|mu|)-2(epsilon J+|mu|).

Thus, if

    D0=2J0 C Delta+2J+4J0,
    epsilon D0 <= Delta/2,                               (5)

then g_mu>=Delta/2 and

    H-mu N >= E_ref(mu)-2M epsilon^2 J^2/Delta.           (6)

Here E_ref is minimized over the entire code reference space, not
over the given N. This deliberate grand-canonical comparison can
then be optimized in mu to recover a canonical density cost.

For the physical quotient, use the lifted space with x in Z/LZ and
all total cell lengths temporarily allowed. The same local transfer
rules define an extension conserving total length. Only the fixed-L
sector is quotiented by the reference action. Fourier transformation
in x turns the reference seam into a length-hopping flux. The proof
of (6) also allows this finite spectator x: filling missing cell
codes leaves x untouched, and all removed-interaction norm sums are
still counted on their incident cells. No Lieb-Robinson argument or
locality assumption on this spectator is used in the energy bound.
Restriction to the physical fixed-L subspace preserves (6).

## 3. Projected pressure and a canonical lower bound

The bulk length ground energy at length-bit density rho in [0,1] is

    e_len(rho)=-(2J0/pi)sin(pi rho),
    p_len(alpha)=min_(0<=rho<=1)[e_len(rho)-alpha rho].     (7)

[proved] Uniformly in the reference seam flux and in alpha,

    E_ref(epsilon alpha)
      >= epsilon M[kappa+e_pay+p_len(alpha)]-2epsilon beta.
                                                               (8)

One finite-volume boundary argument suffices. For an open M-site
chain with bond norm B and arbitrary onsite terms, its ground energy
E_open satisfies |E_open-M e_infinity|<=B. The lower inequality
follows by tiling the infinite line with independent open-chain
ground states and bounding every joining bond by B. For the upper
inequality, cut a long-chain ground state into M-site blocks, bound
the removed bonds by B, and use E_open as the lower energy of each
block. Take the length of that long chain to infinity. Restoring a
possibly twisted seam costs at most another B. Apply this separately
to length and payload. The density field alpha and kappa are onsite,
so they add no boundary cost. The existence and value of the bulk
energies in (7) and e_pay also follow directly from their free bands.

[proved] Put rho=(L-bM)/M. Optimizing (6)-(8) over |alpha|<=2J0 gives

    E_phys(M,L) >= epsilon M G(rho)-2epsilon beta
                                   -2M epsilon^2 J^2/Delta,   (9)

where the convex function is

    G(rho)=2J0/pi-2J0 rho,                         rho<=0,
           (2J0/pi)[1-sin(pi rho)],                0<=rho<=1,
           2J0/pi+2J0(rho-1),                      rho>=1.

Indeed the derivative e_len'(rho)=-2J0 cos(pi rho) ranges over the
whole permitted alpha interval. Its Legendre recovery equals e_len
inside [0,1] and its endpoint supporting lines outside. Adding (2)
proves (9), including sectors that have no all-code state. G is
nonnegative and has its unique minimum at rho=1/2. The inequality is
an energy comparison for the interacting microscopic chain, not just
the first-order projected spectrum.

For small M where a cyclic bond convention needs special treatment,
increase the harmless boundary constant in (8) by a fixed amount.
The density argument is thermodynamic. Alternatively the explicit
C8 conventions for M=1,2 can be treated as finitely many exceptional
sectors. Their wrong-length atomic cost grows linearly in L, while
their perturbation norm is bounded by M epsilon J, so they cannot
be the large-L minimum. All formulas below use M>=3.

## 4. A physical upper trial near half length filling

Let a=b+1/2, and choose M_*=round(L/a), N_*=L-bM_*.
For sufficiently large L, M_*>=3 and 0<=N_*<=M_*. The rounding error
obeys |N_*-M_*/2|<=a/2.

[proved] There is a physical legal trial with energy at most

    E_legal(L) <= epsilon B_L,
    B_L=2beta+pi J0 a^2/(4M_*).                           (10)

Use an open payload ground state, an open length-chain ground state
at N_*, and the pinned-marker embedding of canonical_length_bounds.
The payload open energy is at most M_* e_pay+B_pay. The length open
energy is at most M_* e_len(N_*/M_*)+J0. To see the latter, remove
one bond from the zero-flux hard-core ring, whose canonical energy is

    -2J0 sin(pi N_*/M_*)/sin(pi/M_*)
           <= M_* e_len(N_*/M_*).

Removing the bond costs at most J0. The physical seam restored in
the pinned-marker trial costs at most beta. Finally

    M_* G(N_*/M_*) <= pi J0 (N_*-M_*/2)^2/M_*
                        <= pi J0 a^2/(4M_*),

using 1-cos x<=x^2/2. All trial states are in the exact atomic kernel,
so this is also an upper bound on the unprojected microscopic energy.

## 5. The actual legal ground marker density is interior

[proved] Fix 0<delta<1/2 and define

    g_delta=(2J0/pi)[1-cos(pi delta)]>0,
    s_delta=2J0 sin(pi delta)>0.

In addition to (5), assume

    epsilon <= g_delta Delta/(4J^2).                      (11)

Then every ground sector among legal positive-marker sectors satisfies
|rho-1/2|<delta for all sufficiently large L. One explicit size
condition, together with the trial admissibility above, is

    g_delta L/[2(a+delta)] > 2beta+B_L.                   (12)

Proof. Convexity and reflection about 1/2 give, for |rho-1/2|>=delta,

    G(rho)>=g_delta+s_delta(|rho-1/2|-delta).

The second-order subtraction in (9) is at most g_delta M/2 by (11).
For rho>=1/2+delta, M<=L/(a+delta), and

    M[G(rho)-g_delta/2]
      >= s_delta L+[g_delta/2-s_delta(a+delta)]M
      >= g_delta L/[2(a+delta)].

The coefficient of M in the middle line is negative because
g_delta<=delta s_delta and a>0. For rho<=1/2-delta the analogous
linear expression is minimized at M=L/(a-delta), giving the larger
lower bound g_delta L/[2(a-delta)]. The same comparison applies to
all exterior densities, because the supporting line remains below
the convex extension G. Equation (9) now exceeds the physical trial
bound (10) whenever (12) holds. This excludes every exterior sector.
The exceptional fixed small M sectors were excluded after (9).

[proved] Once the density is in [0,1], (9)-(10) give the more useful
quantitative estimate

    |rho-1/2| <= sqrt{ pi/(4J0) *
                 [2epsilon J^2/Delta+(2beta+B_L)/M] }.    (13)

This uses 1-cos(pi x)>=2x^2 for |x|<=1/2. The displacement is
O(sqrt(epsilon)) in the thermodynamic limit. It is not shown to be
O(epsilon), unique, differentiable in the couplings, or governed by a
positive interacting compressibility.

## 6. Consequence for the length-tail theorem

Take delta=1/10 and eta=1/5 in canonical_length_bounds.md. Impose
(5), (11), and

    epsilon J <= min[eta/(4C), eta gamma/16].              (14)

All right sides are positive for the one fixed finite candidate.
For sufficiently large L, every legal ground marker sector then has
rho in (2eta,1-2eta) and M>=1/eta. Its ground state has the proved
factorial cell-length tail. Its within-sector low-energy window also
has the operator-norm tail and the spectral cutoff bound established
there. This removes the need to impose the reference density by hand.

It does not compare the invariant illegal selector sectors or the
marker-free sector to these states. It does not turn the effective
projected Gaussian theory into a CFT of the full microscopic chain.
The latter still needs a controlled low-frequency local elimination
and infrared identification at one fixed epsilon. Energy and length
control alone do not meet that obligation.

## 7. Verification

marker_density_checks.py checks the chemical-tilt energy bound against
literal physical scalar-ring diagonalization, tests every allowed
marker number at the stated lengths, and checks the pressure formula
and supporting-line inequalities independently. These controls have
no payload, J0=1 and kappa=2/pi. They test the mechanism and quotient,
not the high-channel C8 CFT. Inputs, finite sizes and residuals are
stored in marker_density_checks.json.
