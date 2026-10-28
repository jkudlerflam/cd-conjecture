# Exact spectrum of the two-length marker fluid

Status: proved finite-volume spectrum and asymptotic towers for the
projected marker Hamiltonian. Identifying the resulting free oscillator
towers with a compact boson uses the standard free-fermion continuum
input. Independent entropy and Casimir estimates are numerical. This
is NOT a finite-coupling theorem for the microscopic selector or its
payload. It does not settle the c-d conjecture.

## 1. Physical Hilbert space and the missing ring constraint

Fix integers b>=1, M>=2, 0<r<M and L=bM+r. A basis state is a set
of M marker positions on a ring of L physical sites, with every
consecutive spacing b or b+1 and exactly r long spacings. A marker
moves right with amplitude -J when its left spacing is b and its
right spacing is b+1, exchanging the two. The reverse move is its
adjoint. Here J>0. A chemical term kappa M is allowed when different
M sectors are combined at fixed L.

This is the exact first-order graph in selector_wrapper.md, with
J=u A_b^2 beta. Its physical Hilbert dimension is

    dim H_(L,M)=(L/M) binomial(M,r).                      (1)

To count, distinguish one reference marker at x in Z_L and record
the successive gap bits eta=(eta_0,...,eta_(M-1)), sum eta=r.
Each physical configuration has M choices of reference. Thus the
cover has L binomial(M,r) elements, with the free identification

    (x,eta) ~ (x+b+eta_0,S eta),                         (2)

where S shifts the gap sequence one step left. This also explains
why a periodic XX chain of M sites alone is not the full spectrum.

The physical hopping graph is connected. On the reference cover,
nonreference marker moves give all adjacent exchanges of 0 and 1
on an OPEN gap chain, hence connect all gap strings of weight r
at fixed x. Arrange eta_(M-1)=0, eta_0=1, which is possible since
0<r<M. Moving the reference marker right changes x to x+1. Open-chain
exchanges can then restore any desired gap string. Repeating reaches
every x, so the cover and its physical quotient are connected. No
numerical connectivity assumption is needed.

## 2. Complete free-fermion spectrum

Use the coefficient convention F(x,eta)=exp(iPx) phi(eta), with
P=2pi ell/L, ell=0,...,L-1. The quotient condition (2) gives

    phi(S eta)=exp[-iP(b+eta_0)] phi(eta).                 (3)

Define Q_eta=sum_(j=0)^(M-1) j eta_j and write

    phi(eta)=exp[-iP Q_eta/M] psi(eta).

Since Q_(S eta)=Q_eta-r+M eta_0, (3) becomes

    psi(S eta)=exp[-iP L/M] psi(eta).                    (4)

The same gauge distributes the reference-hop twist evenly among all
M gap bonds. Jordan-Wigner on these r hard-core particles gives
distinct single-particle wave numbers

    q_a=(2pi/M)[I_a+(r-1)/2],   I_a in {0,...,M-1},       (5)

with the conventional parity-dependent fermion boundary condition.
Use Slater coefficients exp(-i q_a j). Then cyclic shift of the
configuration multiplies psi by exp(i sum_a q_a), so (4) imposes

    m+ell=0 mod M,
    m=sum_a I_a+r(r-1)/2 mod M.                          (6)

With this sign convention the single-particle energy is
-2J cos(q_a+P/M). The COMPLETE spectrum, with multiplicity, is

    E=kappa M-2J sum_(a=1)^r cos(q_a+2pi ell/(LM)),        (7)

for every r-element subset of {0,...,M-1} and every ell satisfying
(6). There are no further choices or constraints. The Fourier
decomposition in x, gauge transformation, and Jordan-Wigner map are
unitary in their indicated sectors. They prove completeness, rather
than just producing candidate eigenvalues. They also reproduce (1).

At P=0 the positive ground vector of the periodic gap XX chain is
cyclically invariant and lifts through (2). Connectedness and
Perron-Frobenius make it the unique physical ground state at fixed
M. Filling the symmetric Fermi interval gives

    E0(L,M)=kappa M-2J sin(pi r/M)/sin(pi/M).              (8)

Every nonzero-momentum state is therefore above this state, even
though the twist in (7) varies with momentum.

## 3. Thermodynamic density and every ring length

Let nu=r/M in (0,1), rho=M/L=1/(b+nu). Up to a constant e_base L,
the bulk energy density is

    f(rho)=kappa rho-(2J/pi)rho sin[pi(1/rho-b)].          (9)

Its derivatives are

    f'(rho)=kappa-(2J/pi)[sin(pi nu)
                                   -pi(b+nu)cos(pi nu)],
    f''(rho)=2pi J (b+nu)^3 sin(pi nu)>0.                (10)

Consequently any prescribed interior nu is the unique bulk minimum
if

    kappa=(2J/pi)[sin(pi nu)-pi(b+nu)cos(pi nu)].          (11)

For kappa=0 this reduces to tan(pi nu)=pi(b+nu), the earlier
mobile-selector equation. To select nu=1/2, use kappa=2J/pi.
The microscopic selector chemical potential can supply this term
at first order by a small change in its marker cost.

At fixed L the allowed M are all integers in [L/(b+1),L/b]. For
large L, minimizing (8) rounds the interior density to M=rho L+O(1).
Strict convexity excludes other macroscopic densities. The finite-size
correction shifts the continuous minimizer by only O(1/L), so for
large L at most two adjacent integers compete for the ground state.
Each fixed-M ground state is unique. Thus ground degeneracy stays
bounded on EVERY sufficiently large ring, with the usual charge-sector
offset from rounding rho L. There is no required isolated defect and
no defect band with an L-dependent degeneracy.

## 4. Conformal towers and the physical Luttinger parameter

The gap-chain Fermi velocity is v_gap=2J sin(pi nu). Because its
length is M=rho L, the velocity in physical site units is

    v=v_gap/rho=2J(b+nu)sin(pi nu)>0.                    (12)

Expanding (8) at fixed nu gives

    E0(L,M)=L f(rho)-pi v/(6L)+O(L^-3).                 (13)

This is the c=1 Casimir coefficient. It alone is not the tower
argument; the zero modes and oscillator levels can also be read
from the complete spectrum (7).

First change M by delta M at fixed L near the minimizing density.
Equations (10),(12) give

    Delta E_charge = [pi v/(2K L)](delta M)^2+o(L^-1),
    K=rho^2=(b+nu)^(-2).                                (14)

If rho L is not integral, replace delta M by the integer's offset
from rho L and subtract the smallest such energy when quoting gaps.
No new mode is created by this rounding.

For the current mode, boost the filled gap Fermi interval by j
wave-number slots. Its total q changes by 2pi rj/M. Equations
(6),(7) give allowed physical momenta

    P=2pi(tM-rj)/L mod 2pi,

for an integer t. The displacement of the center in the cosine is

    (2pi j+P)/M=2pi(bj+t)/L mod 2pi/M.

Put D=bj+t. As the boost and the momentum branch vary, D runs through
the integers. Expanding the filled-sea energy gives

    Delta E_current=2pi v K D^2/L+o(L^-1),
    P=2pi rho D mod 2pi + O(L^-1).                     (15)

For example D=1 can be taken with the unboosted symmetric sea and
P=2pi M/L. Its exact gap is

    [2J sin(pi nu)/sin(pi/M)] [1-cos(2pi/L)],             (16)

which is 4pi J sin(pi nu)/[(b+nu)L]+O(L^-3). This small
current gap explains the previously measured first marker-graph gaps.
Using the XX value K=1 in physical coordinates would miss a factor
(b+nu)^2.

Particle-hole excitations a fixed number of slots from the two Fermi
edges have the usual partition multiplicities and energies
2pi v(N_R+N_L)/L. Their change in (6) adjusts P by O(1/L), changing
the current displacement only at O(1/L^2). Including charge and
current modes, the leading finite-size energies are therefore

    E-L f(rho) = (2pi v/L)[-1/12
                  +(delta M)^2/(4K)+K D^2+N_R+N_L]
                  +o(L^-1).                            (17)

The physical momentum also has the usual charge-current and oscillator
terms of order 1/L. Exchanging the names R and L compensates the
Slater sign convention used in Section 2. Thus (17) is one compact
boson, with central charge one and the parameter (14), rather than
additional quadratic defect states. Here the free-fermion continuum
identification of the partition towers is the stated standard input.

As nu approaches 0 or 1, the velocity vanishes. Those endpoints
are excluded from the conformal statement, consistently with the
isolated-defect quadratic spectrum of mobile_selector.md.

## 5. Physical-site entropy, not entropy of the gap variables

Let phi_XX(eta)>0 be the normalized periodic gap-chain ground
wavefunction at r particles. It is invariant under cyclic shifts.
The normalized physical marker wavefunction is

    Psi({x_j})=sqrt(M/L) phi_XX(eta).                    (18)

Any reference marker gives the same value. Its norm follows from
the M-to-one cover in (1). The gap amplitude is the absolute value
of the filled-sea Slater determinant, so (18) gives finite data
without diagonalizing the physical graph.

The script computes reduced density matrices using actual subsets
of markers in a physical interval. It does not substitute the gap
entropy. Sectors with different marker numbers in the interval
form separate Schmidt blocks. This allows a physical bipartition
of 72 sites with 57915 allowed full configurations and a largest
Schmidt block of only 500 by 500.

For b=4, nu=1/2, J=1, the half-chain results are

| M | physical L | S(L/2) | XX control S(M/2) |
| :-- | :-- | :-- | :-- |
| 8 | 36 | 2.00379104576 | 1.04373803047 |
| 12 | 54 | 2.14188009444 | 1.17547582153 |
| 16 | 72 | 2.23878096696 | 1.27019133100 |

The successive slopes 3 Delta S/Delta log L give c=1.0217085 and
1.0104996. The known XX controls give 0.9747161 and 0.9877102 at
the same cell counts. These are finite-size estimates approaching
one from opposite sides, not exact equalities of the two entropies.

Independently, Casimir fits at M=40,48,64,80,96,128,160 give
c=0.999999704 for b=1,4,9, with velocities (12) and an L^-4
correction included in E/L. At M=160 the current gaps reproduce
K to relative error below 7e-5, and the neighboring-M charge
stiffness reproduces 1/K to below 0.3 percent. The complete
physical graph spectra at six (b,M,r) values, dimensions 24 to
315, agree with (7) within 1.7e-14, including asymmetric fillings
where the twist sign matters.

Reproduce with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 length_fluid_checks.py --max-M 16

The output is length_fluid_checks.json. These two central-charge
estimators test the exact PROJECTED marker model. They cannot be
reported as two estimators for the full microscopic candidate.

## 6. Implication and limit

The all-ring marker mechanism now has a complete first-order spectrum,
not just a compatible-length subsequence or a ground-energy formula.
The local composition in selector_wrapper.md produces this graph with
the logical identity on the transported payloads. This removes a
specific finite-size bookkeeping gap.

The first-order payload couplings still depend on r_p r_p', and finite
u,lambda generate higher-order terms and cell leakage. The present
calculation does not prove that their infrared limit is the direct
sum of this boson and the desired payload CFT. It also does not
construct the missing high-channel initialized core. Both tasks
remain necessary for a microscopic counterexample.
