# The joint ring constraint for a moving payload chain

Status: proved finite-volume spectral reduction for the equalized
first-order model and checked with a free XX payload. No finite-coupling
microscopic infrared theorem is claimed.

## 1. Arbitrary translation-invariant logical payload

Use the physical marker Hilbert space of length_fluid.md and attach
a D-dimensional logical state to every cell. The reference cover is

    (x,eta,z_0,...,z_(M-1)),

with x in Z_L and sum eta=r=L-bM. Its quotient is

    (x,eta,z) ~ (x+b+eta_0,S eta,S z).                   (1)

The dimension is (L/M) binomial(M,r) D^M. A marker move transports
the cell states by the identity, and the logical bond is independent
of adjacent lengths, as supplied by helix_equalizer.md. Let H_pay(M)
be any translation-invariant logical ring Hamiltonian and choose
simultaneous energy/translation eigenvectors with

    H_pay chi_a=E_a chi_a,
    chi_a(S z)=exp(2pi i m_a/M) chi_a(z),
    m_a in Z_M.                                         (2)

Degeneracies are included by repeating a. At fixed physical momentum
P=2pi ell/L, the same gap gauge as in length_fluid.md gives a uniform
gap hopping flux P. H_pay is unaffected by this gauge. The quotient
(1) now requires

    m_gap+m_a+ell=0 mod M.                              (3)

Consequently the COMPLETE first-order spectrum is

    E=kappa M+E_a-2J sum_(j=1)^r
                         cos(q_j+2pi ell/(LM)),          (4)

where q_j=2pi[I_j+(r-1)/2]/M for distinct I_j in {0,...,M-1},
m_gap=sum_j I_j+r(r-1)/2 modulo M, and ell ranges over
0,...,L-1 subject to (3). This is an exact finite-dimensional
statement, with no continuum approximation.

Proof. The marker reference shift rotates BOTH gap and payload
coordinates, which explains the extra phase m_a in (3). The
marker Hamiltonian in the cover acts only on eta and x, while the
logical Hamiltonian acts only on z. They commute, including at the
reference bond, since a global cyclic shift commutes with H_pay.
Fourier transform in x, the diagonal gap gauge, and Jordan-Wigner
on the gaps diagonalize the first part exactly as before. Applying
the quotient condition selects (3). These are unitary decompositions
of the entire invariant Hilbert space, proving completeness. In
physical coordinates a marker crossing the origin cyclically
rotates the canonically ordered payload labels. Ignoring this
rotation would give a different spectrum.

Equation (4) is more informative than saying the two systems are
decoupled. Their local Hamiltonians separate, but their allowed
global momentum sectors are linked. In particular a generic
payload excitation cannot be combined with an arbitrary gap
eigenstate at the same chosen physical momentum.

## 2. One free XX payload

Take a spin-1/2 logical XX ring with hopping amplitude g>0. In its
sector of n up spins the free wave numbers are

    p_j=(2pi/M)[A_j+(n-1)/2],
    E_a=-2g sum_j cos(p_j),
    m_a=sum_j A_j+n(n-1)/2 mod M.                       (5)

Using (5) in (3),(4) gives every energy in a physical Hilbert space
of dimension (L/M) binomial(M,r) binomial(M,n). The positive filled
sea in both sectors has m_gap=m_a=0, so it lifts at P=0 and attains

    E0=kappa M-[2J sin(pi r/M)+2g sin(pi n/M)]/sin(pi/M). (6)

Each term individually has its lowest energy there. Thus (6) is
the joint ground energy at those charges, without an assumption
that momentum minimization leaves a product state available.

At r/M=n/M=1/2, the two physical velocities are

    v_len=2J(b+1/2),    v_pay=2g(b+1/2).

The Casimir term is -pi(v_len+v_pay)/(6L), the sum of two c=1
sectors. The two independent sets of Fermi-edge particle-hole
partitions are present in (4): their total momentum changes the
allowed ell through (3), but ell can always be chosen. The induced
change of gap flux only shifts zero-mode sectors and does not
remove either oscillator family. The charge and current energies
are quadratic near interior fillings. Thus the standard free
continuum identification gives two bosons, with their zero-mode
quantization coupled by (3).

At J=g and kappa=2(J+g)/pi, both physical velocities coincide and
the bulk minimum has r/M=n/M=1/2. On compatible rings

    E0=4JM/pi-4J/sin(pi/M),

so the single-velocity Casimir coefficient measures c_total=2.
The same mechanism adds one length oscillator to an n-flavor free
logical target. This claim concerns the first-order Hamiltonian;
it is not evidence that microscopic higher-order corrections
preserve all n+1 modes.

## 3. Independent physical-graph checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 joint_fluid_checks.py

The physical implementation labels marker positions in increasing
order and stores one payload bit per cell. Marker moves that wrap
from L-1 to 0 rotate the entire payload tuple. The XX term exchanges
payload bits on neighboring cells. This directly builds a matrix
in the physical basis, independently of (3),(4).

For (b,M,r,n)=(4,4,2,2), (3,5,2,2), and (4,6,3,3), with J=1 and
g=0.7, the full Hilbert dimensions are 162,340,1800. EVERY eigenvalue
agrees with (4),(5), with maximum errors 1.1e-14,1.2e-14,3.4e-14.
The asymmetric r=2,M=5 case is useful for detecting the twist sign.
The ground energies agree with (6) within 1.6e-14.

For J=g=1, b=4, M=40,48,64,80,96,128,160, a fit of E0/L against
1,L^-2,L^-4 at known velocity v=9 gives c=1.999999409. This is a
Casimir check of the projected joint model, not a second independent
estimator for any microscopic counterexample. Output is saved in
joint_fluid_checks.json.

## 4. Limit of the reduction

The equalizer and this calculation resolve the leading length/payload
mixing and the global ring bookkeeping. A fixed nonzero microscopic
coupling still produces leakage, additional local terms, and possibly
relevant or marginal perturbations of the joint theory. The two missing
pieces remain the initialized high-channel uniform core and an infrared
argument that survives those corrections. No whole-class bound or
verified counterexample follows from the reduction alone.
