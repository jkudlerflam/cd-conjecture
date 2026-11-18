# Exact structure and analytic anchors

All spin operators are Pauli matrices, j starts at zero, and the field is
counted once at **every** site. Set a=9/16, b=7/16, delta=a-b=1/8. These
derivations and their independent dense/sparse implementations are by the
same AI author; they are not an independent human verification.

## Symmetries and negative fields [proved algebraically]

On the infinite chain or an even periodic ring let P=product Z_j, T be a
one-site translation, K complex conjugation in the Z basis, and
S=product_(j odd) Z_j. At every h:

* [H,P]=[H,T]=0. The continuous rotation about Z is **not** a symmetry:
  the bond changes the number of down spins by two.
* A=S K is an antiunitary symmetry, A^2=1, with T A T^-1=P A. In a parity
  p=+1/-1 sector A sends momentum k to -k or pi-k, respectively. The
  latter relation enforces equal energies at k and pi-k when distinct;
  it is not evidence for two independent gapless fields.
* With site inversion I and R=product exp(-i pi Z_j/4), B=I R K is an
  antiunitary inversion symmetry. Ordinary I or K alone is not a symmetry.
* C=(product X_j)K obeys C H(h) C^-1=H(-h). C sends (X,Y,Z) to
  (X,Y,-Z), hence M and staggered Z change sign, while the literal bond,
  distance-two XX-YY, and kappa_j=X_j Y_(j+2)-Y_j X_(j+2) are unchanged.
  This relation holds on open chains and on odd as well as even rings.
* On a bipartite chain there is also a **unitary** field reversal
  F=product_(even j) X_j product_(odd j) Y_j. It changes Z and kappa to
  their negatives and leaves XX-YY unchanged. At h=0 it is an additional
  symmetry. The two field-reversal maps differ by A and a harmless phase;
  their different choices of chirality vacuum are consistent.

Consequently E(h)=E(-h), M(h)=-M(-h) on corresponding states, and gaps,
entropy, order magnitudes and central charges are even in h. Negative
fields require no new scan. Open chains retain P,A,B,C,F (with the stated
field reversal for C,F), but not T; odd rings do not retain the bipartite
definitions of A,F. At h=0 the thermodynamic Neel vacua can be exchanged
by F or T. No completeness theorem for all possible nonlocal symmetries
is claimed.

The order quantum numbers in the original frame are:

| Extensive staggered observable | P | A | T |
| --- | --- | --- | --- |
| N=sum (-1)^j Z_j | even | even | odd |
| D=sum (-1)^j b_j, literal bond b_j | even | even | odd |
| B2=sum (-1)^j (X_j X_(j+2)-Y_j Y_(j+2)) | even | even | odd |
| K2=sum (-1)^j kappa_j | even | **odd** | odd |

Chiral order breaks A and one-site translation separately, preserving
their product and P. It is distinct from Neel or ordinary bond order.
An antiunitary symmetry does **not** give an ordinary even/odd block of
energy eigenstates: A-odd matrix elements may be imaginary. The usable
spectral selection rule for K2 is **same parity, momentum shifted by pi**.
Low opposite-parity levels alone cannot diagnose this Ising transition.

## Real frame and boundary conditions [proved algebraically]

Let U=product_j exp(-i j pi Z_j/4) and H_R=U^dagger H U. Away from a seam,

    H_R = sum_j [delta (S+_j S-_(j+1) + S-_j S+_(j+1))
                    + (-1)^j (S+_j S+_(j+1) + S-_j S-_(j+1))
                    + Z_j Z_(j+1)] - h sum_j Z_j.

Here S+=(X+iY)/2, not a Pauli X. Thus Jxx=(delta+(-1)^j)/2 and
Jyy=(delta-(-1)^j)/2. The field and Z correlations are unchanged.

For L divisible by four this is the periodic real Hamiltonian. For
L=2 mod 4 **both** transverse coefficients on the closing bond get an
extra minus sign. Ignoring that sign changes the ring. An open chain has
no twist and retains full fields on both ends. The ED constructs the
unitary phase explicitly, so the seam is not guessed. The transformed
translation operator is U^dagger T U; translation in the two-site written
Hamiltonian must not be confused with physical one-site translation.

    U^dagger kappa_j U = -kappa_j,
    U^dagger (XX-YY)_(j,j+2) U = -(-1)^j (XX-YY)_(j,j+2).

Thus the source table's uniform real-frame B is minus the original
staggered B2 per site; its chiral m is minus the original staggered
chirality. Their **magnitudes** are unchanged. A becomes K in this frame.

## Exact fermion description and high-field anchor

With n_j=(1-Z_j)/2 and c_j=(product_(l<j) Z_l) S+_j, on open chains:

    H = sum_j [ i delta (c†_(j+1)c_j-c†_j c_(j+1))
                 + i(c†_j c†_(j+1)-c_(j+1)c_j) + 4 n_j n_(j+1)]
          + sum_j (2h-2 degree(j)) n_j + (# bonds)-hL.

The degree is two in the bulk, one at either open end. On a ring the
bulk number coefficient is 2h-4 and the constant L(1-h). The closing
fermion hop/pair has the additional boundary multiplier -P: even spin
parity gives antiperiodic fermions, odd parity periodic fermions, in the
untwisted original frame. The spin computation does not replace this by
one fixed fermion boundary condition. The real-frame gauge has its own
additional L-dependent seam described above.

This is an interacting pairing chain: V=4, pairing magnitude 1, imaginary
hop magnitude 1/8, and chemical potential 4-2h. It is the V4 cut with full
field h=7/4 when h=1.75. It is **not** the phase-1.53/p153c Hamiltonian.
No density, mode covariance, exponent or momentum mismatch from p153c
is used here.

For h positive and sufficiently large, divide H by 2h. The unperturbed
on-site projector sum has a unique product vacuum and unit gap, and the
finite-range bond perturbation tends to zero. The small-perturbation
gap theorem, e.g. [Yarotsky, Theorem 1](https://arxiv.org/pdf/math-ph/0411042),
therefore guarantees a unique gapped thermodynamic phase for all
sufficiently large h. This argument supplies **existence**, not a numerical
threshold such as h=1.8 or 1.9.

The all-up product is not an eigenstate at any finite h: each bond creates
an adjacent down-spin pair with nonzero amplitude i(a+b)=i. These distinct
states cannot cancel in its image. On a large ring the large-h expansion is

    e0(h) = 1-h - 1/[4(h-1)] + O(h^-3),
    M(h)  = 1 - 1/[4(h-1)^2] + O(h^-4),
    Delta_odd(h) = 2h-4-1/4 + O(h^-1).

The pair denominator 4(h-1) retains its ZZ cost exactly at second order;
these are asymptotic expressions, not formulas near the transition. The
leading one-down-spin dispersion is 2h-4+2 delta sin k. The gap is above
the dressed vacuum, not an exactly saturated vacuum. There is no finite
saturation field for the literal all-up/all-down product states.

## A nearby exact reference, not a substitute model

For delta=0 only, distribute half of the field onto each adjacent ring
bond. Its parallel block has eigenvalues 1 +/- sqrt(1+h^2), and its
antiparallel block has eigenvalue -1 twice. Therefore, for |h|<=sqrt(3),
each shifted bond is positive semidefinite and both Neel products attain
energy -L. For |h|<sqrt(3), the inequality

    H(delta=0)+L >= [2-sqrt(1+h^2)] sum_j (1+Z_j Z_(j+1))/2

gives a uniform gap above those two vacua on compatible rings. At
delta=1/8 the hopping acts on these products, so this certificate does
**not** prove the actual model's boundary or its gap. It explains why a
low-field Neel regime and a sensitive vicinity of sqrt(3) are natural
starting hypotheses. The actual low-field claim is supported by the new
finite-size doublet, bulk-gap and structure-factor data, not this changed
Hamiltonian's exact threshold.

## Candidate continuum theories and their limitations

The high-field phase and the Neel phase are massive; an entropy central
charge is not assigned to their finite-size vacuum doublets. A liquid
with one emergent compact boson has Virasoro c=1. There is no microscopic
U(1), so its stability requires all symmetry-allowed nonoscillatory locking
vertices to be irrelevant; this is not supplied by the entropy slope.

For the proposed chiral transition a useful **conditional** theory is

    H_b = (v_b/2pi)[K(∂theta)^2 + K^-1(∂phi)^2],
    H_I = -(i v_I/2)(xi_R ∂xi_R-xi_L ∂xi_L) - i m_I xi_R xi_L,
    H_int = lambda (∂phi) i xi_R xi_L + other allowed operators.

The density-times-Ising-energy coupling is symmetry allowed and marginal
at the decoupled point. Different boson and Majorana velocities are
allowed. [Alberton et al.](https://arxiv.org/html/1609.02599v2) establish in
their model that such coupling can give continuous or fluctuation-driven
first-order behavior, depending on velocities; their phase diagram is
not transferred to this Hamiltonian. The target velocities and coupling
remain unmeasured here.

If K2 has a nonzero overlap with the Ising spin times the boson identity,
its critical correlations contain (-1)^r A r^-1/4, its scaling dimension
is 1/8, and its same-parity momentum-pi level has leading energy
2pi v_I/(8L). The neutral Ising energy field has dimension one. However,
the lattice operator can also contain terms cos(Q r+phase) O_b O_I with
different dimensions. One must fit those oscillatory pieces explicitly
and establish the unmodulated overlap before reading a pure Ising eta
from an envelope. Our short-ring connected kappa correlators are retained
without an imposed power law; the source's long-distance correlators are
not in Git.

For several independently regulated critical components the entropy has
the conditional form S=sum_a (c_a/6) log xi_a+constant. Hence

    6 dS/dlog xi_max = sum_a c_a dlog xi_a/dlog xi_max.

Unequal constant velocities alone do not change the asymptotic total c,
but different finite-chi cutoffs, crossovers and branch changes can change
the measured slope. Transfer lengths in parity sectors are operator
lengths, not labeled boson/Ising regulators or measured velocities. The
observed ratio xi1/xi0 about 4.2 cannot be read as v_b/v_I. A two-boson
interpretation of c about 2 would additionally require two positive-speed
mode families and stability against allowed locking interactions. Neither
has been measured in the committed target ladders.

The c-theorem constrains suitable RG flows, not an arbitrary traversal
of h through a microscopic phase diagram with possibly unequal velocities.
It does not by itself rule out an extended c=2 region adjacent to another
critical regime. No such extension is established by the present data.
