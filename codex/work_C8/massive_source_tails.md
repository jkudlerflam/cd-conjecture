# A massive source can have an algebraic connected correlation

Status: proved exact control. It excludes an incorrect definition of a
fast kernel, not local massive-sector integration itself.

## 1. A model whose massive sector can be removed exactly

On each four-state site, use a logical qubit s and an ancillary qubit a.
Let

    H_dec=-(J/2)sum_j(X_s,j X_s,j+1+Y_s,j Y_s,j+1)
                              +(Delta/2)sum_j(1-Z_a,j),
    U_theta=product_j exp[(i theta/2)Z_s,j Y_a,j],
    H_theta=U_theta dagger H_dec U_theta.                      (1)

J,Delta>0. Equation (1) is a short explicit repeated nearest-neighbor
Hermitian bond, obtained by conjugating each bond and each split on-site
term by the same on-site unitary. It preserves total Z_s. It has a
critical XX sector and exactly independent massive ancillas after the
known on-site change of basis. Its spectrum and spatial entanglement
across site cuts agree with those of H_dec. In particular its massless
sector is the usual free c=1 theory, and the ancillary excitation cost
is Delta.

One may set J=epsilon and theta=epsilon/Delta. Then H_theta is a weak
bounded nearest-neighbor perturbation of Delta sum n_a, with the XX
chain as its first-order retained Hamiltonian. No spectrum of a large
cell or unspecified effective interaction is used in this control.

## 2. The connected bare-source correlation is not fast

[proved] The elementary on-site identity is

    U_theta X_a U_theta dagger
                    =cos(theta) X_a+sin(theta) Z_s Z_a.        (2)

The bare operator X_a is off diagonal between the zero and excited
spaces of the unperturbed ancillary Hamiltonian. Nevertheless its
connected Euclidean correlation in the exact thermodynamic ground
state of H_theta is

    C_X(tau)=cos^2(theta) exp(-Delta tau)
                         +sin^2(theta) C_Z(tau), tau>=0,       (3)

where C_Z is the connected local Z_s correlation of the half-filled
XX chain. To prove (3), transform both insertions using (2). The mixed
terms have zero ancillary vacuum expectation. The X_a term creates
one ancilla and hence has factor exp(-Delta tau). The Z_a term acts
as +1 in the ancillary vacuum. At half filling <Z_s>=0, so subtracting
the one-point functions does not remove the second term.

[proved] Wick contraction for the XX chain gives

    C_Z(tau)=4 a(tau)^2,
    a(tau)=(1/pi)int_0^(pi/2) exp[-2J tau cos(k)] dk.           (4)

Both particle and hole propagators equal a(tau) by the half-filled
particle-hole symmetry. Z_s=1-2n_s supplies the factor four.
Writing x=pi/2-k, the only unsuppressed endpoint as tau grows is x=0.
Since sin x=x+O(x^3), endpoint Laplace asymptotics yields

    a(tau)=1/(2pi J tau)[1+O((J tau)^(-2))],
    C_X(tau)=sin^2(theta)/(pi^2 J^2 tau^2)
                    [1+O((J tau)^(-2))]+O(exp(-Delta tau)).    (5)

The asymptotic can also be bounded directly by splitting the integral
at a fixed small x, Taylor bounding sin x there, and bounding the other
interval exponentially. Thus no positive exponential rate bounds the
full connected bare-source correlation when sin(theta) is nonzero.
The microscopic model nevertheless has an exactly removable massive
sector with gap Delta.

## 3. The necessary subtraction is an operator subtraction

In the decoupled frame, let P_a be the ancillary vacuum projector.
The retained part of the source (2) is sin(theta) Z_s P_a.
Subtract that operator, rather than just its expectation. Its remaining
source creates an ancilla with amplitude cos(theta); the extra term
proportional to Z_s(1-P_a) annihilates the ground state. The resulting
two-source function is cos^2(theta) exp(-Delta tau).

Equivalently use the known dressed local projector
U_theta dagger P_a U_theta in physical variables. This construction
is completely local in the control. It illustrates what a valid
irreducible massive kernel must remove.

[proved limitation] Additive defect cost, small defect occupation, and
ordinary connected correlations of bare off-code sources do not imply
exponential temporal decay. A proposed massive-sector proof has to
separate retained propagation and the operator mixing in (2). The
word connected alone does not specify that separation. The earlier
global bare-Q resolvent obstruction and this local source example are
different failures of shortcuts. Neither is a no-go for the desired
low-frequency reduction.

## 4. Reproduction

[numerical finite and quadrature checks] Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 massive_source_tails.py

At J=.2, Delta=1, theta=.1, the literal repeated spin Hamiltonian on
L=4, dimension 256, agrees with the full on-site conjugation below
1.5e-15. Its exact total-Z_s commutator vanishes. Spectral sums for the
source correlator at tau=.1,1,10 agree with (3), using the independently
diagonalized finite XX reference, below 1.5e-15.

Thermodynamic quadrature of (4) at tau=1,10,30,100,300,1000 gives

    (pi J tau)^2 C_slow(tau)/sin^2(theta)=1.0000125007

at the largest time. There C_slow=2.52463e-8, while the massive
exponential is below floating-point range. These computations check
(2)-(5); they are not estimates of the C8 candidate's central charge.
