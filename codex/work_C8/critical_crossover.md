# An exact massive crossover with a critical first-order projection

Status: proved soluble control, with two calibrated numerical diagnostics.
This is a four-state nearest-neighbor chain. It is not the C8 candidate
and is not a counterexample to the c-d conjecture.

## 1. Literal microscopic bond

Each physical site has two qubits. Let XI=X tensor I and similarly for
other two-qubit products. For Delta>0 and 0<epsilon<=Delta/4 put

    a=(Delta/2)(II-IZ)-(epsilon/2)ZI
                                  +(epsilon/2)(YX+XY),
    h=(a tensor II+II tensor a)/2
                                  -(epsilon/2)XZ tensor XI.    (1)

Use the same Hermitian h on every bond of a periodic ring. The on-site
massive Hamiltonian is Delta(II-IZ)/2. Its complete local zero space
is the first qubit tensored with |0> on the second qubit, and its gap
is Delta. The interaction has a bounded norm independent of ring length.
Product ZZ parity is an exact symmetry.

[proved] The first-order compression is the critical transverse-field
Ising chain

    P H P=-(epsilon/2)sum_j (X_j X_(j+1)+Z_j).                 (2)

The mixing terms in (1) have zero P compression; the IZ in XZ restricts
to +1. The projected theory has one massless Majorana, velocity epsilon,
and c=1/2 under the standard free-Ising CFT identification.

## 2. Exact bands at fixed nonzero coupling

[proved] Choose four on-site Majorana matrices

    a_j=XI, b_j=YI, alpha_j=ZX, beta_j=ZY,

and attach the product ZZ parity string on all previous sites. Equation
(1), apart from the constant M Delta/2, is the quadratic Hamiltonian

    H=(i/2)sum_j [epsilon a_j b_j+epsilon b_j a_(j+1)
                 +Delta alpha_j beta_j
                 +epsilon a_j alpha_j-epsilon b_j beta_j].    (3)

For example i b_j a_(j+1)=-XZ_j XI_(j+1), including its string.
The on-site terms follow from the Pauli multiplication table. The ring
seam has the usual parity-dependent Jordan-Wigner twist.

The Hermitian single-particle Bloch matrix is iA(k), with

    A(k)=[[0,t,epsilon,0],[-t*,0,0,-epsilon],
          [-epsilon,0,0,Delta],[0,epsilon,-Delta,0]],
    t=epsilon(1-exp(-ik)).

After ordering the two bipartite sets as (a,beta) and (b,alpha), its
positive energies are the singular values s_-,s_+ of

    T(k)=[[t,epsilon],[epsilon,-Delta]].

Set x=sin^2(k/2). Then

    s_-^2+s_+^2=Delta^2+2epsilon^2+4epsilon^2 x = F(x),
    s_-^2 s_+^2=epsilon^4
                    +4Delta epsilon(Delta epsilon+epsilon^2)x = G(x).

In particular G(x)>0 for every real k at every epsilon>0. There is
no zero-energy bulk mode.

[proved] For epsilon<=Delta/4 the minimum is at k=0 and the bulk gap is

    gap=(sqrt(Delta^2+4epsilon^2)-Delta)/2
        =epsilon^2/Delta+O(epsilon^4/Delta^3).                  (4)

For completeness, let lambda(x)=s_-^2. Differentiating
lambda^2-F lambda+G=0 gives

    lambda'=(G'-F'lambda)/(F-2lambda).

The denominator is positive in this range. The numerator is positive
because G'/F'=Delta^2+Delta epsilon, while
lambda<=|t|^2+epsilon^2<=5epsilon^2<=5Delta^2/16.
At x=0 the two singular values can be read from the real symmetric
T(0), giving (4). Thus the gap is the global minimum, not merely the
value at a sampled momentum.

This is the disordered massive side of the Ising point. The low
Schur complement is t+epsilon^2/Delta at zero frequency, so the
on-site Ising field increases while the inter-site coefficient stays
fixed at leading order. There is no massless conformal sector at any
of the fixed couplings in (4). The antiperiodic free vacuum belongs to the even fermion sector. Its
energy is compared directly to the spin ground energy at the finite
lengths tested below. No finite-size uniqueness assertion is needed
for the absence of a massless bulk mode.

The mass is also visible directly in second order. YX+XY couples
|00> only to |11>, with matrix element 2i, and annihilates |10>.
The virtual energy at epsilon=0 is Delta. Consequently

    H_eff=-(epsilon/2)sum(XX+Z)
           -(epsilon^2/(2Delta))sum(I+Z)+O(epsilon^3).          (5)

The coefficient in (5) agrees with the exact low-band gap (4).
The O(epsilon^2) relevant mass is compatible with uniform convergence
to (2) at fixed slow time epsilon*t. Resolving it takes much longer
physical times, of order Delta/epsilon^2.

## 3. Two diagnostics and microscopic checks

[numerical] Reproduce with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 critical_crossover.py

The script first diagonalizes the actual 4-by-4 Bloch matrix on 4097
momenta at epsilon/Delta=.25,.1,.02,.005, Delta=1. The two positive
bands agree with the singular-value formulas to below 1.6e-15. Their
minimum agrees with (4). Independent dense diagonalization of the
literal spin bond (1) on L=3,4,5, dimensions 64,256,1024, matches the
free-fermion ground energies below 3.1e-16 at epsilon=.1.

The entropy calculation uses the ground-state covariance symbol
C(k)=-sign(iA(k)), on 32768 half-shifted momenta, and diagonalizes
its restriction to ell consecutive physical sites. If its eigenvalues
are nu_r, then

    S(ell)=(1/2)sum_r h((1+nu_r)/2),
    c_S(ell)=3[S(2ell)-S(ell)]/log 2,
    h(p)=-p log p-(1-p)log(1-p).

The energy calculation uses the antiperiodic free-vacuum energy E(L)
and the bulk-canceling combination

    c_E(L)=4L[E(2L)-2E(L)]/(pi epsilon).

For a massless mode with velocity epsilon this tends to its c. In the
massive model these are apparent central-charge diagnostics, not actual
CFT central charges. They do not assume that the massive model has a
conformal velocity.

Set the mixing terms to zero for the calibration. This leaves the
critical Ising chain plus an inert massive qubit. At epsilon=.02,
c_S(128)=.4999461 and c_E(512)=.5000007, matching c=1/2.
Restoring the mixing at the same epsilon gives

| size | c_S(size) | c_E(size) |
| :-- | :-- | :-- |
| 8 | .388948 | .504089 |
| 16 | .277357 | .491574 |
| 32 | .135105 | .452822 |
| 64 | .031449 | .338250 |
| 128 | .001854 | .139201 |
| 256 | not evaluated | .015465 |
| 512 | not evaluated | .000134 |

The entropy at ell=256 is already .975540, approaching a constant.
At epsilon=.005, c_E(8)=.503842 and c_S(8)=.481175, but both decrease
as sizes grow. The crossover length grows as Delta/epsilon. The full
input and output values are in critical_crossover_checks.json.

## 4. Scope of the obstruction

[proved] A positive additive on-cell gap, bounded neighboring-cell
couplings, an exactly critical first-order Hamiltonian, and convergence
of dressed local dynamics at fixed slow time do not imply a CFT at
fixed nonzero coupling. Equations (1)-(5) are an explicit example.
Even two small-size estimates near the reference c would not fix this
order-of-limits problem.

This does not refute the C8 candidate. That candidate was designed to
exclude the low-dimension resonances that are allowed in this Ising
control. Its Gaussian stability conditions must still be combined with
a valid low-frequency reduction. The example only excludes replacing
that step by first-order criticality or the one-step dynamics estimate.
