# Multivelocity Gaussian currents: local modes and a conditional c bound

Status: **[proved covariance algebra and conditional reconstruction/CFT
embedding theorems; microscopic Gaussian scaling remains unproved]**.

Author: the Codex GPT-6-family exact-alternative subagent, 2026-10-28.
The coordinator independently checked the matrix constants. No human
review or counterexample is claimed.

The stable multichannel reference covariance already determines a real
current commutator matrix of fixed signature and a positive quadratic
Hamiltonian. A *constant local linear change of fields* diagonalizes it.
No differential operator, inverse derivative, or nonlocal momentum
projection is needed. Different velocities therefore present no obstacle
to constructing the reference Gaussian current algebra. Transferring that
algebra to the microscopic chain, and placing it inside the conjecture's
full conformal scaling theory, remain explicit hypotheses.

## 1. Input and conventions

Use the conditional source comparison of
[IR_SOURCE_REVIEW.md](IR_SOURCE_REVIEW.md). It gives a limiting invertible
probe matrix U and a reference density covariance, up to fixed contact
terms. Let V=diag(v_b), with all v_b real and nonzero; let Zbar be the
invertible positive diagonal *bare* reference normalization. Let Lambda be
a real symmetric forward-coupling matrix. All quantities refer to one
fixed finite number N of branches.

Write

    D=(4 pi |V|)^(-1),   S=sign(V),
    B(omega,p)=D(-i omega+Vp)(i omega+Vp)^(-1).

After absorbing Zbar into the probes, the reference covariance is

    C(omega,p)=(I+B(omega,p)Lambda)^(-1) B(omega,p).       (1)

Indeed the original convention with Lambda_Z=Zbar^(-1)Lambda Zbar gives
C_ref=Zbar^(-1) C Zbar^(-1). Thus the physical probes are related to these
reference fields by the constant invertible matrix T=U Zbar^(-1).

Assume the explicit stability condition

    ||D^(1/2) Lambda D^(1/2)||_2 < 1.                     (2)

This condition is sufficient, not asserted necessary. It is available
at small enough coupling when the velocities stay away from zero.
The reference formula (1) is the input read in Proposition 5.8 of
[Porta--Scola--Singh](https://arxiv.org/html/2509.08665#S5.SS5), as recorded
in IR_SOURCE_REVIEW. The following matrix transformations are derived here.

## 2. Exact local Hamiltonian form of the covariance

Define

    M=I-D Lambda,       Nmat=I+D Lambda,
    K=M^(-1) S M^(-T)/(2 pi),
    H=(D^(-1)-Lambda D Lambda)/2,
    A=M^(-1) V Nmat.                                    (3)

**Theorem 1 [proved].** M is invertible, H is real symmetric positive
definite, K is real symmetric invertible with the same inertia as S,
and A=KH. The temporal-axis contact subtraction gives

    C(omega,p)-C(omega,0)
                  =p(i omega I+A p)^(-1)K.               (4)

Here C(omega,0)=-M^(-1)D for every nonzero omega. This subtraction fixes a
convenient contact convention; it does not remove a propagating mode.

**Proof.** Set E=D^(1/2)Lambda D^(1/2). The matrix D Lambda is similar
to E. Condition (2) makes I-E and hence M invertible. Also

    H=(1/2)D^(-1/2)(I-E^2)D^(-1/2)>0.                   (5)

Symmetry and inertia of K follow by congruence from S/(2 pi).
The identity M^(-1)D=D M^(-T), and VD=S/(4 pi), give

    K=2 M^(-1)V D M^(-T).

Moreover

    H=(1/2) M^T D^(-1) Nmat,

so KH=M^(-1)V Nmat=A. No commuting of Lambda with V or D has been assumed.

Multiplying I+B Lambda from the left by i omega+Vp gives

    (i omega+Vp)(I+B Lambda)=i omega M+p V Nmat.

Using this identity in (1), adding M^(-1)D, and combining the two
numerators gives

    C+M^(-1)D
       =2p(i omega M+pV Nmat)^(-1)V M^(-1)D
       =p(i omega I+A p)^(-1)K.

This proves (4). QED.

Thus K is the candidate equal-time current matrix and H is its positive
quadratic energy matrix. The velocities themselves are the eigenvalues
of KH, not in general the original diagonal entries v_b.

## 3. A constant transformation isolates every velocity

Let O be real orthogonal with

    O^T H^(1/2) K H^(1/2) O=diag(u_1,...,u_N)=Uvel.

The middle matrix is real symmetric and invertible. Its eigenvalues are
therefore real and nonzero, and Sylvester inertia implies

    number(u_j>0)=number(v_b>0),
    number(u_j<0)=number(v_b<0).                          (6)

Repeated velocities cause no difficulty: choose any orthonormal basis
inside their eigenspaces. There is no loss of residue rank there.

For reference fields n, define

    J=diag[(2 pi |u_j|)^(-1/2)] O^T H^(1/2) n.            (7)

For the physical scaling probes O_phys, replace n by T^(-1)O_phys.
Every matrix in this transformation is constant in space, time and scale.
Consequently (7) preserves locality and scaling dimension. It also avoids
the dimension change caused by applying a differential pole filter.

**Theorem 2 [proved].** In the contact convention (4), the transformed
covariance is diagonal:

    <J_j J_k>(omega,p)
       =delta_jk [sign(u_j)/(2 pi)] p/(i omega+u_j p).     (8)

**Proof.** A=KH is similar to H^(1/2)KH^(1/2). Conjugating (4) by
O^T H^(1/2) therefore yields

    p(i omega I+Uvel p)^(-1) Uvel.

The rescaling in (7) gives (8). QED.

In a periodic finite-band lattice reference, the number of positive and
negative transverse Fermi crossings is equal: the number of occupied
one-particle bands changes by minus one at a positive crossing and plus
one at a negative crossing, then returns to its initial value after one
Brillouin period. Thus N=2n and both counts in (6) equal n. The stable
interaction preserves this signature, even if magnitudes of left and
right velocities are not paired.

## 4. What covariance and positivity do determine

Let s_j=sign(u_j). With the inverse time Fourier factor exp(i omega tau),
the positive-time covariance obtained from (8) is

    <J_j(tau,p)J_j(0,-p)>
       = [|p|/(2 pi)] 1_(u_j p>0) exp(-|u_j p| tau),
       tau>0.                                           (9)

This follows by taking the pole at omega=i u_j p when u_j p>0.
It displays a positive spectral weight on the propagation ray
E=u_j p>0. In particular, its reflected-time quadratic form is an
integral of

    [|p|/(2 pi)] 1_(u_j p>0)
              |integral_0^infinity exp(-|u_j p|tau)f(tau,p)dtau|^2,

so it is nonnegative. Constant contacts at tau=0 do not change this
positive-time spectral statement.

[proved conditional Gaussian reconstruction] Suppose the physical
scaling probes exist as local Hermitian fields with tempered Euclidean
correlation distributions, all their moments obey Wick's rule with the
above noncontact covariance, and the usual positive-energy/positivity
reconstruction is available. Then their generated vacuum theory is the
product of N independent chiral Heisenberg current theories with speeds
u_j and positive level one after orienting each coordinate by s_j.

One direct proof constructs the oscillator Fock fields with covariance
(9). Wick's rule makes all inner products of polynomially smeared fields
equal to those of the given scaling fields. After quotienting null states,
this is an isometry on a dense generated vacuum space and extends to a
unitary equivalence. The oscillator commutators are therefore those of
the generated fields. With the space Fourier convention exp(-ipx),

    [J_j(x),J_k(y)]
        =-i delta_jk s_j delta'(x-y)/(2 pi),              (10)

as distributions. On a circle, the modes defined with coordinate s_j x
obey [j_m,j_n]=m delta_(m+n,0). The Sugawara modes

    L_m=(1/2) sum_l :j_(m-l)j_l:

have central charge one in each chirality. The speed enters time evolution,
not this central term. The Gaussian covariance alone would not imply the
operator commutator in an arbitrary non-Gaussian theory; the Wick/GNS step
is essential.

This construction determines the oscillator current algebra and its
positive signature. It does not determine the compactification lattice,
zero-mode spectrum, charged fields, full partition function, or the
complete microscopic scaling theory. Those distinctions matter for the
conjecture's compact-unitary-CFT premise.

## 5. Conditional embedding in a product of unitary conformal sectors

A full central-charge lower bound still needs an ambient conformal
interpretation. Here is one explicit sufficient version.

Assume independently that the full scaling theory has a unitary compact
conformal decomposition into finitely many chiral sectors. Sector alpha
has signed propagation velocity w_alpha != 0, positive energy and chiral
central charge c_alpha. Its time and space translations act with that
velocity. Assume local scaling fields have the usual discrete finite
expansion at fixed scaling dimension in tensor products of sector fields,
with positive two-point Gram forms. These are hypotheses about the full
limit, not consequences of the action recursion or (8).

The conventional nonanomalous total central charge used here is

    c_total=(sum_(right alpha)c_alpha
                 +sum_(left alpha)c_alpha)/2.            (11)

For an ordinary product of nonchiral CFTs this is exactly the sum of their
usual c. Allowing different left/right speeds also covers a tilted
Hamiltonian for a nonchiral CFT. If the intended definition permits only
paired equal speeds, that additional property must be checked separately;
it is not supplied by the generic reference covariance.

Suppose the physical limiting currents from (7) have scaling dimension
one, covariance (8), and Gaussian four-point functions as separated-point
distributions. Then

    sum_alpha c_alpha >= N,
    c_total >= N/2.                                      (12)

**Proof, first step: a spectral ray selects one velocity group.**
A nonzero positive-weight field in a chiral sector has positive spectral
measure on its chiral ray, with momentum magnitude ranging over an open
positive half-line. For a tensor-product field these measures convolve.
If q_alpha>=0 is the momentum magnitude contributed by its active sector,
then

    E=sum_alpha |w_alpha| q_alpha,
    p=sum_alpha sign(w_alpha) q_alpha.

For this convolution to be supported on one ray E=u p, each active factor
must satisfy |w_alpha|-u sign(w_alpha)=0, hence w_alpha=u. Otherwise
independently varying positive q_alpha produces spectral weight away from
that ray. Positive Gram forms prohibit cancellation of such off-ray
weights between orthogonal field components. A weight-zero vacuum factor
is harmless. Thus J_j belongs to the grouped chiral sector of velocity u_j.
This also covers a product of several fractional-weight fields with the
same velocity: the group has total weight one and its usual chiral stress.

The spectral-support argument requires the stated conformal decomposition
and expansion of fields. It is not being inferred from two-point poles in
a completely unspecified scaling theory.

**Second step: Gaussian four points make the visible currents abelian.**
Group all sectors with the same signed velocity u. Their total stress
has central charge C_u=sum_(alpha:w_alpha=u)c_alpha. Let r_u be the
multiplicity of u in (8). The r_u independent positive-norm weight-one
fields are chiral currents of this grouped theory.

For two such Hermitian currents A,B, the simultaneous simple-pole
coefficient in the four-point function of pairs (A,B),(A,B) is minus
the norm squared of their current-algebra bracket, divided by the squared
separation of the pair centers. A Wick four-point function has no such
coefficient. Positivity therefore makes the bracket zero, even if it
would point outside the selected current span. This is the argument proved
in section 8 of IR_SOURCE_REVIEW. Distinct velocity groups are already
separate tensor factors under the ambient hypothesis.

The r_u abelian currents have Sugawara central charge r_u. The difference
between the grouped stress tensor and their Sugawara tensor has central
charge C_u-r_u and nonnegative vacuum norm, so C_u>=r_u. Summing gives
sum_alpha c_alpha>=sum_u r_u=N, proving (12). QED.

This proof does not assume a single common stress tensor or common
velocity for the entire model. It uses a stress tensor for each declared
conformal velocity group. It therefore cannot be applied before that
ambient decomposition, or an adequate substitute for it, is established.

## 6. Precisely what is left for the microscopic candidate

The preceding algebra removes an obstacle at the reference level: under
(2), velocities are real, nonzero and of fixed signature, residues are
positive after local normalization, and every branch has a local mode.
It does not establish the microscopic hypotheses needed to use them.

1. The normalized marked-GN estimates in IR_SOURCE_REVIEW must still be
   proved for the full retarded action and induced sources. The source
   matrix must have the proved nonzero limiting singular value.
2. A fixed-momentum comparison must be strengthened to the local scaling
   distributions needed for Gaussian four-point identities or full Gaussian
   reconstruction. Exceptional momenta, contact terms and positivity cannot
   be omitted from that step.
3. For a statement about the conjecture's *full* c, prove the ambient
   unitary conformal decomposition or another theorem that embeds the
   current subtheory into the full low-energy theory with additive positive
   conformal central charges. The existence of a Gaussian current subtheory
   alone does not establish that all remaining low-energy modes are CFTs.
4. Compactification/zero-mode information is not contained in the current
   covariance. It must be supplied if compactness of the full limit is
   part of the asserted counterexample.

With those inputs, N=46 visible branches in the frozen d=23 construction
would imply c>=23>22. The inputs are not all proved. The present result is
a conditional reduction and positive-mode decomposition, not a solution.

## 7. Finite checks

Run `python3 contrib/decisive_attack_20261028/exact/multivelocity_checks.py`
with numpy available. It checks (3)--(8) for noncommuting symmetric Lambda,
unequal signed velocities, and nonscalar Zbar/source mixing, at N=2,4,6,8.
It also checks inertia and a deliberately repeated free velocity. The saved
`multivelocity_results.json` records errors and smallest positive eigenvalues.
These checks supplement the matrix proof and do not test microscopic RG.
