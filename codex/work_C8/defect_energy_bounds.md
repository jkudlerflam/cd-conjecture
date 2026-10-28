# Ground-energy and cell-length bounds from local leakage control

Status: proved energy and occupation estimates. They do not identify
an infrared CFT. These estimates use no upper bound on the on-cell
energy and no finite truncation of the direct sum over cell lengths.

## 1. A block-diagonal interaction may be arbitrarily critical

Use the atomic setting of local_dressing.md. Let

    H=A+W+F, A=sum_j A_j, A_j>=Delta(1-P_j),
    W=sum_X W_X, F=sum_X F_X,
    P_X=product_(j in X)P_j, Q_X=1-P_X.

The bounded Hermitian interactions may have long supports, provided

    w=sup_j sum_(X contains j)||W_X||<infinity,
    f=sup_j sum_(X contains j)||F_X||<infinity.

Assume each W_X is block diagonal with respect to P_X plus Q_X.
F is the possible leakage. Define the complete reference Hamiltonian
H_ref=P(W+F)P on P=product P_j, and E_ref=inf spec H_ref.
No gap, uniqueness, or frustration-free property of H_ref is assumed.

[proved] Put g=Delta-2(w+f). If g>0, then

    E_ref-M f^2/g <= E0(H) <= E_ref.                           (1)

For a normalized state with energy at most E_ref+E_ex, E_ex>=0,

    <A>/M <=4 Delta f^2/g^2+2 Delta E_ex/(g M),
    <sum_j(1-P_j)>/M <=4 f^2/g^2+2 E_ex/(g M).                (2)

In particular E_ex=0 applies to every ground state. The estimates
hold for density matrices as well. Finite-volume ground states exist
in the C8 fixed-length sectors; the infimum formulation also covers
other atomic models.

## 2. Proof, including the comparison with a critical reference

Decompose F_X into its two diagonal blocks and its off-diagonal part.
Let H_keep=sum_X P_X(W_X+F_X)P_X. It commutes with each P_j.
In a sector whose defect set is B, it is the reference interaction
with all terms meeting B removed. Fill the missing sites with any
code state. The total norm of the removed terms is at most

    sum_(X meets B)(||W_X||+||F_X||)<=(w+f)|B|.

The variational principle for the full reference therefore implies

    H_keep >= E_ref-(w+f)N, N=sum_j(1-P_j).

This argument also applies to entangled states of the remaining code
sites, by filling with a density matrix and taking its expectation.
The Q_X diagonal blocks obey

    sum_X Q_X(W_X+F_X)Q_X >= -(w+f)N,

because Q_X<=sum_(j in X)(1-P_j). Thus the entire diagonal part is
at least E_ref-2(w+f)N. This proves the necessary comparison without
inverting a critical reference Hamiltonian or assuming a global band.

For delta>0 and f>0, Cauchy-Schwarz and Young's inequality give

    |<Q_X F_X P_X+P_X F_X Q_X>|
       <=(delta/f)||F_X|| <Q_X>
                          +(f/delta)||F_X|| <P_X>.

Sum over X. The first sum is at most delta<N>, and the second at
most M f^2/delta. The f=0 case is immediate by continuity. Since
N<=A/Delta, for 0<delta<=g we have the quadratic-form inequality

    H >= E_ref+[(g-delta)/Delta]A-M f^2/delta.                 (3)

Set delta=g for the lower energy bound in (1). The upper bound uses
a reference ground state as trial state. Set delta=g/2 in (3) and
use the stated energy upper bound to obtain (2). Every step involves
bounded interactions and positive forms; it remains valid for
unbounded A and infinite-dimensional cell spaces.

## 3. Conserved-charge constraints must be retained

[proved] If C is a conserved sector and its all-code reference space
is nonempty, let E_ref^C be that reference ground energy and put

    dE_C=E_ref^C-E_ref >=0.

The global lower bounds above still hold after restriction to C.
For a ground state in C, replace E_ex in (2) by dE_C. For a state
E_ex above that constrained ground energy, replace it by dE_C+E_ex.

One cannot silently set dE_C to zero. Filling a missing site in the
proof can change the conserved charge, so E_ref in that proof is the
unconstrained reference energy. This distinction matters for fixed
physical length in C8.

## 4. Consequences before and after the one-step dressing

[proved] Applied directly to A+epsilon V, take W=0, F=epsilon V.
Then f=|epsilon|J and g=Delta-2|epsilon|J. The unconstrained ground
energy differs from the first-order reference by at most

    M epsilon^2 J^2/(Delta-2|epsilon|J),

and the mean bare defect occupation is at most

    4 epsilon^2 J^2/(Delta-2|epsilon|J)^2.                     (4)

This is a local occupation statement. The total overlap with the
all-code subspace can still vanish exponentially with M.

[proved] Apply the same theorem to the exact dressed Hamiltonian
from local_dressing.md. There W=epsilon D, F=epsilon^2 R(epsilon),
w<=|epsilon|J, f<=epsilon^2 J_R, J_R=||R(epsilon)||_0. Its complete
all-code reference is

    H_ref,dressed=epsilon PVP+epsilon^2 P R(epsilon)P.

At g_d=Delta-2(|epsilon|J+epsilon^2 J_R)>0,

    0<=E_ref,dressed-E0(H)<=M epsilon^4 J_R^2/g_d,
    <N>_dressed/M<=4 epsilon^4 J_R^2/g_d^2.                    (5)

Thus a single local dressing gives a fourth-order error in the ground
energy density relative to its exact dressed P block. The P block
includes all the powers contained in this one conjugation. Equation
(5) is not an all-order block diagonalization.

For a constrained sector, add the dE_C terms of Section 3. Without a
bound on the constrained-versus-unconstrained dressed reference
energy, the fourth-order occupation estimate is not a claim for every
fixed-length physical ring. The bare estimate below has an explicit
such comparison for the selected reference density.

## 5. Bare C8 estimates at the physical reference density

Use the unperturbed atomic gap Delta_at and J<=2B_flat. Fix M and L
with r=L-bM in [0,M]. The exact first-order spectrum in joint_fluid.md
allows a direct estimate of the constraint cost. The payload ground
state on M cells can be chosen with cell momentum m_pay in Z_M.
At zero length flux the length ground sea has m_gap=0 and energy

    E_len(M,r)=-2 sin(pi r/M)/sin(pi/M), M>=2.

Choose ell in {0,...,M-1} such that ell+m_pay=0 modulo M. This satisfies
the physical quotient, because L>=bM>=M. The length momenta then acquire
phi=2pi ell/(LM), while the payload state is unchanged. The symmetric
length sea has zero sine sum, so its energy increases by at most

    phi^2/sin(pi/M)<=2pi^2/(b^2 M).

The minimum length energy over r differs from E_len(M,r) by at most

    pi^2(r-M/2)^2/(2M).

Indeed compare to -2/sin(pi/M), use 1-cos x<=x^2/2 and
sin(pi/M)>=2/M. A flux cannot lower the hard-core hopping ground
energy below zero flux, by applying the variational principle to the
absolute values of its components. This supplies the matching
unconstrained lower reference bound in every physical momentum sector.
The marker chemical potential and payload ground energy cancel in
the difference.

Consequently the actual first-order constraint excess is bounded by

    dE_C <= epsilon C_b(M,r),
    C_b(M,r)=pi^2(r-M/2)^2/(2M)+2pi^2/(b^2 M).                (6)

Here epsilon>0 and E_ref includes its factor epsilon. For
M=round[L/(b+1/2)], sufficiently large L makes r admissible and
|r-M/2|<=(b+1/2)/2. Hence

    C_b(M,r)<=[pi^2(b+1/2)^2/8+2pi^2/b^2]/M.

Using (2), the bare mean defect occupation in the ground state of
these legal fixed-M physical sectors is therefore at most

    4 epsilon^2 J^2/(Delta_at-2epsilon J)^2
       +2 epsilon C_b(M,r)/[(Delta_at-2epsilon J)M].           (7)

The last term is O(epsilon b^2/M^2), for the one fixed candidate b.
The estimate concerns these chosen marker sectors. It does not prove
that the full microscopic ground state chooses this M after dressing,
or compare it to illegal-selector sectors.

## 6. Quantitative control of arbitrarily long cells

In the C8 atomic representation, for ell outside {b,b+1},

    A_ell>=gamma_wrong(1+|ell-b|).

For any R>b+1 this implies the operator bound

    1_(ellhat>=R)<=A/[gamma_wrong(R-b)].

Thus any state to which (2) applies obeys

    (1/M)sum_j <1_(ellhat_j>=R)>
       <=[4 Delta f^2/g^2+2 Delta E_ex/(gM)]
                                    /[gamma_wrong(R-b)].      (8)

Before dressing, use f=epsilon J and the constraint cost in (6).
After dressing, (5) gives the stronger bound when its reference
constraint cost is controlled. No dimension factor appears in (8).

For a cell-translation-invariant state, let rho_Z be its reduced state
on a fixed finite cell interval. Let Pi_R project every cell there to
length less than R, and p_R=tr(Pi_R rho_Z). The union bound gives
1-p_R<=|Z| times the right-hand side of (8). If p_R>0, the normalized
truncated reduced state satisfies

    ||rho_Z-Pi_R rho_Z Pi_R/p_R||_1<=2 sqrt(1-p_R).

This follows by projecting a purification and using the trace distance
of the two pure states, then taking a partial trace. For states without
cell-translation invariance, the occupation estimate is spatially
averaged; no individual-site bound is asserted.

This proves local-state tightness as the length cutoff increases.
It is not a global-norm truncation, a spectral approximation at O(1/L),
or a bound on low-frequency response kernels. The critical crossover
control shows that even much stronger small ground-energy errors do
not alone identify the infrared phase.

## 7. A stronger grand-canonical length tail from one-site transfers

[proved, unconstrained atomic ground state] There is a factorial tail
when the perturbation changes an individual cell length by at most one.
This uses an additional property of the C8 marker transfer, beyond
bounded interaction norm. Suppose [A_j,ellhat_j]=0, and for ell>=b+2,

    A_j|ell >= gamma(ell-b),

while A_j has a zero vector at a selected length. Let v be the actual
local interaction norm of H-A, including its coupling. If v<=gamma/2,
then any unconstrained ground state obeys, at each cell j and every
integer R>=b+2,

    ||1_(ellhat_j>=R)psi||
          <= (v/gamma)^(R-b-1)/(R-b-1)!.                      (9)

This is independent of volume and local internal dimension.

Proof. Remove cell j and every interaction touching it. Call the
remaining Hamiltonian H_rest and its ground energy E_rest. The sum
of removed interactions has norm at most v. Filling site j with a
zero vector gives E0(H)<=E_rest+v. With Q_R=1_(ellhat_j>=R), positivity
of the local energy and the interaction norm imply

    Q_R(H-E0)Q_R >= [gamma(R-b)-2v]Q_R.

All unperturbed terms preserve ellhat_j. The perturbation's length
bandwidth one therefore gives

    Q_R H(1-Q_R)=Q_R(H-A)1_(ellhat_j=R-1).

Apply Q_R to the ground-state equation and invert its strictly
positive compressed operator. The result is

    ||Q_R psi|| <= v/[gamma(R-b)-2v] ||Q_(R-1)psi||.

For v<=gamma/2 the denominator is at least gamma(R-b-1). Iterating
from R=b+2 and using ||Q_(b+1)psi||<=1 proves (9). The same argument
applies to every vector in a degenerate ground space. It requires
no spectral gap in H_rest.

The complete C8 atomic Hamiltonian has this length bandwidth: fields
and ports preserve length, and marker moves transfer one site. Thus
(9) applies to its unconstrained atomic version, for v<=gamma_wrong/2.
It gives a stronger cutoff estimate than (8) in that setting.

It cannot be copied to the physical fixed-L quotient. Filling cell j
in the proof can change total length, and the best remaining-cell
energy in a fixed charge sector depends on the removed length. A
uniform canonical chemical-potential bound would repair that proof,
but it has not been established here for the dressed physical model.
For physical fixed-length sectors the proved statement remains (8),
with the explicit reference constraint cost of Section 5. The two
statements have different hypotheses.
