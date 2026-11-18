# The genuine loop derivative and its finite-delay coefficient

2026-11-18. **[Proved below; internal numerical checks, no independent review.]**
This concerns the leading quadratic coefficient of the eleven-state
microscopic bond, not its interacting infrared limit. The all-n distinctness
criterion in Section 6 remains open; it is not included in the proved claims.

## 1. Physical conventions and exact parity

Use the literal modified C84 bond from
[the preceding contribution](../minimal_d_unbounded_audit_20261117/PROOF.md).
There are n active conversions and p=n+1 data vertices. The carrier is

    D_0,...,D_n, X_0,...,X_n,M,

with data edges -w, express/junction edges -1, and diagonal a on data
and a+D on express/marker vertices. Here a=3-e_v(p)>0, with the actual
biased vacuum energy e_v(p)=E_p-1-q^2. Delta remains exactly 3.
The inactive D_n and X_n remain. The actual physical phases are
phi_A=pi/(2p), phi_B=0, psi_A=pi, psi_B=pi/2 and ports -pi/4,+pi/4.
After the established gauge both species have

    H(t)=H_0+i t J, J_(M,D_n)=-1, J_(D_n,M)=1,
    t=t_d real.                                      (1)

Only the two existing data-to-marker transitions are restored. Neither
the physical L1 field nor the scalar source phase is replaced by a
freely chosen logical operator. The microscopic coefficients are fixed
across wrong-length cells as well as the selected length.

Remove the common positive factor lambda^2 g^2 zeta, where
zeta=|v_0 v_p|^2. With R=(H tensor I+I tensor H)^(-1), x=R|M,0>,
let U_i=P_i^(1)x and V_i=P_i^(2)x for i<n. Then

    T(t)=-2 Re[U^dagger R V],
    M(t)=-2i Im[U^dagger R U+V^dagger R V].            (2)

The exact physical field contributes h_L diag(gamma_i) before division
by the common factor, with gamma_i=sum_(j=i+1)^p |v_j|^2. All assertions
below keep that field fixed when differentiating t.

Since H(-t)=conj(H(t)), T is real and even in t, while M is imaginary,
Hermitian, antisymmetric under transpose, and odd in t. Consequently

    h(k,-t)=conj(h(-k,t)),
    h(k,t)=T(t)e^(ik)+T(t)^T e^(-ik)+M(t)+F.          (3)

For a simple reciprocal pair at t=0, label its continuations k_j,+
and k_j,- from +k_j and -k_j. Their sum sigma_j(t) is odd in t.
Also det[T z^2+(M+F)z+T^T] has equal constant and leading coefficients.
When its 2n roots remain simple on the unit circle, their total angle
sum is zero modulo 2pi. In particular sum_j sigma_j'(0)=0. This trace
constraint does not force any two individual pair sums to coincide.
For n=1, M is identically zero and the sole pair stays reciprocal for
every t; there is no comparison between two different Cooper pairs.

## 2. Exact finite-D first derivative

All matrices in this paragraph are real tree matrices. Set

    L=H_0 tensor I+I tensor H_0, R=L^(-1),
    J_2=J tensor I+I tensor J,
    x=R|M,0>, y=-R J_2 x,
    U_1,i=P_i^(1)y, V_1,i=P_i^(2)y.

Differentiation of the inverse gives R'(0)=-i R J_2 R and x'(0)=i y.
The exact antisymmetric matrix B_D is

    B_D=U^T R U_1-U_1^T R U-U^T R J_2 R U
       +V^T R V_1-V_1^T R V-V^T R J_2 R V,
    M'(0)=-2i B_D, T'(0)=0.                          (4)

This computes the imaginary derivative directly, without subtracting
large real onsite energies. It is valid for every finite n,D,w,q.
Let z_j be any normalized null vector of h(k_j,0), and let
v_j=z_j^dagger partial_k h(k_j,0) z_j be its nonzero signed speed.
Implicit differentiation and (3) give the exact pair-sum coefficient

    c_j=sigma_j'(0)=-2 z_j^dagger M'(0) z_j/v_j.      (5)

The corresponding half-centre derivative is c_j/2. Equation (5) must
use the physical kernel vector, not a Euclidean eigenvector of a
congruent pencil without transforming both numerator and denominator.

If the tree hopping is replaced by its symmetric Gram limit, z_j can
be real and z_j^dagger M'(0)z_j=0. Thus a nonzero M'(0) alone does
not prove a nonzero first centre shift. This is a real cancellation
in that limiting argument, not a finite-D symmetry obstruction.

## 3. The finite express delay supplies the missing imaginary vector

Let A=aI-w adjacency(path_p), retaining all p data vertices, and put

    f_i(s)=integral_0^s [e^(-A(s-r))]_(0i)
                         [e^(-Ar)]_(i0) dr,
    Gamma_ij=integral_0^infinity f_i(s) f_j(s) ds,
    Q_ij=integral_0^infinity f_i'(s) f_j(s) ds.        (6)

Gamma is positive definite on the active principal block; Q is real
antisymmetric, by integration by parts and vanishing endpoint values.

The exact tree transfer identity is
G_(M,i)(z)=beta_D(z) G_(0,i)(z), where beta_D(z)=1/det(C_D+z) and
C_D is the (p+1)-vertex express block. Thus the source-M return filter
is beta_D(0)^2 times the source-0 filter convolved with a positive
delay kernel b_D. Its Laplace transform is

    product_r [c_r/(c_r+z)]^2,

where c_r are C_D's eigenvalues. Its mean is
tau_D=2 sum_r 1/c_r=2(p+1)/D+O(D^-2), and its second moment is O(D^-2).
The resolvent identity and delay expansion are the same ones used in
[C33, Sections 3–4](../../codex/work_C33/PROOF.md), now with p data
vertices and only the first p-1 conversions retained.

The full source-0 filters differ from (6) by O(D^-1); their Gram
correction is symmetric. Expanding b_D*f=f-tau_D f'+O(D^-2) gives

    D^(2p+2) T(0)=-Gamma+S/D+2(p+1)Q/D+O(D^-2),     (7)

where S is real symmetric. The expansion is at each fixed finite n,w,q;
it does not assert a uniform-in-n error bound. It also includes the
scalar expansion of D^(2p+2) beta_D(0)^2 in S.
The normalized coefficients are rational functions of D with a regular
expansion at 1/D=0, by finite block elimination. Thus the displayed
first-order delay coefficient specifies an O(D^-2) remainder without
assuming uniform high-rail bounds on every time derivative.

## 4. The actual loop source and the low-rail matrix B

Use the scaled loop eta=t D^p. Eliminating the fast rail from a source
at M gives, to leading order,

    D^(-p-1)(|0>-i eta |n>).

The sign follows from H_(D_n,M)=+it. The long route to 0 and direct
route to the inactive endpoint n are both retained; deleting n destroys
this calculation. Additional high-rail visits have another inverse
power of D. On the low product graph set

    R_0=(A tensor I+I tensor A)^(-1),
    x_0=R_0|0,0>, x_n=R_0|n,0>,
    U_s,i=P_i^(1)x_s, V_s,i=P_i^(2)x_s, s=0,n,
    B=U_0^T R_0 U_n-U_n^T R_0 U_0
      +V_0^T R_0 V_n-V_n^T R_0 V_0.                 (8)

B is real antisymmetric. The V difference is zero on this real low
graph by time reversal in the unconverted first coordinate; keeping
it in (8) is useful for checking the resolvent calculation. Fast-rail
terms not represented in (8) are O(D^-1) after the normalization below.
The direct resolvent expansion gives

    partial_eta [D^(2p+2)M](0)=2i B+O(D^-1).         (9)

In particular the limiting loop matrix has no independent phase
freedom beyond the physical phases in (1).

## 5. First nonzero mixed coefficient

In the D^(2p+2) normalization, set the physical field to the previous
schedule h_L=W nu/q^2, W=lambda^2 g^2 zeta D^(-2p-2). Thus

    F=nu diag(gamma_i/q^2).

Let rho_j be the distinct generalized eigenvalues of (F,Gamma),
ordered decreasingly, and choose real x_j with
F x_j=rho_j Gamma x_j and x_j^T Gamma x_l=delta_jl. Assume
0<rho_j<1/2. The limiting roots are +/-acos(rho_j/2).
Write Q_jl=x_j^T Q x_l and B_jl=x_j^T B x_l.

Second-order eigenvalue perturbation, once in the imaginary part of
(7) and once in (9), gives

    partial_eta sigma_j(0)=A_j/D+O(D^-2),
    A_j=-16(p+1) sum_(l != j) Q_jl B_jl/(rho_j-rho_l).
                                                        (10)

Indeed the two perturbations are 4i(p+1)Q sin(k)/D and 2i eta B.
Their mixed energy shift is
16(p+1)eta sin(k) sum_l Q_jl B_jl/[D(rho_j-rho_l)].
Dividing by the limiting speed 2sin(k) and summing the two root
shifts gives (10). The real symmetric S in (7) cannot contribute
at this order: its kernel correction is real and its expectation
against the antisymmetric B is zero. The O(D^-1) correction to (9)
has the same property on an unperturbed real vector.

The coefficient with respect to the physical t is D^p times (10).
All these coefficients are held at fixed q,w,nu,D when the loop is
varied. Taking D=infinity before computing the mixed coefficient
would erase the leading answer. sum_j A_j=0, consistently with (3).

## 6. A rational all-n criterion and an analytic partial family

The following reduction is exact at the level of asymptotic
coefficients. It does not assume the desired all-n sign conclusion.
Use q=(wu)^2, as before. For fixed small w and u tending to zero,
the generalized eigenvalues have distinct valuations u^(4j).
The limiting vectors are reverse Gram-Schmidt vectors: x_j has no
components before j and is Gamma-orthogonal to every later coordinate.

To state the subsequent w->0 coefficient economically, put
c=2a, omega=w/c, U=c u and N=nu c^3. Then q=(omega U)^2 exactly.
Time y=ct gives Gamma=c^-3 Gamma_bar, Q=c^-2 Q_bar and B=c^-3 B_bar,
where the barred low Hamiltonian has diagonal 1/2 and hopping omega.
This is a change of variables in the proof; physical Delta remains 3
and a tends to 4. Removing the diagonal powers omega^(2i) gives

    C_ij=(2i+2j+2)!/[(2i+1)!(2j+1)!],
    Q^0_ij=C_ij (i-j)/[2(i+j+1)].                    (11)

For i<j the shortest walk contributing to (8) visits 0,i,j,n
in order. It has n data hops; the reverse ordering has at least
n+2(j-i) hops. Integrating its three positive time intervals gives

    B_ij=w^n/(2a)^(n+3)+O(w^(n+2)), i<j.             (12)

After the corresponding diagonal congruence, the unique largest entry is the last
off-diagonal pair (n-2,n-1). Up to a common positive factor,

    B_bar,scaled=omega^(6-3n)[J+O(omega^2)],
    J_(n-2,n-1)=1, J_(n-1,n-2)=-1.                  (13)

Let v_j have v_j,j=1, v_j,i=0 for i<j, and (C v_j)_i=0 for i>j.
Put d_j=v_j^T C v_j>0. These rational vectors are mutually
C-orthogonal. Define the rational coefficients

    L_j=sum_(l>j) (v_j^T Q^0 v_l)(v_j^T J v_l)/d_l,
    j=0,...,n-2.                                    (14)

The barred normalized vectors are v_j/sqrt(d_j), and the leading
generalized eigenvalues are N U^(4j)/d_j. In the physical Gamma metric,
Q_jl gains a factor c and B_jl a factor omega^(6-3n).
Thus the sum in (10), for j<n-1, has leading hierarchy

    [omega^(6-3n)/(nu c^2)] L_j U^(-4j),             (15)

provided L_j is nonzero. Terms with l<j have smaller inverse powers.
The last branch has coefficient

    -[omega^(6-3n)/(nu c^2)] L_(n-2) U^(-4(n-2)).

The last two branches therefore have opposite leading signs. For
full distinctness it suffices that L_j !=0 for j=1,...,n-2; L_0
need not be nonzero, since its response divided by the common factor
stays bounded on this scale. The ordering is: first choose small
positive omega so the nonzero limiting coefficients persist, then
choose small U so their different valuations dominate the remainders.

There is one closed all-n coefficient. Write r=n-2. Then
v_(n-1)=e_(n-1), v_(n-2)=e_r-(C_(r,r+1)/C_(r+1,r+1))e_(r+1),
and

    C_(r,r+1)/C_(r+1,r+1)=(r+1)/(4r+5),
    L_(n-2)=Q^0_(r,r+1)/C_(r+1,r+1)
            =-1/(16n-12).                            (16)

Consequently **for every finite n>=2**, sufficiently small w and then
u give two distinct extreme centre responses, separated also from
all remaining responses. For n=2,3 this already proves that all
responses are distinct. It does not separate the remaining n-2
responses when n>=4.

Exact rational elimination at n=4 gives

    (L_0,L_1,L_2)=(-80879/328737072,-5033/2227500,-1/52).

In particular L_1 !=0 completes the n=4 case analytically. The other
small values are n=2: (-1/20), and n=3: (-95/32368,-1/36).
moments.py reproduces these fractions using rational arithmetic.
An explicit witness for the only extra n=4 condition is

    v_1=(0,1,-92/225,49/825),
    v_2=(0,0,1,-3/13), v_3=(0,0,0,1),
    d_2=900/13, d_3=3432,
    L_1=-49/25740-2576/7239375=-5033/2227500.

They establish finite cases of the asymptotic construction, not an
all-n extrapolation from numerical roots.

**[Open.]** Prove L_j !=0 for all n and 1<=j<=n-3, or find a structural
zero and analyze its next coefficient. All displayed L_j are negative;
negativity is a conjecture here, not an invoked total-positivity theorem.
Even a zero in (14) would be a failure of this leading asymptotic
criterion, not a proof of an exact physical pair-centre identity.

For n=1,...,4 choose w, then u, then the previous sufficiently small
nu, then a sufficiently large finite D. Equations (10)-(16) give
distinct finite-D c_j, except the vacuous n=1 comparison. The already
proved simple-root and non-Cooper margins persist. Finally take a
sufficiently small nonzero t. Analyticity and distinct c_j then give
strict separation of **all** unordered pair sums in the leading
quadratic symbol, including repeated roots. This is an analytic finite
family on the same eleven-state alphabet, not an arbitrary-n theorem.
