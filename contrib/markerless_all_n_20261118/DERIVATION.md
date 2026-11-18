# An arbitrary-n logical band family on nine physical states

2026-11-18. **[Proved finite-dimensional statement; internally checked,
not independently reviewed.]** For every finite n>=1 the modified literal
nine-state bond below has a leading quadratic logical symbol with exactly
2n simple unit-circle roots and distinct sums of all different unordered
pairs of roots, including repeated-root pairs. Parameters may depend on n,
but are fixed as the physical chain length L grows. This is not a theorem
of a microscopic unitary compact Virasoro CFT.

The modification is necessary to this proof: it uses a biased vacuum clock,
an independently slower active rail, and an L1-only field. It does **not**
prove the original unbiased, uniform-charge-field certificate schedule at
every n. Its physical alphabet and express-rail mechanism are unchanged.

## 1. Literal bond and physical controls

Start with [bond9.py](../markerless_d9_d10/bond9.py), whose nine letters are

    L0,L1,R0,R1,W,A,B,XA,XB.

Keep every express, switch, conversion and port rule. Make these changes
to its repeated nearest-neighbour bond:

1. Replace the two W moves by

       -q sum_a |L_a W><W R_a| + h.c.,  0<q<=1/4.

2. Replace the four active moves by

       -w sum_(s,a) z_s(a)|L_a s><s R_a| + h.c.,  0<w<=1/2.

3. Add the diagonal nearest-neighbour bias

       -sum_(y in {R0,R1,W,A,B,XA,XB}) |y W><y W|.   (1)

4. Set Delta=3 and replace the old conserved-charge field by the
   nonnegative physical on-site field h_L P_(L1).

Use K=0, or any K>=0, and on-site head energy nu_head. Express and switch
amplitudes and the bias magnitude define the unit energy scale. For the
selected n fix the original matched phases

    phi_A=pi/(2n), phi_B=0, psi_A=pi, psi_B=pi/2,
    port phases -pi/4,+pi/4.                         (2)

No coefficient or phase is changed on wrong-length cells. `model.py`
defines all 81-by-81 entries, and recovers the original bond exactly when
q=w=1, bias=0, and both fields vanish. The bond remains Hermitian, one-site
translation invariant and microscopically nearest-neighbour. It conserves
the original charge q(L_a)=q(R_a)=2a, q(A)=q(XA)=-1, q(B)=q(XB)=1,
q(W)=0, and the head number.

The original exact grammar proof still applies: changing amplitudes does
not change any off-diagonal move, and (1) and the field are diagonal.
Every forbidden LR-pair count and every boundary projector

    B_i=[letter_i is not L][letter_(i+1) is not R]

commutes with every term. No finite-K compression error is introduced.
For a legal cell L* H R*, (1) is exactly -|0><0| when H=W and is zero
on a charged cell, independently of its neighbours' bit or head states.
This is a physical boundary operator, not a freely assigned logical clock.

## 2. Clock, readout and zero-port selector

A cell with m data sites has physical length m+1 and W-clock positions
0,...,m. At zero field its clock is

    H_v(m)=-q adjacency(path_(m+1))-|0><0|
          =Q_m-(1+q^2)I,
    Q_m=diag(q^2,1+q^2,...,1+q^2)-q adjacency.        (3)

Let E_m be the lowest eigenvalue of Q_m. Including m=0, where E_0=q^2,

    (1-q^2)^2 q^(2m+2) <= E_m <= q^(2m+2),
    gap(Q_m)>=1-2q >=1/2,  m>=1.                    (4)

Indeed Q_m=B^T B for the upper bidiagonal B with q on its diagonal and
-1 on its superdiagonal. The Frobenius norm of B^(-1) is at most
q^(-m-1)/(1-q^2), giving the lower bound. The trial q^j gives the upper
bound. Interlacing after deleting the first coordinate gives the gap.
The normalized positive ground vector obeys v_j=q^j[1+O(q^2)] at fixed m.
There are exactly 2^m degenerate vacuum bit states at h_L=0.

The *actual* projected L1 field on a selected n-data cell is

    h_L sum_(i=0)^(n-1) gamma_i occupation_i,
    gamma_i=sum_(j=i+1)^n |v_j|^2
           =q^(2i+2)[1+O(q^2)].                     (5)

The upper summation endpoint is n, not the n+1 of the eleven-state clock.
The endpoint port visibility is zeta=|v_n v_0|^2, not the unbiased alpha^4.
`checks.py` projects the literal word operator and verifies both identities.

For completeness, the clock retains a zero-port ground-length selector.
Fix target n and choose

    nu_head=1+q^2-(n+2)E_n,
    e_*=-E_n,
    gamma_sel=E_(n+1)/(n+2).                         (6)

The vacuum cell energy is E_m-(n+2)E_n. Its excess over e_*(m+1) is

    delta_m=E_m+(m-n-1)E_n.

Write r=q^2/(1-q^2)^2<=16/225. The bounds (4) imply
E_(m+1)/E_m<=r. For m<n, with k=n-m, delta_m>=[r^(-k)-(k+1)]E_n>E_n.
For m=n+1 it equals E_(n+1); for m>=n+2 it is at least (m-n-1)E_n.
Dividing by m+1 gives a wrong-length density gap >=gamma_sel.
This argument includes zero-data W cells, unlike the marker construction.

At |g|,|t_d|<=1/4, a charged-cell hopping row sum is at most
2+|g|+|t_d|. With Delta=3 its energy is positive at (6), so it also loses.
Headless segments have nonnegative energy; they lose to e_*<0. In sectors
with LR violations, cutting at conserved boundaries leaves either the same
one-head cells or headless L*R* segments. Their only extra energy is the
nonnegative K penalty and readout. Thus at lambda=0 (6) selects the target
length on compatible rings, including comparisons with different numbers
of heads. Positive h_L does not spoil the comparison: it is a positive
operator and the all-zero vacuum trial achieves the previous minimum.
Nonzero-port, bulk-subtracted sector selection is a separate question.

## 3. Carrier, phases and exact quadratic formula

The selected cell vacuum energy, excluding nu_head, is
e_v=E_n-1-q^2. Set a=3-e_v>=4. After the A parity gauge and the express
phase gauge u(X_j)=exp(-ij phi_A), both species see the same carrier:

    D_0,...,D_n, X_0,...,X_n.                        (7)

There are exactly 2n+2 vertices. Diagonals are a on D and a+D on X.
Edges are -w along D, -1 along X, -1 from D_0 to X_0, and

    H_(X_n,D_n)=-i t_d, H_(D_n,X_n)=+i t_d.

The loop flux is psi_s-n phi_s=pi/2 for both species. The sources are
X_n on the left cell and D_0 on the right cell. Only D_0,...,D_(n-1)
convert; the inactive D_n stays in every inverse. Neither the source nor
the closing edge is moved to a fictitious marker vertex.

Write R=(H tensor I+I tensor H)^(-1), x=R|X_n,D_0>,
U_i=(|D_i><D_i| tensor I)x and V_i=(I tensor |D_i><D_i|)x. In units
lambda^2 g^2 zeta the leading quadratic blocks are

    T=-2 Re(U^dagger R V),
    M=-2i Im(U^dagger R U+V^dagger R V).             (8)

The CAR signs and the overall phase are those of the literal matched-flux
source, not adjustable logical phases. The changes in (1)-(5) only change
the scalar propagator and source visibility; the two conversion insertions
retain the same Jordan-Wigner strings. On the entire two-cell Fock code
their coefficient is scalar plus a number-conserving quadratic operator.
The literal n=1,2 calculation tests this on all bit states, not just the
one-particle block, and also verifies zero first port order, a scalar
lambda^2 g^0 coefficient and zero lambda^2 g^1 coefficient.

For all real t_d, T is real and even in t_d, and M is imaginary Hermitian,
antisymmetric under transpose and odd in t_d. Thus

    h(k,-t_d)=conj(h(-k,t_d)),
    h(k,t_d)=T e^(ik)+T^T e^(-ik)+M+F.              (9)

The total sum of its root angles is zero modulo 2pi: the constant and
leading coefficients of det[T z^2+(M+F)z+T^T] are both det T. This does
not require equality of individual opposite-branch pair sums.

An exact finite-D derivative follows without differencing. Put H=H_0+i t_d J,
J_(X_n,D_n)=-1, J_(D_n,X_n)=1, J_2=J tensor I+I tensor J, and y=-R J_2 x
at t_d=0. Restrict y to U_1,V_1 just as x was restricted. Then

    B_D=U^T R U_1-U_1^T R U-U^T R J_2 R U
       +V^T R V_1-V_1^T R V-V^T R J_2 R V,
    M'(0)=-2i B_D, T'(0)=0.                         (10)

If z_j is a tree-root null vector and v_j^speed is its signed nonzero
Bloch speed, the derivative of the sum of the two continued opposite
roots is -2 z_j^dagger M'(0)z_j/v_j^speed. Both numerator and denominator
must be transformed under any whitening congruence.

## 4. The markerless delay and loop scalings

Let p=n+1 and A=aI-w adjacency(path_p). Define the low-rail filters

    f_i(s)=integral_0^s [exp(-A(s-r))]_(0i)[exp(-Ar)]_(i0) dr,
    Gamma_ij=integral_0^infinity f_i(s)f_j(s) ds,
    Q_ij=integral_0^infinity f_i'(s)f_j(s) ds.        (11)

Gamma>0 on the n active coordinates: their lowest time powers are distinct,
w^(2i)s^(2i+1)/(2i+1)!. Integration by parts gives Q^T=-Q.

For the full tree, eliminating the p-vertex express block C_D gives the
exact endpoint identity

    [(H_0+z)^(-1)]_(X_n,D_i)
       = beta_D(z)[(H_0+z)^(-1)]_(D_0,D_i),
    beta_D(z)=1/det(C_D+z).

The endpoint cofactor is one because there are n express edges and one
unit switch. Squaring this identity turns the left return filter into
a convolution with a positive delay kernel of transform
product_r [c_r/(c_r+z)]^2, where c_r are C_D's eigenvalues. Its mean is

    tau_D=2 sum_r 1/c_r=2p/D+O(D^(-2)).

Expanding this delay and the low-rail self-energy gives

    D^(2p) T(0)=-Gamma+S/D+2p Q/D+O(D^(-2)),        (12)

where S is real symmetric. The expansion is at each fixed n,w,q.
This is the exact delay method of [C33](../../codex/work_C33/PROOF.md),
applied with p express vertices and the actual n+1 low vertices. The
normalized entries are rational in the finite carrier parameters and
regular in 1/D near zero, so the stated remainder can equivalently be
obtained by finite block elimination. The coefficient 2p is 2(n+1),
not the eleven-state 2(n+2).

The loop scaling also changes. Define

    eta=t_d D^n=t_d D^(p-1).

Eliminating a source at X_n gives, to leading order,

    D^(-p)(|D_0>-i eta |D_n>).

The first path crosses the n express edges and the switch; the second
uses the direct closing link and one high-energy denominator. On the
low product graph put R_0=(A tensor I+I tensor A)^(-1),
x_0=R_0|0,0>, x_n=R_0|n,0>, and define U_0,V_0,U_n,V_n as above. Then

    B=U_0^T R_0 U_n-U_n^T R_0 U_0
      +V_0^T R_0 V_n-V_n^T R_0 V_0,
    partial_eta[D^(2p)M](0)=2i B+O(D^(-1)).         (13)

The V difference is zero by time reversal on the real low graph; retaining
it makes (13) directly checkable. For i<j a shortest path in the first
U term visits 0,i,j,n in order and has n hops. Its three resolvent factors
give w^n/(2a)^(n+3). The reversed ordering needs 2(j-i) extra hops.
Thus B_ij=w^n/(2a)^(n+3)+O(w^(n+2)) for every i<j.

## 5. First nonzero opposite-pair splitting

Normalize energies by

    W=lambda^2 g^2 zeta D^(-2p),
    h_L=W epsilon/q^2,
    F=epsilon diag(gamma_i/q^2).                    (14)

This is a choice of the physical L1 coefficient. Let rho_j be the distinct
generalized eigenvalues of (F,Gamma), and x_j their real Gamma-orthonormal
vectors. Suppose 0<rho_j<1/2. At D=infinity the roots are
+/-acos(rho_j/2). The two imaginary perturbations in (12)-(13) are

    4ip Q sin(k)/D,  2i eta B.

Their mixed eigenvalue shift, divided by the speed 2sin(k) and summed
over the two opposite roots, gives

    partial_eta sigma_j(0)=A_j/D+O(D^(-2)),
    A_j=-16p sum_(l != j)
             (x_j^T Q x_l)(x_j^T B x_l)/(rho_j-rho_l).   (15)

Here sigma_j is the **sum**, not the half-centre. The symmetric S in (12)
cannot contribute at this order: its first vector correction is real and
has zero expectation against imaginary antisymmetric B. The same is true
of the O(D^(-1)) correction in (13) on an unperturbed real vector.
The strict Gram limit alone would miss (15), because a real null vector
has zero expectation against B. Finite express delay is essential.

The exact derivative (10) and direct finite-D roots verify the factor
-16(n+1), as well as both powers of D, at n=2,3. No higher-n root scan
is used in the proof.

## 6. Nonvanishing and separation for every n

Set q=(wu)^2. With D_w=diag(w^(2i)),

    D_w^(-1) Gamma D_w^(-1)=C_a(w)>0,
    D_w^(-1) F D_w^(-1)=epsilon diag(u^(4i)[1+O(q^2)]).

At fixed w>0, principal-minor expansion of the generalized characteristic
polynomial gives distinct eigenvalue valuations u^(4j), j=0,...,n-1.
The leading coefficients are positive because all principal minors of
C_a(w) are positive. The limiting vectors are its reverse Gram-Schmidt
vectors. This also follows successively by eliminating the largest
diagonal field and taking a positive Schur complement.

For the subsequent w->0 coefficient use c=2a, omega=w/c, U=c u.
Then Gamma=c^(-3) Gamma_bar, Q=c^(-2) Q_bar, B=c^(-3) B_bar and
q=(omega U)^2. After the congruence diag(omega^(2i))^(-1), the limits are

    C_ij=(2i+2j+2)!/[(2i+1)!(2j+1)!],
    Q^0_ij=C_ij(i-j)/[2(i+j+1)],
    B_bar,scaled=omega^(6-3n)[J+O(omega^2)].        (16)

The unique largest skew entry is the last pair (n-2,n-1); its coefficient
is one by the shortest-walk calculation. The inactive D_n remains present.
Let v_j,d_j,L_j be the exact objects of [MOMENT_SIGN.md](MOMENT_SIGN.md).
For j<n-1 the sum in (15) has successive leading behavior

    [omega^(6-3n)/(epsilon c^2)] L_j U^(-4j),        (17)

and the last branch has the same expression with -L_(n-2) and exponent
-4(n-2). The limits are interpreted in this order: derive the U->0
coefficient at fixed omega, then expand that coefficient at omega->0;
choose a sufficiently small fixed omega first, and then U. Smooth
dependence a=4+q^2-E_n does not alter these leading coefficients.

The new sign theorem gives L_j<0 for every j. Choose small enough w so
all these leading coefficients keep their signs. Then choose small enough
u so the powers in (17) dominate their remainders. The first n-1 A_j
are positive and successively separated by their different valuations;
the last is negative. Thus all A_j are distinct for every finite n>=2.
Their sum is zero, as required by the exact determinant identity.

## 7. A finite, L-independent parameter schedule

This is an existence schedule with terminating finite algebraic tests;
it is not a moderate-coefficient estimate or a closed optimized formula.
All small/large choices below depend on n, never on the chain length L.

1. For n>=2 choose a sufficiently small dyadic w<=1/2 for the nonzero
   leading coefficients in (17). They are computable from finite low-rail
   resolvents and reverse Gram-Schmidt, and their limiting signs are proved.
2. Choose a sufficiently small dyadic u with q=(wu)^2<=1/4 so the A_j
   are distinct and the unscaled generalized eigenvalues r_j of
   (diag(gamma/q^2),Gamma), in decreasing order, satisfy r_(j+1)/r_j<1/32.
3. Choose a positive dyadic epsilon with

       epsilon <= min[1/(2r_0), sqrt(r_(n-1)/(4r_0^3))].

   At D=infinity the tree roots are simple. Every different unordered
   pair-sum comparison except two opposite-branch pairs has angular margin
   at least delta_nc=epsilon r_(n-1)/4. To check this, expand
   acos(epsilon r_j/2)=pi/2-asin(epsilon r_j/2). In a nonzero combination
   with coefficient one-norm <=4, the largest r term dominates the rest
   by >=29/32. The cubic arcsine remainder is bounded by
   epsilon r_(n-1)/8 using the displayed choice. Comparisons with an odd
   half-period are bounded away from zero separately. The same estimate
   bounds distinct-root separations. Only opposite-pair sums are all zero.
4. Choose a sufficiently large finite dyadic D. The simple roots and
   non-opposite-pair margins persist, Herm(-T)>0, and the exact finite-D
   derivatives c_j=partial_eta sigma_j(0) are distinct by (15). The clock
   gap and carrier gap remain positive. Fix D before adding the loop.
5. Let kappa=min_(i<j)|c_i-c_j|>0. On a small eta interval bound root
   first derivatives by B_1 and |sigma_j'''|/6 by B_3. Choose a positive
   dyadic eta small enough that

       4 B_1 eta < min(delta_2,delta_nc)/2,
       2 B_3 eta^2 < kappa/2,
       t_d=eta D^(-n)<=1/4.

   The opposite sums are odd in eta, so their separation is at least
   kappa eta/2. All other margins remain at least half their tree values.
   Exactly 2n simple unit-circle roots persist. Hermitian reciprocal
   symmetry and degree 2n exclude a lost or additional root. For n=1
   skip the empty pair-centre comparison and choose a nonzero small field;
   simple roots already separate its three unordered pair sums.
6. Choose nonzero g and then lambda small enough for the finite-cluster
   quadratic corrections to fit below the resulting positive band
   tolerance. Set h_L by (14), nu_head by (6), Delta=3 and K=0.

For an entirely mechanical version of steps 1-5, enumerate the stated
dyadic candidates and use increasing exact/interval precision to test the
strict finite-dimensional inequalities; all matrices are algebraic at
rational q,w,D,eta,epsilon, including the clock eigenprojection. The
asymptotics prove termination of the nested choices. These are scalar and
finite-matrix tests, not new large scans carried out in this contribution.
The derivative bounds follow by implicit differentiation on isolated roots;
they are finite and effectively enclosable. An explicit implementation of
the terminating interval schedule is still a useful next deliverable.

## 8. What has and has not been proved

The strongest established statements are the all-n strict moment sign,
the resulting arbitrary-n nonresonant **leading quadratic** band family
at physical d=9, and the exact-grammar biased zero-port selector. Literal
n=1,2 checks connect the carrier calculation to the Hamiltonian, and
independent dense high-precision roots test the finite-delay effect.

For each fixed n the strict finite band margins survive sufficiently small
quadratic corrections on finite cell clusters. Higher conversion and port
orders also produce nonquadratic operators. Neither their norm smallness
nor momentum mismatch proves an interacting compact unitary CFT. Uniform
thermodynamic elimination, nonzero-port sector selection, interacting
stability, and conformal identification are still required. The new
readout is not the conserved charge, so a mobile lift would also need its
own field-dressing analysis. No microscopic unbounded-c theorem follows.
