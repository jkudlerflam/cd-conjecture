# What transfers to eleven physical states

Started 2026-11-17; completed 2026-11-18. OpenAI Codex, GPT-6 family.
**[Proved below, not independently reviewed.]** The new statements are
an eleven-state biased code and uniform zero-port selector, a faithful
inactive-endpoint Gram limit, and an arbitrary-n simple-band/non-pair
margin. The **full** pair-sum margin of the tree remains exactly zero.
No all-n nonresonant loop theorem or microscopic unbounded-c result is
claimed. The source and finite-certificate distinctions are in the ledger.

## 1. Literal alphabet and modification

Keep C84's physical states, with no encoded-dimension substitution:

    L0,L1,R0,R1,DA,DB,XA,XB,E,MA,MB.                  (1)

There are eleven, independent of target n. A cell contains p=n+1 data
sites and one marker, hence b=n+2 physical sites and n payload bits.
The marker set is M={E,MA,MB}. The allowed pairs, U(1) charges,
conversions and ports are exactly those in C84. model.py independently
rebuilds its bond and recovers it entry for entry at the old settings.

Make the following explicit changes to the repeated nearest-neighbor
bond. None depends on the position in a cell.

- Give each vacuum move R0 R_a <-> L_a R0 and R0 E <-> L0 E amplitude
  -q, with 0<q<=1/4, instead of -1.
- Add the **actual boundary interaction**

      - sum_(m in M) |m R0><m R0|.                    (2)

- Optionally give each charged data-rail move D_s R_a <-> L_a D_s
  amplitude -w z_s(a), with 0<w<=1. Express and junction amplitudes
  retain unit magnitude. Keep the small closing edge t_d when studying
  the original loop; set it to zero only for the stated tree theorem.
- Use a nonnegative readout h_L on L1 alone, rather than equally on
  L1 and R1. The coefficient can be set to zero for selection/code work.
- Take Delta=3, D>=0, |g|<=1/2 and |t_d|<=1. Choose the marker
  coefficient below. Delta does not require a tiny E_n calibration.

All original phases can be kept: for the selected p, phi_A=pi/(2p),
phi_B=0, psi_A=p phi_A+pi/2, psi_B=pi/2, and port phases -pi/4,+pi/4.
They are fixed for the entire microscopic Hamiltonian, including
wrong-length cells. Onsite terms are split evenly between bonds.
The full finite K penalty is retained on forbidden pairs.

The bond remains one-site translation invariant, Hermitian, U(1)
conserving, and commutes with each on-site marker projector. At lambda=0
it also preserves each cell's excited-token count, including on illegal
words. Legal grammar is **not** exactly invariant at finite K: some
vacuum moves create a predecessor violation. This is not the exact-grammar
model used in the larger alphabet, and its splitting must be budgeted.

## 2. Why the clock bias is legitimate, and the naive replacement is not

In C84's legal vacuum words the first data site is R0 precisely when
the vacuum clock is at position zero. At later clock positions it is
an L letter. Thus (2) is -|0><0| on the clock, identity on all n bits.
It vanishes on every legal charged word. The preceding marker's type
does not matter, so it does not couple neighboring payloads at zero port.

The legal vacuum clock has positions 0,...,p and matrix

    H_v(p)=-q adjacency(path_(p+1))-|0><0|
          =Q_p(q)-(1+q^2)I,
    Q_p=diag(q^2,1+q^2,...,1+q^2)-q adjacency.         (3)

Q_p is the edge-projector clock with its terminal q^2 penalty. Its
lowest eigenvalue satisfies, for every p>=1,

    (1-q^2)^2 q^(2p+2) <= E_p <= q^(2p+2),
    gap(Q_p) >= 1-2q >= 1/2.                         (4)

For completeness, Q_p=B_p^dagger B_p with q on B_p's diagonal and
-1 on its superdiagonal. Its inverse has entries q^(-(j-i+1)) for
j>=i, so its squared Frobenius norm is at most
q^(-2p-2)/(1-q^2)^2. This gives the lower bound. The trial vector q^j
gives the upper bound. Deleting Q_p's first row and column gives a
Toeplitz matrix with minimum at least (1-q)^2; interlacing gives the
gap bound. See also the original
[biased-clock argument](../smaller_d_unbounded_20261116/PROOF.md).
The physical clock energy is e_v(p)=E_p-1-q^2<=-1. At each fixed p,
the normalized positive eigenvector satisfies v_j=q^j[1+O(q^2)].
The complete legal vacuum ground space remains exactly 2^n dimensional.

This construction does **not** assign a new physical head letter.
It also does not replace every old hop by a positive edge projector.
That naive replacement fails because R0 is a stored zero as well as
the clock head. For p=2, the bit-zero and bit-one vacuum restrictions
of those naive projectors are respectively

    Q_ff + q^2|0><0|,       Q_ff.                      (5)

The first has positive ground energy and the second has a zero mode.
At q=1/4 their ground energies are about 0.05835 and zero. This is
payload splitting at infinite K, not a finite-penalty correction.
The marker context in (2) avoids this error. It changes physical diagonal
coefficients; it is neither an irrelevant constant nor a relabeling.

## 3. Uniform zero-port selection, including finite K

Fix target p_*=n+1 and set

    mu_total=1+q^2-(p_*+2)E_(p_*),
    e_*=-E_(p_*), gamma=E_(p_*+1)/(p_*+2).             (6)

At lambda=h_L=0 a legal p>=1 vacuum cell has energy
mu_total+e_v(p)=E_p-(p_*+2)E_(p_*). Its excess over e_*(p+1) is

    delta_p=E_p+(p-p_*-1)E_(p_*).                     (7)

The ratio E_(p+1)/E_p is at most r=q^2/(1-q^2)^2<=16/225.
For p<p_*, setting m=p_*-p gives
delta_p>=[r^(-m)-(m+1)]E_(p_*)>E_(p_*).
For p=p_*+1 the excess is E_(p_*+1). For p>=p_*+2 it is at least
(p-p_*-1)E_(p_*). Therefore every wrong p>=1 has density at least
e_*+gamma. A zero-data cell has energy mu_total, not the fictitious
p=0 continuation of (3), and is also strictly above this bound.
Legal marker-free rings have energy density zero.

The statement survives finite K with a uniform bound that includes
all illegal strings and all token counts. In a cell with p>=1 data
sites, split off the offending vacuum moves as -q(a_j+a_j^dagger).
The inherited [C80 local-context calculation](../../codex/work_C80/PROOF.md),
also used by [C84](../../codex/work_C84/PROOF.md), gives partial isometries

    a_j a_j^dagger <= Pi_bad,(j-1),
    sum_j a_j^dagger a_j <= (p-1)I.                   (8)

All remaining terms preserve every bad pair. The new bias and readout
are diagonal, so do not alter this decomposition. Call that remainder
H_good, excluding K N_bad. For no bad pairs its minimum is
mu_total+e_v(p). Charged legal fragments have nonnegative energy after
subtracting marker energy: their hopping row sums are at most two per
head and their conversions at most |g|, absorbed by Delta=3.

For m>=1 bad pairs, cutting at them leaves at most m+1 legal fragments.
Each vacuum fragment has kinetic row sum at most 2q. The boundary bias
has norm at most one in the whole cell. Consequently

    H_good >= mu_total-1-2q(m+1)
             >= mu_total+e_v(p)-4q m,                (9)

where the last inequality uses e_v(p)<=-1 and m>=1. Combining with
the exact legal block proves (9) as a form inequality on all strings.
Young's inequality with K>4q now gives

    mu_total+e_v(p)-(p-1)q^2/(K-4q)
       <= E_K(p) <= mu_total+e_v(p).                 (10)

The upper trial is the legal zero-bit vacuum. For p=0 the minimum is
mu_total directly. On a marker-free ring the same bad-move argument,
without the bias, bounds the density below by -q^2/(K-4q).

Thus any finite

    q^2/(K-4q) < gamma/2                             (11)

selects the length b_*=n+2 among **all marker patterns** on compatible
rings at zero port. The selected marker sector is a physical ground
sector of the full finite-K model, not just a chosen legal compression.
Its wavefunctions can contain illegal components. A nonnegative L1
field preserves this ground-length comparison by the zero-bit trial.
At nonzero port, a common-bulk dressed-defect comparison is still
required; (11) does not establish it.

## 4. The finite-K code is a split band

The selected legal vacuum code at h_L=0 has gap at least 1/2.
Each charged head costs at least Delta-2-|g|>=1/2 before marker
energy, whereas the vacuum clock energy is at most -1. The extra
express mass only helps. The legal code thus has 2^n states and a
uniform auxiliary gap at fixed n.

At finite K, use H_K=K N_bad+B on a fixed one- or two-cell cluster.
For one cell, a sufficient finite bound is

    ||B|| <= B_p=(p+1)(Delta+D+|mu_total|+|h_L|+24).

For K>=8B_p and |z|<=2B_p, the illegal block inverse is at most 2/K
and its Schur correction is at most 2B_p^2/K. Taking that smaller
than one eighth of the legal gap isolates a 2^n-dimensional band.
Its splitting is O(B_p^2/K); it is not zero. The same finite-dimensional
resolvent and polar-overlap estimates control its source maps and
second-order coefficient. For two cells replace B_p by 2B_p.

This splitting is an **absolute physical energy**, already present at
lambda=0. It must be smaller than the intended lambda^2 g^2 payload
scale times its certified tolerance. It must not be multiplied by
lambda^2 when budgeting it. One may choose K last, at each fixed n,
after the nonzero couplings and finite margins have been fixed.
This gives finite K, potentially very large; it is not a moderate-ratio
or volume-uniform nonzero-port theorem. Old C84 numerical K budgets
cannot simply be copied to a changed clock/readout.

## 5. Literal two-cell coefficient and the inactive endpoint

For p=n+1 retain the scalar graph

    M, X_n,...,X_0,D_0,...,D_n.                       (12)

There are 2n+3 vertices. Only D_0,...,D_(n-1) convert. D_n and X_n
remain in the propagator, including its closing edge. Set
a=Delta-e_v(p)>0. The data diagonal is a; X,M have diagonal a+D.
Data edges have magnitude w, and the other path edges magnitude one.
After the old C84 gauge, the closing edge is -i t_d for both species.

With R=(H tensor I+I tensor H)^(-1), x=R|M,D_0>, define U_i by
retaining row D_i of x and V_i by retaining column D_i, i<n. Put

    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j,
    zeta=|v_0 v_p|^2.

The CAR calculation from C84 applies unchanged: the lambda^2 g^2
coefficient on the **whole two-cell Fock code** is scalar plus a
number-conserving quadratic operator, with

    T=-2 zeta Re C,       M=-2i zeta Im A.             (13)

Only its source visibility changed from the unbiased alpha^4 to zeta.
The field on L1 has exact compressed coefficients

    gamma_i=sum_(j=i+1)^p |v_j|^2
            =q^(2i+2)[1+O(q^2)], i=0,...,n-1.        (14)

Full literal two-cell calculations at n=1,2 verify the Gaussian identity
and this scalar mapping. The inactive endpoint is kept in those tests.
Deleting it changes even the limiting Gram matrix; an active principal
block is not the same operation as shortening the physical carrier.

## 6. Arbitrary-n positive Gram limit on this same alphabet

Set t_d=0 for this section. Let A_p(w)=aI-w adjacency(path_p), with
all p=n+1 data vertices retained, and

    f_i(t;w)=integral_0^t [exp(-A_p(t-s))]_(0i)
                          [exp(-A_p s)]_(i0) ds,
    Gamma_ij(w)=integral_0^infinity f_i(t;w)f_j(t;w) dt,
    0<=i,j<n.                                        (15)

The large-D traversal calculation gives

    2 D^(2p+2) C -> Gamma(w).                         (16)

Both express traversals finish before the low-energy propagation.
The two conversion-time orderings are the two equal triangles of
the product integral in (15), by time reversal of the real symmetric
low rail. This is precisely the applicable
[Gram lemma from the larger construction](../smaller_d_unbounded_20261116/PROOF.md),
now with the actual eleven-state carrier and its extra
inactive endpoint. The finite resolvent expansion is at each fixed n;
no fixed D is asserted to work for all n.

The small-t leading terms of the filters have different powers,
w^(2i)t^(2i+1)/(2i+1)!, so they are linearly independent and Gamma>0.
Let D_w=diag(w^(2i)). Then C(w)=D_w^(-1)Gamma(w)D_w^(-1) extends
positively to w=0, with

    C(0)_ij=(2i+2j+2)!/
      [(2i+1)!(2j+1)!(2a)^(2i+2j+3)].                (17)

The extra endpoint affects higher walk powers but not this leading
moment matrix. Its effects at finite w are retained. Any positive
finite-n Hermitian margin persists at sufficiently large finite D.
This explains why failure at one fixed express energy is not a
dimension obstruction.

## 7. An analytic all-n margin, and the exact remaining resonance

Fix n and 0<w<1/2. Set q=(wu)^2 and choose the physical readout

    h_L=W nu/q^2,
    W=lambda^2 g^2 zeta D^(-2p-2).                    (18)

After D_w congruence and removal of W, the limiting tree symbol is

    h(k)=-2 C(w) cos k+nu E(q),
    E(q)_ii=gamma_i/(q^2 w^(4i))=u^(4i)[1+O(q^2)].    (19)

Here a=Delta-e_v(p) tends to 4 when Delta=3 and q tends to zero,
so C(w) has a positive limiting matrix. The minor expansion of
det(E-r C) has a unique least-power term at each coefficient: the
minor using the first j diagonal entries of E has valuation 2j(j-1).
Positivity gives nonzero coefficients. The generalized eigenvalues,
ordered decreasingly, consequently obey

    r_j=kappa_j u^(4j)[1+o(1)], kappa_j>0.            (20)

This is the same elementary Newton-polygon step used in the larger
Gram proof, with a diagonal second pencil. Choose u so
r_(j+1)/r_j<1/32. Then choose

    0<nu<=min[1/(2r_0), sqrt(r_(n-1)/(4r_0^3))].      (21)

The 2n roots are

    k_(j,+)=acos(nu r_j/2),  k_(j,-)=-k_(j,+).

They are simple. Distinct root separation and every distinct pair-sum
comparison **except between two zero-sum opposite pairs** have the
explicit positive lower bound

    delta >= nu r_(n-1)/4.                           (22)

To see this, any nonzero integer combination of the r_j with coefficient
one-norm at most four has magnitude at least (29/32)r_(n-1), by the
largest participating term. For |x|<=1/4, |asin x-x|<=|x|^3.
The four-root remainder is at most 4(nu r_0/2)^3, which (21) bounds
by nu r_(n-1)/8. Combinations with a nonzero odd half-period remain
separated by at least pi-1. These estimates imply (22); the two-root
case is easier. Physical leading speeds are at least
(sqrt(15)/2)W lambda_min(Gamma)>0 before small finite corrections.

The excluded comparisons are not a technicality. For n>=2,

    k_(i,+)+k_(i,-)=0=k_(j,+)+k_(j,-), i!=j.          (23)

Thus the **full** pair-sum margin is exactly zero. At finite D on the
tree, T is still real and M=0; with a real readout,
h(-k)=h(k)^*. Reciprocity preserves the same obstruction. No choice
of the geometric readout removes it. Simplicity and the non-pair
margin of the quadratic symbol do persist at sufficiently large finite D
and sufficiently small quadratic corrections, at each fixed n. Higher
orders also generate interactions, so this continuity statement is not
a theorem about the interacting many-body spectrum. It does not establish
the missing Cooper-transfer nonresonance or replace interacting RG analysis.

## 8. Why the complete larger-alphabet proof does not transfer for free

The larger all-n proof uses two **independent bus labels**, each with
both conversion charges, different data speeds, and independently
phased ports. Its leading hopping is a complex combination of two
independent positive Gram matrices. C84's A/B labels instead are the
two charges converted into each other. They are not independent buses.
Renaming them cannot reproduce that block-orthogonality argument.

The faithful duplicated-bus shifted alphabet would have four backgrounds,
eight charged data/express letters, and five markers: seventeen physical
states. This is a state-count audit of that transfer, not a proposal to
replace the eleven-state benchmark by a larger construction. More general
contextual implementations are not excluded by this count.

C84's genuine loop can evade the tree reciprocity in (23), and its
existing finite certificates demonstrate that it sometimes does. What
is missing is an all-n nonidentity or margin theorem for that **single
loop with the inactive endpoint**, not the ability to encode arbitrary
n, nor the rank of the limiting Gram matrix. The failed exponent-66
bracket prescription and fixed-D sign failures do not supply such an
obstruction. Conversely, one n=11 certificate does not prove it.

The current result therefore transfers the clock, selector, readout,
Gram positivity and simple-band ingredients at d=11, and exhibits an
exact obstruction to transferring the **tree nonresonance conclusion**.
It neither rules out an all-n eleven-state loop construction nor the
reported ten-state construction, whose source remains unlocated.

## 9. Physical conclusions still conditional

For the original certified C84 member, nonzero-port sector selection,
volume-uniform auxiliary elimination, interacting multichannel stability
and conformal-observable identification remain required. The modified
family additionally needs an all-n full pair-sum margin, or a controlled
resonant interacting alternative. The old n=11 numerical certificate
does not certify these modified couplings verbatim.

No mobile lift is constructed in this round. Two colours suffice only
under the established pinless grammar hypotheses; the twelve-state
one-vacancy route has different slow dynamics and extra interactions.
Adding one hole to the nine-state seed is not automatically a d=10
mobile model. The complementary chat studies that smaller seed.

No claim of unbounded physical c, or a new microscopic c at any n, is
made. Any eventual c=n or c=n+1 conclusion must still identify a unitary
compact infrared theory and its total nonchiral Virasoro charge.
