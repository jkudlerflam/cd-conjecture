# Tail-sensitive mobile length comparisons

2026-11-17. Author and AI assistance: OpenAI Codex, GPT-6 family.
**[Proved below, awaiting independent review.]** The results are a
nonperturbative scalar reduction and a favorable second-order clock
comparison. They are **not** a uniform mobile length-selection theorem.
The precise remaining scalar-response and higher-order questions are in
§8. Numerical evidence is separated in README.md and checks.json.

## 1. Literal sector, flags, and ensemble

Use the parent's unchanged 22-state nearest-neighbor mobile bond, with
lambda=f=0. Fix the number N of occupied material sites and a cyclic
ordered word of cell lengths k_c>=1, with sum_c(k_c+1)=N>=3. Fix all
forbidden-arrow flags to zero and every cell's charged/vacuum flag to
vacuum. Fix a payload bit string; all-zero bits are sufficient since
all payload strings give identical vacuum clock matrices at f=0.
Marker positions are fixed in **ordered material coordinates**, not at
fixed physical sites. The vacancies and the clocks remain quantum and
can be entangled. No classical head-position pin is imposed. Vacancy
colors are uniquely determined by the legal arrows.

The charged buses, Delta,D,g,w,v,eta and the finite grammar penalty K
have no matrix elements in this invariant block. Setting lambda=f=0
does not prove that the full model's ground state lies in this block.
The earlier charged-rail theorem remains a separate pin result with
its own allowance; it is not implicitly upgraded to J/t=1/80 here.

Let g_j>=0 be the vacancies after material site j, A_j=1_(g_j=0), and
S_j the unilateral gap-raising shift. For a cell starting at a and
ending at marker m=a+k, the **literal** Hamiltonian is

    H_w = H_ch + J sum_c [sum_(r=0)^(k_c-1) A_(a_c+r) E_(c,r)
                          + beta C_(m_c) T_c + mu C_(m_c)],
    H_ch = -t sum_j(S_(j-1) S_j^* + h.c.),
    E_r = |q r-(r+1)><q r-(r+1)|,
    beta=q^2, T_c=|k_c><k_c|,
    C_m=(A_(m-1)+A_m)/2.                              (1)

In E_r the ket notation means q|r>-|r+1>. Both contact halves of
the terminal and marker onsite terms are present. Only -mu_occ N is
a scalar at fixed N and has been removed. In particular mu C_m is
not a constant, even at fixed marker count. checks.py independently
applies the literal 400x400 material bond to both tested length words.

We use the positive-length tilt H_w+sV, V=sum g_j, with
t/2<=s<=3t/4. Write F_w(s)=min_V[E_w(N,V)+sV]. This is a ground
energy, not a thermal free energy. A comparison proved for F_B-F_A
also holds at any V_B minimizing the **competitor B**:

    E_B(N,V_B)-E_A(N,V_B) >= F_B(s)-F_A(s).             (2)

The direction matters: this does not automatically prove the same
bound at an A-minimizing V or an arbitrary prescribed density. When
both words minimize at the same V, their differences agree. The n=1
checks certify that situation at s/t=0.6 by excluding every other V
with free-sea lower bounds. Some n=2 checks use other canonical blocks
only as diagnostics; they are not passed off as tilted minima.

## 2. The weak first-order comparison

Let Q_k=sum_r E_r+beta T and E_k its lowest eigenvalue. The established
frozen bounds are

    (1-beta)^2 q^(2k+2) <= E_k <= q^(2k+2).             (3)

Set mu=-(n+2)E_n. Compare

    A: n+2 cells of n bits,
    B: n+1 cells of n+1 bits,
    N=(n+1)(n+2).

At fixed N,V the free gap ground state Omega is unique and has common
contact expectation

    alpha = N/(N+V)
      - sin^2(pi N/(N+V))/(N(N+V)sin^2(pi/(N+V))).      (4)

Compression onto Omega is exactly E_ch+J alpha sum_c(Q_k+mu).
Thus the right derivative at J=0 gives

    E_B-E_A = J alpha Delta_n + O_(N,V)(J^2),
    Delta_n=(n+1)E_(n+1).                             (5)

Indeed (n+1)(mu+E_(n+1))-(n+2)(mu+E_n)=Delta_n.
The common charge energy and occupied chemical potential cancel
exactly. The marker-count change is essential. Repeating both words
r times multiplies the first-order margin by r at the same N,V.
For a unique free tilted minimizing V, (5) also gives the tilted
derivative. At a tilt crossing, minimization among tied blocks must
be kept explicitly; no differentiability is assumed there.

## 3. A nonperturbative bound retaining the terminal tail

For each cell define its frustration-free vector and tail by

    u_k=(1,q,...,q^k)/sqrt(Z_k), Z_k=sum_(j=0)^k q^(2j),
    P_c=|u_k><u_k|, p_k^0=q^(2k)/Z_k,
    tau_k=beta p_k^0, d_k=1+beta-2q cos(pi/(k+1)).      (6)

Remove only beta C_m T from (1), obtaining H_0. Every edge annihilates
u_k; every P_c commutes with H_0. For scalar coefficients v_c define

    F_w^sc(v;s)=inf spec[H_ch+sV+J sum_c v_c C_(m_c)].   (7)

The exact ground energy of H_0+sV is F_w^sc(mu;s), with product clock
u_k and the fully interacting scalar vacancy ground state. This
reference retains the entire marker potential and the moving fluid.
No contact expectation or bath gap enters this assertion.

Assume

    0<q<=1/4, n>=1, mu=-(n+2)E_n,
    0<J/t<=1/80, t/2<=s<=3t/4.                        (8)

The coherent-gradient pin proof of the previous contribution, with
the terminal term absent, gives delta_k=d_k/257. Its source bound
is c_k+|mu|/2, with c_k<=25/16 and |mu|/2<=3/512;
the same allowance (8) suffices for every k>=1. It works with any
other cells simultaneously pinned, since their terms and projectors
are unchanged by the target-cell channel. Successively replacing
pinned clocks by u_k lowers the ground energy by at least J delta_k
per pin. Since all P_c commute with H_0, block decomposition proves
the operator inequality

    H_0+sV-F_w^sc(mu;s) >= J sum_c delta_(k_c)(1-P_c).  (9)

This is coercivity against FF-clock excitation in the terminal-free
reference. It is not a charge gap or a spectral gap of the full model.

For 0<theta<1 put

    kappa_k=theta delta_k/(theta delta_k+beta).

For any 0<=C<=1 commuting with the clock, the two-dimensional
determinant/trace bound for the nontrivial block gives

    theta delta_k(1-P_c)+beta C T_c >= kappa_k tau_k C.

The remaining clock eigenvalues equal theta delta_k and obey the same
inequality. Split off theta of (9), use this inequality cell by cell,
and observe that the remaining FF edges still have the same exact
kernel. One obtains the **all-order, tail-sensitive scalar sandwich**

    theta F_w^sc(mu;s)
      +(1-theta) F_w^sc(mu+kappa_k tau_k/(1-theta);s)
        <= F_w(s) <= F_w^sc(mu+tau_k;s).               (10)

The upper bound is compression onto the product u_k. Every added
potential in (10) retains q^(2k+2). There is no absolute J^2 N/t error,
no charge cut, and no charged-mass factor or growing k in J/t.
Both sides use the original gap fluid with the actual marker pattern.

A better upper bound uses the frozen eigenvector v_k of Q_k. Define
local positive sources b_j=E_j for j<k-1,
b_(k-1)=E_(k-1)+beta T/2, b_k=beta T/2. Let

    a_(k,j)=<v_k,b_j v_k>,
    sigma_(a+j)=a_(k,j)+(mu/2)1_(j=k-1 or k).

Then

    F_w <= U_w := inf spec[H_ch+sV+J sum_j sigma_j A_j]. (11)

Consequently L_B-U_A>0, where L_B is (10)'s lower side, is an explicit
sufficient condition for the weak comparison. These are **scalar**
vacancy problems on the same N sites, not independent extensive-error
estimates. No new Hamiltonian term or retuning is required to formulate
this condition. Its uniform positivity is not proved here. At n=1,
N=V=6,J/t=1/80, numerical evaluation gives L_B-U_A about 4.08e-8>0.
This is a finite numerical evaluation of a rigorous sufficient bound,
not an interval-arithmetic certificate or a thermodynamic conclusion.

For completeness, if P is the product of the u_k projectors, the exact
off-diagonal block has orthogonal single-cell excitation ranges:

    ||(1-P)H_w P||^2 <= J^2 beta sum_c tau_(k_c).       (12)

This also retains the tail. It cannot alone be divided by a single-cell
pin cost to obtain a full interacting product-code gap.

## 4. Exact finite-volume second-order decomposition

At fixed finite N,V define the free response matrix

    x_j=(A_j-alpha)Omega,
    G_ij=<x_i,(H_ch-E_ch)^(-1)_perp x_j>.               (13)

G is real, positive semidefinite, and circulant. This is the exact
finite-volume coefficient; **no estimate using the reciprocal charge
gap is claimed uniform in N**. Its eigenvalues gamma(theta)>=0 at
theta=2pi l/N will only be used for a sign comparison below.

Use frozen v_k, not u_k, in perturbation theory. Set

    xi_(k,j)=(b_j-a_(k,j))v_k,
    F_k(theta)=||sum_(j=0)^k exp(i j theta) xi_(k,j)||^2,
    R_k=(1/N)sum_theta gamma(theta) F_k(theta).         (14)

Because sum_j b_j=Q_k, sum_j xi_(k,j)=0 exactly. Centered responses
on distinct cells are orthogonal. Degenerate perturbation theory in
the free-sea subspace, whose first-order ground vector is the product
v_k, therefore gives

    E_w = E_ch + J alpha sum_c(mu+E_(k_c))
          - J^2 [sigma^T G sigma + sum_c R_(k_c)]
          + O_(N,V,w)(J^3).                           (15)

The scalar part sees the positions of all markers. The clock part is
additive in the length counts, independent of their arrangement at this
order. The constants in the remainder are not uniform; (15) is not an
all-order thermodynamic matching theorem.

## 5. An all-length favorable sign for the clock correction

For every k>=1, 0<q<=1/4, and real theta,

    F_(k+1)(theta) <= (512/81)q^2 F_k(theta)
                         <= (32/81) F_k(theta).        (16)

Here is a direct proof, including the endpoint estimates needed for
uniformity in k. Write E=E_k, v=v_k>0, p=v_k^2, beta=q^2, and
g_j=q v_j-v_(j+1). The eigenvalue equation gives

    g_j=E q^(-j-1) sum_(l=0)^j q^l v_l > 0,
    g_(k-1)=(beta-E)v_k,
    g_j=q g_(j+1)-E v_(j+1) <= q g_(j+1).             (17)

In particular g_(k-1-r)<=q^r beta v_k. Also

    (25/32)q^(2k) <= p <= q^(2k).                     (18)

The upper bound follows from beta p<=E and (3). To check the lower
bound let Q_ff=sum E_j and d>=9/16 its nonzero gap. From
Q_ff v=Ev-beta v_k|k> one gets

    (d-E)<Q_ff> <= beta(beta-E)p <= beta^2 p.

Thus p>=(E/beta)(d-E)/(d-E+beta). Since E<=1/256 and beta<=1/16,
d-E>=8 beta, so (3) yields p>=(8/9)(15/16)^2 q^(2k), which is (18).

Put z=exp(i theta), a(z)=sum a_(k,j) z^j. Positivity of the b_j
gives a_(k,j)>=0 and sum a_(k,j)=E. The last component of the
centered vector in (14) is

    v_k [z^(k-1)(E+beta(z-1)/2)-a(z)].

Since |a(z)-E z^(k-1)|<=kE|z-1| and
kE<=beta q^2<=beta/16, this proves

    F_k(theta)>=(7/16)^2 beta^2 p |z-1|^2.             (19)

For the upper bound subtract E z^(k-1)v before orthogonal projection:

    [sum z^j b_j-E z^(k-1)]v
      =sum_(j=0)^(k-2)(z^j-z^(k-1))g_j(q|j>-|j+1>)
         +(beta/2)(z^k-z^(k-1))v_k|k>.

Using (17), |z^j-z^(k-1)|<=(k-1-j)|z-1|, and
sum_(r>=1)r q^r=q/(1-q)^2 gives

    F_k(theta) <= beta^2 p |z-1|^2
                 [1/2+sqrt(1+beta)q/(1-q)^2]^2
               <= (35/36)^2 beta^2 p |z-1|^2.         (20)

Divide (20) at k+1 by (19) at k and use (18).
The coefficient is (400/81)(32/25)q^2=(512/81)q^2.
At z=1 both sides vanish, so (16) holds there too.

Since all gamma(theta)>=0, (16) implies R_(k+1)<=r R_k,
r=(512/81)q^2<=32/81, without any bound on the bath susceptibility.
For the weak pair, the **clock part** of the J^2 energy difference is

    D_clock=(n+2)R_n-(n+1)R_(n+1)
        >= [n+2-r(n+1)]R_n >= 0.                     (21)

The same statement holds for any common repetition count and any
finite N,V. Thus the larger O(E_n) clock correction has a favorable
sign, rather than being an arbitrary absolute error competing with
the smaller margin E_(n+1). This proves a coefficient comparison;
it does not make the perturbation series uniform in volume.

## 6. The marker response really does not cancel

The other second-order term is

    D_scalar=sigma_A^T G sigma_A-sigma_B^T G sigma_B.   (22)

Write sigma_A=-E_n 1+r_A and
sigma_B=(-E_n+Delta_n/N)1+r_B, with sum r_A=sum r_B=0.
Circulant symmetry gives the exact cancellation of the common
uniform bulk response:

    D_scalar=gamma(0)[2E_n Delta_n-Delta_n^2/N]
                + r_A^T G r_A-r_B^T G r_B.            (23)

The uniform contribution is favorable. The remaining periodic marker
potentials have different nonzero Fourier modes: periods n+1 and n+2.
There is no identity equating their responses. In the actual free bath
at n=1,N=V=6, checks give

    D_scalar = -3.25703265e-6/t,
    D_clock  = +2.65077799e-5/t,
    D_total  = +2.32507473e-5/t.                       (24)

This is a concrete failure of a proposed exact cancellation, not a
failure of the selector. The favorable clock term wins in these data.

The scalar source is nevertheless smaller in tail order. Each cell's
positive a_j sum to E_k, and its two marker halves are mu/2. Hence

    ||sigma||^2 <= sum_c(mu^2+2E_(k_c)^2).              (25)

For the tested words this is polynomial in n times E_n^2, whereas
the clock response is of order beta E_n. For example the same elementary
Q_ff estimate gives sum_j||xi_(k,j)||^2<=5 beta E_k.
These power counts are exact coefficient bounds before multiplying by
G. A uniform estimate on the relevant spatial response, or an
all-order comparison of (10)-(11), is still needed. Replacing G by
the reciprocal finite-chain charge gap loses uniformity and is not
a resolution. An uncontrolled thermodynamic Taylor expansion is not
justified by the extra E_n factor either.

## 7. Chemical-potential crossing and physical gluing

The frozen weak-pair crossing is

    mu_c(0)=(n+1)E_(n+1)-(n+2)E_n=mu_sel+Delta_n.

At fixed finite N,V the implicit-function expansion from (15) is

    mu_c(J)=mu_c(0)+J D_total(mu_c(0))/alpha+O_(N,V)(J^2). (26)

The full source, including mu, must be reevaluated at mu_c(0) for
this formula. At n=1,N=V=6,J/t=1/80, the observed crossing moves
upward by about 1.0006e-6, making the chosen mu_sel farther from
this competitor. No retuning is required in that finite check. This
does not determine the all-n thermodynamic crossing or certify a
uniform interval of mu for every wrong length.

There is an extra simplification in this restricted vacuum problem.
In its literal basis all off-diagonal entries are real nonpositive:
-t for vacancy motion and -Jq for clock motion. A zero-flux rooted
ground vector can be chosen nonnegative. Sum its rotated copies over
the cyclic orbit of the chosen material word to obtain a rotation-
invariant vector. It satisfies the physical gluing condition at total
momentum P=0 and gives a physical state with the same energy. Conversely
the diamagnetic inequality bounds every twisted rooted block below by
the untwisted rooted ground energy. Therefore, for this sector,

    E_phys(N,V)=E_rooted(0;N,V)                        (27)

exactly, for every N,V, and also after minimizing the positive tilt.
The same reasoning applies to the scalar comparison Hamiltonians.
No O(t) cut or extensive boundary error occurs. If one instead uses
the earlier general gluing theorem, there is a single comparison
allowance 96t/N; (27) improves it to zero here. Once the complex
charged/port processes are restored, this special argument is not
available, and the general global remainder must be reinstated.

## 8. What is proved, and the precise obstruction left

The answer is stronger than an absolute perturbation estimate: (10)
is a nonperturbative tail-sensitive reduction, and (21) proves the
potentially larger second-order clock correction favors the intended
word for every n. The scalar marker response is smaller in microscopic
tail order, but actually nonzero and unfavorable in the tested baths.

There is still no uniform proof of L_B-U_A>0 at J/t=1/80, nor a
uniform all-order replacement for (15). Thus fixed mu_sel survives
the bounded tests; its thermodynamic validity is **unresolved**, not
refuted and not established. No retuning or shrinking J/t is prescribed
merely to cover an unproved remainder.

A sufficient next theorem would establish the scalar inequality
L_B-U_A>=c J alpha Delta_n with c>0 uniform in common repetition
count, in the positive-tilt ensemble, for the explicit source profiles
in §3. Alternatively a relative all-order influence estimate could
retain the favorable clock sign and control (23). Both targets leave
the charge fluid intact and avoid a charge-gap denominator.

Even that would only solve this particular marker-changing competitor.
Arbitrary words require all k (including zero-bit cells and marker-free
words), arrangements, mixtures and isolated defects, with a common
interacting bulk comparison. The earlier charged pin theorem requires
its stated extra bounds and is not uniform for arbitrarily long pinned
cells. Nonzero ports/readout require relative errors below the weak
selection margin, not only small local occupation. Illegal sectors
at finite K still need a separate energy comparison. None of these
results establishes the proposed unitary compact multichannel CFT
or its total nonchiral Virasoro central charge c=n+1.
