# Length fugacity closes the reference ensemble gap

Status: **[proved abstract pressure and thermal comparison theorems]**.
Application to the positively completed C44/C45 reference uses its recorded
local reduction and the root-to-physical comparison in the sibling
`ensemble/` folder. The full-space completion is independently audited in
`COMPLETION_AUDIT.md`. No theorem about the leakage remainder,
interacting Gaussian limit, or central charge is asserted.

Provenance: OpenAI GPT-6-family Codex agent `exact_escape`, 2026-10-28.
AI-assisted derivation; not human-refereed. The task was redirected from
exact encodings to the ensemble gap to avoid duplicating another chat.
The proof extends C46's thermal theorem to a physical-length fugacity and
uses the resulting pressure zero, rather than a fixed-density free-gas
Legendre calculation. It applies to interacting length-dependent bonds
and to arbitrarily long excluded cells.

## 1. Definitions and hypotheses

Let a rooted cell have Hilbert space

    K = direct_sum_(ell>=ell_min) K_ell,  ell_min>=1,
    L_cell|K_ell = ell I.

Every K_ell is finite dimensional. On M rooted cells write

    L_M=sum_j L_cell,j,
    H_M=sum_j A_j+sum_X B_X,
    [H_M,L_M]=0.

Use open-cell Hamiltonians H_M in (2) and the cutting argument below.
The A_j are translates of a semibounded on-cell operator A commuting with
L_cell. The B_X are bounded, translation-covariant interactions. Require
a finite total interaction strength per cell and a finite cut strength

    J=sup_j sum_(X contains j)||B_X|| < infinity,
    C_cut=sup_cut sum_(X crosses cut)||B_X|| < infinity.       (1)

An exponentially weighted interaction norm suffices for both. Periodic
conventions differ from the open ones by an additional bounded seam and
have the same pressure; their finite-size additivity constant may be larger. Very
small M can be specified separately without affecting any limit below.

Let beta>0, T=1/beta, and let mu be a marker chemical potential. Define

    Z_(M,L)=Tr_(L_M=L) exp[-beta(H_M-mu M)],
    Z_M(s)=Tr exp[-beta(H_M+s L_M-mu M)]
          =sum_L exp(-beta s L) Z_(M,L),
    p_cell(T,s)=lim_(M->infinity) (T/M)log Z_M(s).             (2)

The sign convention is that positive s penalizes physical length.
Assume an open interval I on which

    Tr exp[-beta(A+s L_cell)] < infinity.                    (3)

The physical canonical partition functions need not agree with any
single rooted fixed-M sector. It suffices that their sum satisfies

    c_beta sum_(M>=1) Z_(M,L) <= Z_phys,L
       <= C_beta sum_(M>=1) Z_(M,L),                         (4)

for all sufficiently large L, with 0<c_beta<=C_beta<infinity independent
of L. The sibling pin-and-cut argument is designed to supply (4).
Polynomial factors in L would give the same result, by shifting s by an
arbitrarily small positive amount; constants are used here for clarity.

The physical finite-dimensional nearest-neighbor chain has a pressure

    p_phys(T)=lim_(L->infinity) (T/L)log Z_phys,L.             (5)

One can prove existence directly by deleting the single bond joining
two open intervals: the logarithms of their partition functions are
additive up to beta times that bond's norm. Closing a ring changes the
logarithm by another bounded amount. Thus existence in (5) does not
presuppose a CFT or a positive spectral gap.

## 2. The pressure-zero theorem

**[proved]** Suppose (1)-(4) hold and p_cell takes both signs inside I.
There is a unique s_*(T) in I with p_cell(T,s_*)=0, and

    p_phys(T)=s_*(T).                                      (6)

Proof of existence of the rooted pressure. Boundedness of the total
interaction gives finite Z_M(s), using the product on-cell trace (3).
Deleting the terms crossing a cut gives

    |log Z_(M+N)(s)-log Z_M(s)-log Z_N(s)|<=beta C_cut.       (7)

This follows from the eigenvalue min-max principle applied to a bounded
perturbation, and remains valid for trace-class semibounded operators.
For fixed k, write M=qk+r and iterate (7). Dividing by M bounds the
limsup and liminf by (T/k)log Z_k(s) plus or minus C_cut/k.
The bounded remainder r<k disappears as M tends to infinity. Letting
k tend to infinity proves the limit in (2).

The pressure is convex and finite on I. For finite M this follows from
the positive sum over L in (2); convexity passes to the limit. A finite
convex function on an open interval is continuous there. Moreover, if
s'>s, L_M>=ell_min M gives

    p_cell(T,s')<=p_cell(T,s)-ell_min(s'-s).                 (8)

Thus the sign change gives a unique interior zero.

Now put S_L=sum_M Z_(M,L). Positivity and Tonelli's theorem give the
exact generating-series identity

    sum_(L>=1) exp(-beta s L) S_L = sum_(M>=1) Z_M(s).       (9)

If p_cell(T,s)<0, the summands on the right eventually decay
exponentially in M, so the series converges. If p_cell(T,s)>0, its
summands grow exponentially and it diverges. Behavior exactly at zero
is immaterial to the convergence abscissa. Equations (8)-(9) therefore
identify that abscissa as s_*.

By (4), the series using Z_phys,L has the same abscissa. Finally, (5)
implies convergence for s>p_phys and divergence for s<p_phys: for any
epsilon>0, the large-L coefficients lie between
exp[beta(p_phys-epsilon)L] and exp[beta(p_phys+epsilon)L].
This proves (6).

This proof uses neither a free length band nor a product payload/length
Hamiltonian. It retains every allowed M and every cell length. The
interior sign change is essential: if the zero sits at the boundary
of the trace-class interval, a macroscopic long-cell or marker-free
phase can dominate. Such a boundary case is not covered by the theorem.

### 2.1 An actual selected renewal pressure limit

**[proved]** Even without a physical-space identification (4), the
selected code with lengths b and b+1 has an actual fixed-L renewal
pressure, not merely an abscissa or a limsup:

    p_renewal,P(T)=lim_(L->infinity) (T/L)log S_P,L,
    S_P,L=sum_M Z_P,(M,L).                                (9a)

It equals the unique selected rooted-pressure zero whenever that zero
is interior. This allows (17) below to be read as a comparison to a
well-defined thermodynamic reference without first proving a separate
physical pinning theorem for a constrained code.

To prove existence, fix positive L1,L2. In a K-cell chain, let E_M project
onto states whose first M cells have total length L1 and whose remaining
K-M cells have total length L2. For different M these subspaces are
orthogonal: every cell has strictly positive length, so a length-basis
word cannot have two different prefixes of length L1. This avoids an
otherwise spurious multiplicity factor from summing cell counts.

Compress the full H_K to E_M. Its difference from
H_M restricted to L1 plus H_(K-M) restricted to L2 has norm at most
C_cut. If each local B_X conserves total length on X, the internal
terms preserve the two subspaces and this assertion follows directly.
This local conservation is satisfied by the rooted length-preserving
matrix-unit rules; we include it as a hypothesis of this subsection.

For mutually orthogonal projections E_M, scalar spectral convexity in
an eigenbasis of each compressed block gives the Peierls inequality

    Tr exp(-beta H_K) >= sum_M Tr_(E_M) exp(-beta E_M H_K E_M).

Apply it inside total length L1+L2, and compare each bounded cut term
by the min-max principle. Summing K yields

    S_P,(L1+L2) >= exp(-beta C_cut) S_P,L1 S_P,L2.          (9b)

The marker term is additive and does not alter this bound. Every
sufficiently large L is a sum of b's and b+1's, so S_P,L>0 eventually.
For example, writing L=qb+r with 0<=r<b, replace r of the q b-cells by
b+1-cells whenever q>=r, which holds for L>=b(b-1).
Finite d_P and a bounded Hamiltonian norm per cell give an exponential
upper bound on S_P,L.

Set a_L=log S_P,L-beta C_cut. Equation (9b) says a_(L1+L2)>=a_L1+a_L2.
For any sufficiently large fixed k, write large L=qk+r with r in a fixed
interval of k consecutive attainable lengths. Iteration gives
liminf_L a_L/L>=a_k/k. Taking the supremum over sufficiently large k
and comparing with the limsup proves convergence. The same generating
series (9) then identifies its value with the pressure zero. No free
length spectrum or density differentiability was used.

### 2.2 A thermal coefficient needs only one ground-pressure derivative

**[proved]** Suppose, uniformly on a fixed neighborhood of s0,

    p_cell(T,s)=p0(s)+A(s)T^2+r_T(s),
    sup_s |r_T(s)|/T^2 ->0,
    p0(s0)=0,       p0'(s0)=-ell0<0,                     (9c)

where A is continuous at s0. Suppose the finite-temperature pressures
obey (8), and their interior zeros s_T lie in that neighborhood. Then

    s_T=s0+[A(s0)/ell0]T^2+o(T^2).                       (9d)

Indeed p_cell(T,s0)=A(s0)T^2+o(T^2). Inequality (8), applied on the
appropriate side of the zero, gives
|s_T-s0|<=|p_cell(T,s0)|/ell_min=O(T^2).
Differentiability at the single point s0 now gives
p0(s_T)=-ell0(s_T-s0)+o(|s_T-s0|); continuity gives
A(s_T)=A(s0)+o(1). Substitution in p_cell(T,s_T)=0 proves (9d).
Uniform convergence and an interior sign bracket can replace the stated
zero-location condition, since (8) then forces the zeros into any fixed
neighborhood of s0 for small enough T.

No density Hessian, C^2 ground energy, or temperature derivative of the
error is needed. If A(s0) is the established cell-coordinate Gaussian
coefficient, dividing it by the mean ground-state cell length ell0 is
the physical pressure coefficient. This is consistent with multiplying
cell-coordinate propagation speeds by ell0. The lemma supplies the
ensemble conversion, not the assumed uniform scalar expansion in (9c).

## 3. Length-tilted thermal comparison

Now specialize to a locally block-diagonal reference D_M=A+B. Let P be
a finite-dimensional on-cell code of rank d_P, Q=1-P, with

    AP=0, [P,L_cell]=0, [B_X,P_X]=0,
    P_X=product_(j in X)P_j.

Write H_P=P_all D_M P_all and L_P=P_all L_M P_all. Put

    ell_bar = Tr_P(L_cell)/d_P,
    t_s(beta)=Tr_Q exp{-beta[A+s(L_cell-ell_bar)-2J]}.       (10)

Assume this trace is finite. It need not be bounded by the dimension
of Q, which may be infinite.

**[proved]** For every M and beta,

    1 <= Z_D,M(s)/Z_P,M(s)
      <= [1+t_s(beta)/d_P]^M.                            (11)

Both traces include the same marker term -mu M. Consequently

    0<=p_D(T,s)-p_P(T,s)
      <=delta_T(s):=T log[1+t_s(beta)/d_P].                (12)

Proof. Use the product unital positive extension E from C46: on each
cell, E(O) is O on P and (Tr_P O/d_P)I on Q. In a block with Q exactly
on a set S, E is the normalized partial trace over the erased code
factors. Its positivity and contraction in operator norm follow from
that explicit form, not from tensoring an arbitrary positive map.
Local block diagonality gives

    D_M-E(H_P)>=sum_j(A_j-2J Q_j).

The length term obeys an exact identity,

    L_M-E(L_P)=sum_j Q_j(L_cell,j-ell_bar).                (13)

Therefore

    D_M+s L_M-mu M
      >=E(H_P+s L_P-mu M)
         +sum_j Q_j[A_j+s(L_cell,j-ell_bar)-2J].           (14)

This is the main gain over treating the tilt as an additional bounded
interaction. No spurious penalty proportional to |s|b is added to J.
The unbounded tail in s L_cell is retained inside the on-cell heat
trace where the linear selector cost controls it.

For completeness, the needed partial-trace inequality is

    Tr exp[-beta d_S^(-1)Tr_S K] <= (Tr exp[-beta K])/d_S, (15)

where K is any finite code Hamiltonian and d_S=d_P^|S|.
Choose an eigenbasis v_i of d_S^(-1)Tr_S K. Scalar convexity gives
exp[-beta<v_i,d_S^(-1)Tr_S K v_i>]
<=d_S^(-1)sum_a exp[-beta<v_i,a|K|v_i,a>].
Spectral convexity bounds each final scalar exponential by the
corresponding diagonal matrix element of exp(-beta K). Summing i proves
(15). Apply (15) to K=H_P+s L_P-mu M in each Q-set block of the right
side of (14). The independent auxiliary traces contribute t_s^|S|.
The sum over S is the binomial product in (11). The lower bound follows
because the full code is an invariant summand of D_M+sL_M. This proves
(11)-(12), including noncommuting defect motion within local Q_X blocks.

## 4. Fixed-physical-length exclusion of the reference auxiliaries

Assume the code uses only lengths b and b+1, b>=1. Suppose both D and
H_P satisfy the pressure-zero theorem's physical comparison, and their
rooted pressures are defined on an interval containing

    [s_P(T), s_P(T)+delta_*(T)/b],
    delta_*(T)>=sup_s delta_T(s)                           (16)

over that interval. Then

    0<=p_phys,D(T)-p_phys,P(T)<=delta_*(T)/b.               (17)

Indeed p_D(s_P)>=p_P(s_P)=0. The selected code has L_P>=bM, so
p_P(s_P+x)<=-b x. Equation (12) makes p_D nonpositive at
x=delta_*/b. Its strictly decreasing continuous pressure has a unique
zero between these endpoints. Equation (6) identifies those zeros
with the physical pressures, proving (17).

In particular, if t_s/d_P<=C exp(-beta delta) uniformly near the selected
root, with delta>0, then

    0<=p_phys,D-p_phys,P <= (C/b)T exp(-delta/T).            (18)

This controls the actual fixed-physical-length pressure after summing
marker numbers. It repairs the specific obstacle that the partial trace
in C46 need not preserve a fixed-L subspace. The partial trace is used
only after length fugacity has removed that constraint; (6) then restores
the physical ensemble.

If either physical pressure has an existing low-temperature expansion,
(18) gives the same ground pressure and the same coefficient of every
fixed power of T for the other. It does not create a T^2 expansion when
none has been proved. A central charge still requires a Gaussian or
stress-tensor identification and the relevant propagation speeds.

## 5. Explicit trace window from the C44 selector bounds

**[proved conditional on the stated imported inequalities]** Use d0=4464,
b>=2, code rank d_P=2^(n+1), and equal code ranks at b and b+1, so
ell_bar=b+1/2. The following are the heat-trace hypotheses used in C46:

* At the exceptional lengths {1,b,b+1}, excluded states have A>=Delta.
* At other lengths ell>=2, A>=gamma(1+|ell-b|), gamma>=Delta>0.
* dim K_ell<=d0^ell, and J<=Delta/48.

For the full positively completed bond, take gamma to mean gamma_all
from `COMPLETION_AUDIT.md` (5), and shrink Delta to at most that number
and the complete selected-cell gaps. These hypotheses then cover all
strings. The original selector-legal gamma must not be reused for illegal
strings without that adjustment. The sibling's simpler kappa bound is
an alternative way to get the same trace-class conclusion.

Fix s0>0 such that

    s0(b-1/2)<=Delta/4,       s0<=gamma/4.                 (19)

For |s|<=s0 the exceptional excluded states have tilted energy at least
Delta_s=Delta-s0(b-1/2)>=3Delta/4. On all other lengths,
|ell-ell_bar|<=1+|ell-b|, so their tilted cost is at least
gamma_s(1+|ell-b|), gamma_s=gamma-s0>=3gamma/4>=Delta_s.
The last inequality follows from gamma>=Delta and b-1/2>=1.

Set

    beta0=log(2d0)/gamma_s,
    C_th=[d0+d0^b+d0^(b+1)+b d0^(b-1)+2d0^(b+2)]/d_P.    (20)

The same elementary geometric sum as C46, now with gamma_s, gives

    t_s(beta)/d_P <= C_th exp[-beta(Delta_s-2J)]
                   <=C_th exp(-beta Delta/2),             (21)

for beta>=beta0, uniformly over |s|<=s0. The finite exceptional lengths
give the first three terms in (20). Short wrong lengths give at most
b d0^(b-1). The long tail, ell=b+2+r, is bounded by
d0^(b+2)exp(-3beta gamma_s)/(1-d0 exp(-beta gamma_s)),
whose denominator is at least 1/2. Multiplication by exp(2beta J)
and the gap bounds give (21).

Thus any selected physical-pressure root staying a fixed distance inside
(-s0,s0) has an exponentially accurate full-reference counterpart via
(17). This root-location condition is checkable rather than automatic.
For example, if ||H_P||<=KM and |mu|+K is sufficiently small compared
with b s0, the elementary code dimension bound

    p_P(T,s0)<=|mu|+K+T log d_P-b s0,
    p_P(T,-s0)>=- |mu|-K+b s0                              (22)

brackets a root for low enough T, with a margin. In the intended
application the large selector tangent e L must first be subtracted;
otherwise the relevant fugacity is centered at -e, not at zero.
Subtracting e L increases the pressure by e, so the original pressure
is the centered pressure minus e. All scalar port shifts and the marker chemical
potential belong in K and mu when checking (22).

## 6. Leakage remains a separate problem

If H=D+R has total rooted norm at most M epsilon, then its tilted pressure
differs from p_D(T,s) by at most epsilon. The pressure-zero theorem gives
an O(epsilon) physical-pressure comparison, using the positive minimum
cell length. This error is independent of T and is not o(T^2).

Therefore (18) closes the reference's all-length physical-ensemble
transfer, conditional on the pin-and-cut identification. It does not
upgrade C45's finite-order local remainder to a microscopic central-charge
theorem. The missing low-temperature leakage and interacting retained
pressure estimates are exactly the obligations left in the other folders.

The original C45 application had an additional full-space boundary:
conditional endpoint pins were not on-cell on arbitrary malformed words.
The positive completion and token-count split in the sibling
`ensemble/ILLEGAL_SECTORS.md`, independently checked in
`COMPLETION_AUDIT.md`, repair that structural problem for the modified
bond. They provide the complete on-cell gap, linear heat-trace cost, and
bounded residual interactions without discarding malformed sectors.

The finite-rank projected code is
a constrained physical subspace. A full-local-space pinning dimension
ratio does not automatically apply to its trace. Without a separate
projected-code seam comparison, s_P is a rigorously defined *selected
renewal reference pressure*, and (17) compares to that reference; it
must not be silently called the projected physical pressure.

## 7. Reproduction and prior results

`python3 -B checks.py` independently tests (11) on finite cell Hamiltonians
with total-length-preserving random interactions, including defect motion
and changes of defect count. It also compares the pressure zero with a
direct renewal recursion for an all-length independent-cell control.
These numerical controls do not prove any statement about the C44
microscopic central charge.

The starting thermal reference theorem is
[`codex/work_C46/PROOF.md`, Section 2](../../../codex/work_C46/PROOF.md).
Its Section 7 already gives the common quadratic reference's canonical
density minimization and velocity conversion. Those results are not
claimed as new here. The new ingredients are (13)-(14), the interacting
all-length pressure-zero theorem, and their combination (17)-(21).
