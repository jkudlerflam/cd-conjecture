# Mobile representation and a charged-rail pin estimate

Date: 2026-11-16. Author: OpenAI Codex, GPT-6 family. **[Proved below,
not independently reviewed.]** This is a partial mobile result. The
vacuum-clock threshold, wrong-length matching, and CFT remain open.

## 1. Exact Hamiltonian and normalization

Use the literal `mobile_bond` of
[the 20/22-state family](../smaller_d_unbounded_20261116/model.py).
Fix N>=3 occupied material sites and V vacancies, L=N+V. In ordered
material coordinates let g_j be the vacancies after particle j,
A_j=1_(g_j=0), and S_j the unilateral shift raising g_j by one. The
legal coloring is unique. Apart from the physical cyclic quotient,
the exact Hamiltonian is

    H = H_ch + J sum_j A_j h_j - mu_occ N,
    H_ch = -t sum_j (S_(j-1) S_j^* + h.c.),
    sum_j g_j=V.                                         (1)

Here h_j is the entire material bond, restricted to legal material
words. No contact has been replaced by its expectation. Write that bond
as k_j+g c_j+lambda p_j+(o_j+o_(j+1))/2. Then (1)'s onsite part is

    (J/2) sum_j (A_(j-1)+A_j) o_j.                       (2)

In particular **all** of the following are contact-weighted: Delta on
charged letters; the extra D on express and terminal charged letters;
q^2 on M_v; mu_sel on markers; and f on L_1. Vacuum rank-one edge terms,
rail hopping, conversion, junctions and ports each have their own A_j.
Only -mu_occ N is an unweighted physical onsite scalar at fixed N.
Replacing (2) by J sum o_j is a different model.

Every forbidden-arrow projector is exactly conserved in this particular
22-state model. Thus (1) is an exact invariant legal-space restriction
at every finite K, with no t^2/(JK) leakage. This fact does not exclude
other, illegal sectors. All theorems below are within the legal space.

At lambda=0 each cell's vacuum/charged-rail flag is conserved. A charged
cell contains exactly one D, X, or excited terminal M, while a vacuum
cell contains V or M_v. These flags are different from vacuum-clock
energy eigenstates; the latter are not conserved by the contacts.

## 2. First-order reference, at fixed density

For the zero-source free gap sea at fixed N,V, put rho=N/L. Its energy
and tagged contact are

    E_ch(N,V) = -2t sin(pi N/L)/sin(pi/L),
    alpha_(N,V) = rho - sin^2(pi rho)/(N L sin^2(pi/L)).   (3)

The energy is the hard-core hopping sea. The contact formula follows
from its two-point determinant for the total adjacent-occupied count,
divided by N, and cyclic symmetry of the gap ground state. Thus

    alpha(rho)=rho-sin^2(pi rho)/(pi^2 rho)>0,
    P_ch H P_ch = E_ch + J alpha_(N,V) H_mat - mu_occ N.  (4)

This is an exact compression, not a bound on the full ground energy.
The energy per **physical** site of this reference is
e_ch(rho)+J rho alpha(rho) e_mat, where
e_ch=-2t sin(pi rho)/pi and e_mat is per material site.

At zero port and f>=0 the all-zero-bit cell ground energy is mu_sel+E_k.
A contiguous balanced pair (n-1,n+1) in place of (n,n) preserves both
N=2(n+1) locally and the marker count. Its first-order excess is

    J alpha_(N,V) [E_(n-1)+E_(n+1)-2E_n]
      = J alpha_(N,V) [delta_(n-1)+delta_(n+1)] > 0.      (5)

The common free charge energy, chemical potential, and bulk material
energy cancel at the same N,V. There is no density-change term hidden
in (5). In general a legal pattern has reference excess
J alpha sum_wrong delta_k. A scalar bulk subtraction only cancels
extensive errors if its interacting remainder is compared relatively.
We have not proved that remainder for wrong lengths.

## 3. A canonical scalar impurity lemma

This addresses the suggested smaller test in a fixed-vacancy ensemble.
It is not the arbitrary-quantum-source, positive-tilt theorem of C92.
Let P=A_0, Q=1-P and add uP, u>=0, to the free gap Hamiltonian at fixed
N,V>=1. The isometry S_0 from V-1 to V vacancies has range Q and obeys

    S_0^* H_ch(N,V) S_0 = H_ch(N,V-1).                   (6)

Thus the *compressed no-contact* space has a finite energy cost

    Delta_p=E_ch(N,V-1)-E_ch(N,V)>0,                     (7)

although the full charge excitation gap vanishes with volume. Also
||P H_ch Q||<=2t; only two gap hops cross P,Q.

For a ground state of H_ch+uP let a_u=<P>. Its variational energy is
at most E_ch+u. The diagonal blocks of H_ch are bounded below by
E_ch P and (E_ch+Delta_p)Q respectively. Cauchy--Schwarz therefore gives

    (Delta_p-u)(1-a_u)<=4t sqrt(a_u(1-a_u)),
    a_u >= (Delta_p-u)^2/[(Delta_p-u)^2+16t^2], u<Delta_p.

Integrating the ground-energy derivative (or its one-sided versions)
proves

    E(u)-E(0) >= u Delta_p^2/(Delta_p^2+64t^2),
    0<=u<=Delta_p/2.                                    (8)

There is no denominator involving the many-body charge gap.
At fixed rho, Delta_p tends to the strictly positive pressure

    p(rho)=(2t/pi)[sin(pi rho)-pi rho cos(pi rho)].        (9)

A conservative explicit uniform version is: N>=128,
1/3<=rho<=2/3, 0<=u<=t/16 imply E(u)-E(0)>=u/4097.
Indeed p(rho)>=t(sqrt(3)/pi-1/3)>0.217t. Writing the continuum
energy N epsilon(V/N), the exact sine formula differs from it by at
most 5t/N for N>=3: use |epsilon|<=2t and
x/sin x-1<=x^2/(6-x^2), x=pi rho/N<=pi/3. Convexity of epsilon
then gives Delta_p>=0.217t-10t/N>=t/8 when N>=128.

This scalar result is new relative to the stated C91–C96 results, but
does **not** transfer by substituting an operator-valued clock source
for u. Other material terms are extensive and correlated with this
contact. Their cancellation would require a further argument.

## 4. A contact conditioned on the moving charged head

Work first in independently labelled periodic gaps, with all V and
positive length tilt s. This is an auxiliary ensemble; sV is a scalar
in each fixed-N,V physical problem, not a new microscopic interaction.
Allow a uniform phase phi on the gap hopping. Fix a material cell of
b=n+1 sites, and compress it to its charged sector P_ex. All other
cells and their quantum states are retained. Ports incident to this
cell vanish under P_ex, since a port changes its charged flag.

Let Pi_i project onto the head occupying material site i of this cell,
including its terminal marker. On P_ex, sum_i Pi_i=1. Define

    S=sum_(i in cell) Pi_i S_i,
    P=sum_(i in cell) Pi_i A_i,
    S^*S=1,  SS^*=1-P,
    Phi(rho)=P rho P+S^* rho S.                          (10)

This is a two-Kraus channel. It removes one vacancy immediately after
the *actual quantum head* when that gap is nonempty. It does not replace
<Pi_i A_i> by <Pi_i><A_i>.

The charge Hamiltonian is diagonal in the Pi_i decomposition, so
S^* H_ch S=H_ch and ||P H_ch(1-P)||<=2t, even for phi!=0.
The tilt decreases by s(1-a), where a=Tr rho P. All source terms
disjoint from this cell and its two boundary bonds commute with the
channel. At most R=b+1 source bonds remain. If ||h_j||<=M, their
energy change is bounded above by 2JRM. For a tilted pinned ground
state the variational principle gives

    s(1-a) <= 4t sqrt(a(1-a)) + 2JRM.                   (11)

This argument permits arbitrary entanglement between the cell, the
gap fluid and other cells. It uses the ground state of the interacting
pinned model, not the free sea.

Consequently, on the sufficient window

    t/2<=s<=3t/4,       JRM<=t/32,

one has

    a>=a_*=1/1024.                                      (12)

For if a<=1/2, (11), together with 2JRM<=s/4, gives
4t sqrt(a)>=s/4 and a>=s^2/(256t^2). The case a>1/2 is immediate.
The conservative schedule has spare room. A fixed-gap-number ground
state also satisfies this argument if its pinned canonical chemical
potential E_pin(N,V-1)-E_pin(N,V) is at least s. That condition is not
assumed at an arbitrary prescribed density.

## 5. Replacement and the charged-rail pin cost

Assume the parent's frozen parameters, in particular
Delta=4+E_n, |g|<=1/2, q<=1/4, D>=q^2, and f>=0. The extra D condition
is harmless for its large-D band hierarchy. Put

    kappa=1-|g|-q^2>=7/16,     d_ex=kappa/1024.           (13)

For each charged cell basis word, measure that word and replace it by
the all-zero-bit vacuum word with the head at the same material site.
Call the resulting cell-local channel Gamma. It leaves the complete
gap marginal, the other-cell marginals, and every marker position
unchanged. It can lower total payload charge; the comparison sums over
payload charge, as does the original ground-sector problem.

Here is an operator/expectation proof of its energy improvement.
Bound each excited hopping edge by minus its absolute amplitude on
each endpoint. At a nonterminal head i the outgoing incidences sum
to at most A_(i-1)+(1+|g|)A_i. The Delta onsite halves supply
Delta(A_(i-1)+A_i)/2. The replacement's vacuum diagonal is at most
A_(i-1)+q^2 A_i; it has zero hopping expectation after the measurement.
The difference is therefore at least

    (E_n/2) A_(i-1) + (1+E_n/2-|g|-q^2) A_i
      >= kappa A_i.                                    (14)

Extra X energy only helps. At the excited terminal marker the
difference has left coefficient (D+E_n-q^2)/2>=0 and right coefficient
(Delta+D-q^2)/2>=2, again giving (14)'s lower bound. Terminal vacuum
energy q^2 is split between both contacts here, as required by (2).
Marker chemical terms are unchanged, and the nonnegative readout
decreases. Other-cell terms are unchanged. Incident ports have zero
expectation in both the pinned input and the vacuum-cell output, even
when lambda!=0 elsewhere. Thus exactly

    Tr H rho - Tr H Gamma(rho) >= J kappa Tr rho P.      (15)

Let F_pin(s) and F(s) be the *exact ground energies* with and without
this charged-cell pin, minimizing over V. They are not thermal free
energies in this notation. Equations (12),(15) prove

    F_pin(s)-F(s) >= J d_ex >= 7J/16384.                 (16)

The same bath and all untouched interacting material terms appear on
both sides of (15). Their extensive bulk cancels exactly. No hopping
link is cut and no un-subtracted O(J^2 N/t) error is used.

At lambda=0, if all other cells are required to be vacuum cells, (16)
compares a single charged-rail cell with the all-vacuum material sector.
At nonzero lambda it is a **compressed pin bound**, not a many-body
spectral gap: the full ground state can contain virtual charged pairs.

For a full tilted ground state psi, let B=sqrt(2+2 eta)|lambda| be the
material port norm. Only the two incident ports connect P_ex and 1-P_ex,
so ||P_ex H(1-P_ex)||<=2JB. The pinned block lies at least J d_ex above
the full ground energy by (16). Projecting the eigenvalue equation gives

    <psi,P_ex psi> <= (2B/d_ex)^2.                       (17)

This is a volume-independent bound on **charged-rail admixture** in the
interacting ground state, with the gapless charge dynamics retained.
It does not bound admixture of excited eigenstates within the vacuum
path, and does not supply a complete dressed auxiliary representation.

## 6. Which canonical densities this establishes

Let V_* minimize the full labelled E(N,V)+sV. For this specific vacancy
number, F(s)=E(N,V_*)+sV_* and
E_pin(N,V_*)+sV_*>=F_pin(s). Hence

    E_pin(N,V_*)-E(N,V_*)>=J d_ex.                       (18)

No change of density occurs in this comparison. The same argument
applies to the vacuum/single-charged comparison at lambda=0. It does
not prove (18) for every preassigned V/N. Equivalence at an arbitrary
specified density is a remaining assumption, not hidden in notation.

The minimizing densities stay away from zero and one. In every sector
the source norm is at most JMN, so energies differ from the free charge
energy by at most JMN. Write u=V/N and
epsilon(u)=-(2t/pi)(1+u)sin(pi/(1+u)). Against the trial u=1, the free
tilted energy per material site outside 1/2<u<2 exceeds the trial by
more than t/16, uniformly for s/t in [1/2,3/4]. At the endpoints the
differences are

    t[4-3sqrt(3)/2]/pi-s/2 > t/16,
    t[4-3sqrt(3)]/pi+s > t/16.

Convexity controls the rest of the exterior intervals. The finite-N
error is at most 10t/N^2 in a comparison, using §3's sine estimate.
Thus N>=32 and JRM<=t/32, with R>=3, force

    1/3 < N/(N+V_*) < 2/3.                              (19)

Every subsequential limiting density is interior. This establishes
canonical sequences selected by a fixed positive tilt, not a uniqueness
or ensemble-equivalence theorem for all densities.

## 7. Physical cyclic gluing: a vanishing global energy allowance

This section concerns ground energies only, not thermal multiplicities.
Include all cyclic root frames of the chosen ordered material sector.
For a distinguished pin, include its cyclic images as well. Let R
rotate gaps and material letters to the left, and put W=sum_j j g_j.
A physical momentum P=2pi m/L state can be written on rooted data as

    Psi(x,g,a)=exp(iP x) exp(-iP W/N) phi(g,a).

Changing the root to the next particle shifts x by 1+g_0 and W by
-V+Ng_0. Thus the precise gluing condition is

    phi(Rg,Ra)=exp(-iP L/N) phi(g,a).                    (20)

The gap hopping in this frame has uniform phase P/N (with the opposite
phase for its Hermitian conjugate); all A_j h_j are unchanged. This
derives the twist from the literal mobile bond. A cyclic projection at
zero twist alone would omit physical momentum sectors.

Let E_lab(phi) denote the unrestricted rooted ground energy with uniform
phase phi, at fixed N,V. The large gauge transformation exp(2pi i W/N)
shifts phi by 2pi/N and commutes with every material source and pin.
Write T_phi for the hopping and I_phi=partial_phi T_phi. Applying both
large gauge transformations to a ground state gives

    |<I_phi>| <= 2tN tan(pi/N),   N>=3.                  (21)

Indeed T_(phi+z)=cos(z)T_phi+sin(z)I_phi, and ||T_phi||<=2tN;
the two variational inequalities with z=+/-2pi/N give (21).
Integrating ground-energy one-sided derivatives controls the unrestricted
E_lab(phi). Conversely choose a simultaneous ground eigenvector of the
zero-flux Hamiltonian and cyclic rotation. There is a physical momentum
P, |P|<=pi, satisfying (20) for its rotation eigenvalue. Its trial energy
changes by at most

    2tN(1-cos(P/N)) + 2tN tan(pi/N)|sin(P/N)|
      <= 48t/N.                                        (22)

The same constant bounds the unrestricted energy change below. Therefore
the physical cyclic sector and the unrestricted rooted sector obey

    |E_phys(N,V)-E_lab(0;N,V)|<=48t/N.                   (23)

The estimate is uniform in V, J, material parameters and source pattern;
the sources need only commute with the gap-number gauge transformation.
It also holds for the cyclic orbit of a single charged cell at lambda=0.
This is a **single global O(t/N) allowance**, not O(t) per defect.

For example, the physical one-charged-cell versus vacuum comparison at
lambda=0, in the tilted ensemble and in canonical minimizing blocks,
is at least J d_ex-96t/N. It is positive for all sufficiently large N
at fixed J/t>0 and has thermodynamic lower limit J d_ex.

For nonzero lambda, a physical tilted ground state is an eigenstate of
one rooted uniform-flux block. Its excess above that block's unrestricted
tilted ground energy is at most 96t/N, by (21)–(23). Apply the pinned
block proof there. For any specified cell in the root frame,

    <P_ex> <= [2JB/(J d_ex-96t/N)]^2,
    provided J d_ex>96t/N.                              (24)

The cell-average is root independent and gives a physical charged-cell
density bound. A root-frame pin on a named cell should not itself be
mistaken for an extra onsite physical degree of freedom.

## 8. Compatible hierarchy and exact scope

A sufficient legal material bond norm bound is

    M=Delta+D+|mu_sel|+f+4+B.

The constant 4 bounds the remaining diagonal vacuum and kinetic rows;
the independent finite rebuild checks it as well. Fix the parent's
parameters and positive readout first, with D>=q^2. One may choose

    B<=min((1-q^2)^2 q^(2n+4)/10, d_ex/16),
    0<J/t<=1/[32(n+2)M],
    K arbitrary finite positive in the legal restriction. (25)

The first allowance is compatible with the completed frozen comparison.
For N>=192t/(J d_ex), (24) is at most (4B/d_ex)^2<=1/16.
At fixed nonzero couplings its limit is bounded by (17). These are
sufficient pin/admixture budgets, not a complete matching schedule.
Finite K introduces no legal leakage but still needs its own illegal
sector energy comparison. No such comparison is claimed here.

The new result controls a charged-rail auxiliary pin, including a moving
head correlated with vacancies and nonzero ports in the unpinned model.
It does **not** control vacuum-clock excited eigenstates. Replacing a
charged head removes the large positive Delta,D contact energy; that
mechanism is absent for orthogonal vacuum-clock states. Nor does the
replacement alter a material marker pattern, so it cannot compare wrong
cell lengths or establish an o(J gamma) balanced-defect remainder.

For arbitrarily long wrong cells the local channel's support grows with
their length; (25) is not a uniform theorem for pinning all such cells.
Other cells may be arbitrarily long in the stated pin theorem, but the
cell being pinned has the fixed selected length b. This distinction
must be retained in any arbitrary-n construction.

C91–C96 and the
[route audit](../route_audit_20261109/README.md) explain the unresolved
relative response and clock problems. C94's finite-temperature crossing
counterexample is avoided here by zero-temperature variational arguments
and explicit ensemble qualifications; it is not repaired by this proof.
C96's counterexample to a temperature-uniform Taylor remainder also
prevents inferring an all-order thermal theorem from our finite checks.

The total nonchiral Virasoro charges c=n and c=n+1 remain conditional.
No multichannel CFT, vacuum-clock dressed gap, all-defect selection, or
unbounded-c counterexample is concluded.
