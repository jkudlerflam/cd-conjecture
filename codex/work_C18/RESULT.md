# C18. A hopping selector with idle ports, and the remaining energy compiler

This is a partial construction. The history clock can select the cell
length with no separate walker, including all malformed clock sectors.
It does not yet give a fixed-dimensional implementation of the unbounded
route-A target. No value of d0 for unbounded c is established here.

[Proved] Adapting C8's full-sweep clock to two distinct, permanently idle
qubit ports gives a repeated nearest-neighbor bond of dimension

    d_clock+work = 15w + 4,
    b = n + 2,
    K(b) = 2b 2^(b-2) + b - 7.

Here w is the complete bulk work dimension, n is the number of bulk
sites, b is the physical period, and K is the number of transition
edges. The computation graph has C_n=K+1 vertices. The work dimension
must include every data, tape, lane, instruction and event register.
It is not an unspecified constant that can be omitted from d0.

[Proved] A bare positive Bell-check implementation has an additional
obstruction. If its energy terms annihilate the intended code at every
time, the complete port-source resolvent is scalar on the logical
space. Its second-order interaction is therefore scalar. Positive
span identities alone do not provide an exact route-A code.

[Proved construction, numerical spectrum] Giving a bulk site a partner
qubit and a second qubit avoids this restricted obstruction. Section 8
specifies a literal d0=64, b=4 cell whose clock both runs a SWAP program
and selects its length. It has an exact two-dimensional code, zero
first order, and a nonzero second-order ZZ coefficient. It is a small
diagonal gadget, not a large-central-charge model.

[Proved specification, missing local implementation] The detuned
r-flavor free-fermion target requires only 3r+1 operator slots in C9,
with whole-cell dimension at most

    D_r = 4 ceil((18r+19) 2^r / 4).

A uniform algorithm generating this target from b is supplied. The
finite-state energy compiler, its initialization of the logical region,
and its complete alphabet are still missing. Consequently neither the
estimate d0 approximately 45 nor an alternative fixed d0 for unbounded
c follows.

## 1. Attribution, scope, and conventions

[Assumed inputs, with proofs in the supplied repository] C7 proves the
cell-density comparison for its walker. C8's current
`initialized_binary_clock.md` and `sweeping_clock.md` prove the nine- and
15-letter counter clocks. C8 was already developing their work-machine
interface when this run began. Those clocks are used here, not claimed
as new C18 results. The rule list is repeated so the two-port extension
and its checks are self-contained. C9 supplies the spectral-weight
dilation for any Hermitian target bond; G1 supplies the earlier swap
construction. Neither is re-proved here.

C18's additions are the two-port direct-sum wrapper and count, exact
midpoint margins for its period convention, the hopping initialization
scale, the weighted port-interface analysis, the Bell-check obstruction,
the explicit d=64 small cell, and a restricted uniform target specification.
C19/PROGRESS.md was absent initially. It appeared during this run and
was then read in full. Its uniform Jacobi-filter construction supersedes
the need to develop the target matrix recipe here. No network,
subagents, cluster computations, commits or modifications outside this
folder were used. The parent project STATE.md was read but not modified,
in accordance with the explicit scope restriction. PROGRESS.md carries
this run's memory.

Every mathematical assertion below is marked proved, numerical,
heuristic, or assumed, including inherited results. Numerical checks
refer to finite graphs or cells, never to a thermodynamic CFT. Units
are set by the hopping coefficient J>0. The legality coefficient is
kappa, distinct from the step count K. Ports are P_L and P_R; R0,R1
below are return-head letters and do not denote ports.

## 2. An explicit two-port transition table

[Proved construction] The on-site Hilbert space is the direct sum

    H_site = (C^15 tensor C^w) direct_sum C^2_PL direct_sum C^2_PR.

The 15 bulk letters are

    L, 0, 1, C0, C1, R0, R1, I, Z,
    A0, A1, D0, D1, E0, E1.

Each carries the same work space W=C^w. Each port letter carries only
its qubit. The ports have no clock or work tensor factor. There are
nine head letters, C0,C1,R0,R1,I,D0,D1,E0,E1.

The allowed bulk pairs are exactly

    LL,
    Lh                       h=C0,C1,R0,R1,I,
    hx                       h=C0,C1,R0,R1, x=0,1,
    xy                       x,y=0,1,
    IZ, ZZ,
    Aa Ab, Aa Db, Aa Eb       a,b=0,1,
    Da x, Ea x               a,x=0,1.

The allowed pairs involving ports are exactly

    P_R P_L,
    P_L x     x=L,C0,C1,R0,R1,I,A0,A1,D0,D1,E0,E1,
    x P_R     x=0,1,Z,C0,C1,R0,R1,I,D0,D1,E0,E1.

All unlisted pairs have penalty kappa. In particular P_L P_R,
P_L P_L and P_R P_R are forbidden. A legal cell is
P_L [n bulk sites] P_R, and adjacent cells meet on P_R P_L.
Thus its physical length is n+2. This differs from C8's scalar-marker
period n+1.

The forward rules, always accompanied by their adjoints, are

    C1 x       -> L Cx                  x=0,1,
    C0         -> R1,
    L Rx       -> R0 x                  x=0,1,
    I Z        -> L I,
    I P_R      -> R0 P_R,
    P_L Rx     -> P_L Dx                x=0,1,
    Da x       -> Aa Dx                 a,x=0,1,
    Dx P_R     -> Ex P_R                x=0,1,
    Aa Ex      -> Ea x                  a,x=0,1,
    P_L Ex     -> P_L Cx                x=0,1.

There are 21 rule instances, including the one-site rule. This count
is separate from the 15 bulk letters.

For a rule a->a' put

    h_rule = -J ( |a'><a| tensor U_a + |a><a'| tensor U_a^dagger ).

On a two-bulk-site rule U_a acts on W tensor W. On a boundary rule it
acts only on the adjacent bulk work space and is identity on the port
qubit. On the one-site turn it acts on W. The I Z rule is identity on
work, so the entire initialization sweep precedes nontrivial gates.
The activation I P_R -> R0 P_R may seed a work state after the last
initialization check. All gates are fixed by the rule, not by the
site, cell length or ring length.

The repeated bond is the sum of all two-site rules, kappa times the
forbidden-pair projector, optional clock-diagonal initialization and
energy terms, and half of each one-site term on each endpoint. Assign
chemical energy mu/2 to each port letter. Its half-onsite bond form
has mu/4 times the sum of the two endpoint port projectors. Thus each
complete cell pays mu once. Include the one-site C0/R1 hopping once
per site by the same half-onsite convention.

The intercell interaction, when used, is supported only on P_R P_L.
For route A it can be lambda sum_(a=1)^3 sigma_a tensor sigma_a.
It is absent in the selector and isolated-cell theorems below.
Every propagation term commutes with every physical port operator,
including in malformed sectors. This is an operator identity on the
full direct sum, not a projection onto a presumed legal sector.

## 3. Complete path and malformed-sector proof

[Proved, clock combinatorics inherited from C8] A legal segment has
one of the following forms, with head position i=1,...,n:

    L^(i-1) h x_(i+1)...x_n,    h=C0,C1,R0,R1, x_j in {0,1},
    L^(i-1) I Z^(n-i),
    A_(x1)...A_(x_(i-1)) h_(xi) x_(i+1)...x_n, h=D or E.

The left port allows a prefix or a head, but no suffix. Prefixes can
end only in their specified head. A prefix-only word cannot end at
the right port. A head permits only its specified suffix, and no
suffix can contain another head. This proves necessity of the list.
Every listed word obeys all pair constraints, proving sufficiency.
There is no extra input length, tape split, mode count or second
counter in this classical graph.

The first word is I Z^(n-1). Its head sweeps right, then activates R0
at the last site. The return reaches binary zero at the left. C1
carries through trailing ones, writing L, while C0 -> R1 turns the
first zero into a one. Return through the zero prefix restores its
ordinary binary digits. The all-ones carry ends at L^(n-1) C1.
Before each carry starts, the D and E rules insert a full outward
and return sweep preserving all binary digits. There are 2^n such
sweeps.

Every counter configuration lies on this orbit. For a C head, its
erased prefix specifies trailing ones and its suffix specifies the
remaining integer. R1 is the unique turnaround at its position.
For R0, the first one to its right determines the turnaround from
which it returns; no such one means the initial zero return. Each
D/E configuration is fixed by its binary word and sweep position.
The I configurations occur in their single initial sweep. Therefore
all listed words occur once and there is no other legal component.
Counting them gives

    C_n = n + 4(2^n-1) + 2n 2^n
        = (2n+4)2^n + n - 4,
    K = C_n - 1.

For permutation work gates, every work-basis input also follows a
literal path of this length. General quantum gates can branch in a
chosen work basis, but unitary gauging gives a path operator tensored
with the input identity. The spectral statement needed for selection
therefore holds without calling those basis graphs simple paths.

A shorter computation on the work register does not shorten this
clock. Work-machine cycles also do not create clock cycles. The work
unitaries decorate a fixed open path and can be gauged away on that
path. If several logical input vectors are retained, each has a copy
of the same clock spectrum. "Unique longest path" means the unique
classical clock skeleton, with the intended quantum input multiplicity.
It cannot mean a one-dimensional entire code for route A.

[Proved, all illegal clock sectors] Each rule has identical allowed
predecessor sets on its two sides and identical allowed successor
sets. Its internal bond is legal on both sides. C0 and R1 have the
same neighbors. Boundary substitutions retain the same port and
preserve the adjoining bulk constraint. Hence every individual
forbidden-bond projector commutes with every rule, on the full
Hilbert space. Port positions are also conserved. Work decorations
do not alter these statements.

Cut an open segment at each of its v violated constraints, including
its boundary violations. Each piece has at most one head, so there
are at most v+1 heads. Each head participates in at most two rule
instances. The inequality

    2 Re <u,Uv> <= ||u||^2 + ||v||^2

applied to each edge gives H_hop >= -2J N_head, also for arbitrary
unitary work decorations. In a v-defect sector,

    H_hop + kappa v >= -2J + (kappa-2J)v.

Thus kappa>=4J puts every illegal open-segment sector at nonnegative
energy, while the legal hopping ground energy is negative. Cycles,
multiple heads, missing heads, counter fragments and incorrect
boundary letters are included. The penalty is present at every
configuration of that component. It is not a single visited check
whose effect vanishes as the component gets longer.

For a ring, cut at its v forbidden bonds. Complete P_L...P_R cells
have their full cell energy. There are at most two incomplete end
fragments per cut string. Each has at most one head and costs at
least -2J; unassigned port energies are nonnegative when mu>0.
Propagation cannot cross a forbidden cut, and each rule involving a
port acts only on the bulk on its own side. Positive checks can only
raise the bound. If e_*<0 is the least complete-cell energy density,

    H_ring >= e_* L_ring + (kappa-4J)v.

A legal marker-free ring has no head and zero hopping energy. It can
consist of L's, Z's, binary suffix letters, or A-prefix letters.
Such sectors lose to a tiling with e_*<0. A legal marked ring
necessarily decomposes into complete cells. This exhausts the ring
sectors for the unperturbed selector. For incompatible ring lengths,
the minimizing mixture of cell lengths still needs consideration;
no all-ring conformal conclusion is made.

## 4. Hopping histories and initialization

[Proved] Let U_t be the cumulative work circuit, t=0,...,C-1. The
legal hopping operator is unitarily equivalent to a path adjacency
matrix tensored with its input space. Its energies and ground weights
are

    E_k = -2J cos(k pi/(C+1)),                 k=1,...,C,
    phi_t = sqrt(2/(C+1)) sin((t+1)pi/(C+1)),
    f_t = phi_t^2,
    omega = 4J sin(3pi/(2(C+1))) sin(pi/(2(C+1))).

The sine-product expression avoids loss of precision in omega.
The ideal history is V_phi psi=sum_t phi_t |t> U_t psi.
All f_t are strictly positive, but they are not uniform.

Product workspace pins can be imposed while I passes the appropriate
site. Uniform pins use |I><I| tensor p; additional left/right checks
can use P_L I or I P_R. Assume the actual checks select commuting
input projectors. Because the initialization sweep is idle, gauging
away U_t turns a failed check into eta |s><s| at its check time.
Accepted inputs keep the unshifted E_1 exactly. Every rejected input
has strictly larger energy: equality would require the strictly
positive path ground vector to vanish at a check.

For a quantitative bound let

    f_min = 2 sin^2(pi/(C+1))/(C+1).

Since H_path-E_1 >= omega(1-|phi><phi|), adding a single eta pin and
diagonalizing the span of phi and |s> gives

    g_in >= omega eta f_min / (omega+eta).

Indeed the smaller eigenvalue of the comparison operator is

    (omega+eta - sqrt((omega+eta)^2-4 omega eta f_s))/2,

which is at least omega eta f_s/(omega+eta). Other checks only raise
it. The resulting finite-cell gap above the accepted history space
is at least min(omega,g_in), before adding energetic code terms.
If there are no rejected inputs, omit g_in. For eta proportional to
J this conservative lower bound is order J C^-5.

[Proved sharper endpoint example] With a pin eta=J at t=0, the
boundary equation is represented by a ghost value z_-1=-z_0. The
right endpoint still has z_C=0. The positive ground mode is

    z_t proportional to sin((t+1/2) theta),
    theta=pi/(C+1/2).

It solves every interior and endpoint equation, so Perron positivity
identifies it as the ground mode. Its energy is

    E_pin=-2J cos(pi/(C+1/2)),
    E_pin-E_1 = J pi^2 C^-3 + O(J C^-4).

The variational upper bound for any endpoint pin is eta f_0, of
order eta C^-3. Thus the initialization gap cannot in general be
identified with the C^-2 clock gap. This also answers why putting
one check on an arbitrarily longer illegal path would not suffice.
The clock proof above excludes those paths instead.

[Assumed until implemented] Pins for a work interpreter must select
its actual input tape and the correct logical region. Uniform product
pins alone do not select r(b) free qubits while clearing every other
qubit. C8 describes later computed checks as an interface, but this
run does not supply their route-A transition table. Their gap must
be proved when they are added.

## 5. Length selection, margins and gaps

[Proved] Set b=n+2 and

    a_n=2J[1-cos(pi/(C_n+1))],
    delta=2J-mu,
    e_n=(a_n-delta)/(n+2),
    delta_n=(n+3)a_n-(n+2)a_(n+1).

The energy a_n is computed stably as
4J sin^2(pi/[2(C_n+1)]). It is decreasing and strictly convex.
For completeness, put x>=1,

    u=(2x+4)2^x, y=u+x-3,
    k=log(2)+1/(x+2), beta=log(2)^2+2log(2)/(x+2).

Then y'=ku+1, y''=beta u, and direct differentiation gives

    a''(x) = 2J pi/y^4 *
       [(2y(y')^2-y^2 y'') sin(pi/y) + pi(y')^2 cos(pi/y)].

Both trigonometric factors are positive, since y>=10. Moreover

    2(y')^2-y y''
      = beta u(u-x+3) + 2u^2/(x+2)^2 + 4ku + 2 > 0.

This proves the asserted convexity. Therefore

    delta_(n-1)-delta_n = (n+2)(a_(n-1)-2a_n+a_(n+1)) > 0,
    e_j-e_(j+1) = (delta_j-delta)/[(j+2)(j+3)].

Choosing delta_n<delta<delta_(n-1) selects n uniquely among every
positive integer bulk length. For n>=2, delta<delta_1<4a_1<2J,
so mu>0. Also delta>delta_n>a_n, so e_n<0, as required for the
illegal-sector comparison.

At the midpoint delta=(delta_(n-1)+delta_n)/2, define

    gamma_len = (a_(n-1)-2a_n+a_(n+1))/2.

Then the two adjacent cells have identical excess above the selected
bulk density:

    (n+1)(e_(n-1)-e_n) = (n+3)(e_(n+1)-e_n) = gamma_len.

The excess function g(j)=a_j-delta-(j+2)e_n is strictly discretely
convex, g(n)=0, and its adjacent values are gamma_len>0. Its discrete
slopes are increasing. Moving farther on either side therefore only
increases g. This proves, without a finite search cutoff,

    min_(j!=n) (j+2)(e_j-e_n) = gamma_len.

The density sequence decreases up to n and increases after n. Thus
its least competing density is at an adjacent length. At the midpoint
it is the longer neighbor:

    Delta_density = min_(j!=n)(e_j-e_n) = gamma_len/(n+3).

For large n, C_n approximately (2n+4)2^n and a_(n+1)/a_n tends to 1/4.
Consequently

    a_n ~ J pi^2 / [(2n+4)^2 4^n],
    gamma_len ~ 9a_n/8,
    omega ~ 3a_n,
    delta ~ (15/8)(n+2)a_n,
    Delta_density ~ 9a_n/[8(n+3)].

These are asymptotic statements, not small-n equalities. On compatible
rings L_ring in (n+2)Z the gap above the entire retained tensor-product
history space is bounded below by

    min{gamma_len, omega, g_in, |e_n|(n+2), kappa-4J},

with absent initialization checks omitted. Subsequent energetic
selection of a smaller code adds its own gap requirement. The bound
is positive for each fixed n, not uniform in n.

[Proved limitation of the generic suggestion] Monotonicity of K(b)
alone is insufficient for C7's argument. For three consecutive periods
b=3,4,5, take C=(9,9,29). Then a_3=a_4>a_5 and

    delta_3=4a_3-3a_4=a_3,
    delta_4=5a_4-4a_5>a_3.

There is no delta with delta_4<delta<delta_3, so period 4 is never
uniquely selected. The actual full-sweep K satisfies convexity by
the preceding calculation. A program that inserts uncontrolled
length-dependent delays would require a new calculation.

[Proved positive-energy tolerance] Suppose a clock-diagonal,
cell-confined positive program D preserves the classical defects.
Its energy at every other length is bounded below by the bare
clock energy. Suppose at the selected length an accepted trial
history has additional energy at most v_n. Then selection survives if

    v_n < (n+2) gamma_len/(n+3).

Indeed the actual selected density is at most e_n+v_n/(n+2), below
every bare competing density. It is still negative. Re-run the
fragment proof with this actual minimum density; positive terms on
fragments or cut bonds can be discarded in a lower bound. Thus the
all-illegal-sector bound survives as well. If D annihilates an
accepted code at the selected length, v_n=0 and the original
selection is unchanged, regardless of its positive norm elsewhere.
For signed programs a bound |E'_j-E_j|<epsilon(j+2) preserves the
ordering of complete cells when 2epsilon<Delta_density. This last
statement is only a complete-cell comparison; signed cut-bond terms
need their own fragment estimate.

A small energy perturbation can preserve length selection while still
splitting the intended code or producing nonzero first order. These
are separate obligations.

## 6. The weighted route-A interface

[Proved] The port-idle identities do not require a uniform clock vector.
For any port operator F, every U_t commutes with F, so

    F V_phi = V_phi F,
    V_phi^dagger F V_phi = F.

If accepted input constraints preserve the ports, the same statement
holds on their accepted input space. Thus hopping weights do not dilute
the ports. They do change every energy-program average:

    D = sum_t |t><t| tensor h_t,
    A_bar = V_phi^dagger D V_phi
          = sum_t f_t U_t^dagger h_t U_t.

One cannot use the uniform-history factor C without correction. For
instance at C=5 the weights are (1/12,1/4,1/3,1/4,1/12).
At the abstract scheduled level, replacing a contribution A_t by
h_t=U_t A_t U_t^dagger/f_t reproduces sum A_t. This identity is not a
local rule compiler: computing t and C, representing 1/f_t, generating
the gates, and selecting the corresponding energy event all need
physical implementation. Endpoint f_t can be of order C^-3.

Suppose the shifted initialized hopping Hamiltonian H0 has kernel
P=V_phi V_phi^dagger and gap Gamma on Q=1-P. Let ||D||<=M, let A_bar
have exact code G and positive gap Delta, and require

    h_t U_t G = 0                       for every t.

Then D V_phi G=0. Conversely, since phi_t is nonzero, this timewise
condition is necessary for these particular histories to remain
zero eigenvectors of H0+D. Adding a scalar to the total energy does
not remove the corresponding componentwise eigenvector condition.

If Gamma-M>M^2/Delta, the code is the full ground space. To see this,
write a vector orthogonal to the code as x+y with x in P and y in Q.
Its quadratic form is bounded below by

    Delta ||x||^2 - 2M||x||||y|| + (Gamma-M)||y||^2.

The displayed condition makes this two-by-two form positive definite.
Its smaller eigenvalue is an explicit gap lower bound. For D>=0 and
timewise annihilation, uniqueness also follows directly from the
intersection of kernels, with a quantitative bound obtained similarly.

For z>=0 eliminate Q in the resolvent. The Schur correction has norm
at most

    epsilon(z)=M^2/(Gamma-M+z).

A port source F V_phi G lies entirely in P and, when G^dagger F G=0,
in the complement of the code within P. If epsilon(z)<Delta+z, the
resolvent identity gives

    ||R_sim(z)-R_Abar(z)||
      <= ||F||^2 epsilon(z) /
         [(Delta+z)(Delta+z-epsilon(z))].

For two different port insertions replace ||F||^2 by their norm
product. The denominator includes the initialization gap in Gamma,
not just the bare clock gap. This proves convergence for the weighted
hopping histories under the stated exact-code hypotheses. It is the
same Schur-complement mechanism used in programmable_cell, with the
weights and initialization scale made explicit.

Without timewise annihilation, only approximate code and port
statements follow. A leakage of order M/Gamma can generate first-order
terms of order lambda M/Gamma, while second order has scale
lambda^2 rho^2/Delta. The former must be made smaller than the latter;
the code splitting must also be below that exchange scale. These are
finite-cell conditions, not an infrared stability theorem.

## 7. Why bare positive Bell checks do not synthesize the exact code

[Proved restricted obstruction] Consider a port-idle history code
G_t=U_t G and a shifted hopping Hamiltonian H0-E_1. Add only

    D_t = sum_alpha v_(alpha,t) Q_(alpha,t) + B_t,
    v_(alpha,t)>=0, B_t>=0,

where each Q is the complement of a maximally entangled state of one
qubit port and a partner factor, possibly conjugated by a unitary
acting on the non-port registers. Each B_t acts on neither port.
Assume D_t G_t=0 at every time, and that each port has a Bell check
with positive weight at at least one time. Positivity implies that every term
with nonzero coefficient annihilates G_t separately. If these are
intended as the complete ground histories, the second-order logical
interaction through traceless port operators is scalar.

Proof. First use a single port p and one Bell check Q. A vector in
ker Q has its p,s factors in the Bell state. For any traceless Pauli
F_a on p,

    (1-Q) F_a G_t = 0,
    Q F_a G_t = F_a G_t.

For a conjugated Bell check this remains true because its conjugating
unitary commutes with F_a. A check involving the other port commutes
with F_a and therefore annihilates F_a G_t. Likewise

    B_t F_a G_t = F_a B_t G_t = 0.

Different same-port sources have inner products
G_t^dagger F_a F_b G_t=delta_ab I. A time with an active Bell check
establishes this identity; port-idle propagation carries it to every
other time. Opposite-port source spaces are
orthogonal: the relevant own-port Bell check has eigenvalue one on
one of them and zero on the other. These facts also follow directly
from the two Bell factors when they use distinct partners.

For each port p and Pauli label a, the space

    S_(p,a) = span{ |t> F_a^p G_t |i> : all t,i }

is invariant under the full Hamiltonian. Propagation maps a source
at t to the corresponding source at t+1 because its gate is port-idle.
The energy terms act on it as the scalar
v_p(t)=sum of the active weights of that port's Bell checks. Hence
on S_(p,a) the exact Hamiltonian is

    [H_path-E_1 + diag(v_p(t))] tensor I_logical.

It is independent of the logical label and of a. In particular,

    G_hist^dagger F_a^p (H-E_1+z)^(-1)_exc F_b^p G_hist
      = delta_ab r_p(z) I_logical,

and cross-port sources are orthogonal reducing sectors. For a bond
sum_a F_a^R tensor F_a^L, the full two-cell intermediate operator is
also a pair of clock operators tensored with logical identity.
Its second-order compression is scalar. This proof is at finite J;
it does not discard clock-excited states. QED.

The hypotheses are restrictive and explicit. The theorem does not
rule out a non-frustration-free energy program, signed cancellation,
a term that sees an additional workspace factor beside the Bell
partner, or a code with different timewise conditions. In particular,
it does not refute minimal_d0's positive-span identity. Realizing
H_target+alpha I as a positive sum generally makes the target code a
positive-energy eigenspace of the sum, not the common kernel of its
summands. Exact finite-J histories then require another argument.

[Numerical illustration] A three-time port-idle program with cumulative
gates (I,SWAP_(s,l),I) and one static check Q_(p,s) has the weighted
average

    A_bar = (Q_(p,s)+Q_(p,l))/2.

Its lowest energy is 1/4, its ground space is two dimensional, and
the norm of the compressed port X is 1/3. Finite-J diagonalization
converges to these values, not to zero first order. This example
shows a failure of this particular positive-frame program, not an
impossibility for all non-frustration-free programs.

A useful escape is Q_(p,s) tensor A_c. If c is a second factor on the
partner's physical site, this is still a two-site physical term. It
annihilates Bell code states, but acts as the nontrivial operator A_c
on their port syndromes. The next section implements this escape.

## 8. A complete small d0=64 route-A cell

[Proved construction] Set w=4, with W=C^2_s tensor C^2_c on every
bulk letter. The full alphabet is

    15 letters carrying (s,c), plus two qubit port letters,
    d0=15*4+2+2=64.

Use the transition table of Section 2. On every outward rule
Da x -> Aa Dx, apply SWAP to the two c factors and identity to both
s factors. Every other rule has identity work gate. The port factors
are idle on the entire Hilbert space.

On each P_L x bond, for any bulk clock letter x, add

    Q_(pL,s) tensor (I_c+Z_c/4).

On each x P_R bond add the reflected term

    (I_c+Z_c/4) tensor Q_(s,pR),

with tensor factors ordered according to their physical registers.
The notation means the same Q times A on that bulk site, not an extra
physical site. Both terms are positive since A=I+Z/4>=3I/4.
On I P_R add eta |1><1|_c, pinning the rightmost c input while the
I sweep is still idle. These are all repeated two-site terms. Choose
J=eta=1, kappa=8 and mu=1.7970485670015128, the n=2 midpoint.

At the selected length n=2, the complete clock path has C=30 vertices.
The initialization check is at t=1. The four outward SWAP edges are
t=4,10,18,24, where t labels the edge from t to t+1. The exact code is

    sum_t phi_t |t>
      |Omega>_(pL,s1) |Omega>_(s2,pR)
      U_t ( |psi>_(c1) |0>_(c2) ).

Here psi is an arbitrary qubit and U_t alternates between identity
and SWAP as those four edges are crossed. Every positive boundary
term annihilates the code at every time. Conversely a ground vector
must be a hopping ground history, must pass the right-c input check,
and must have both Bell factors, since A is strictly positive.
These conditions leave exactly two states. They also prove zero
compression for every traceless operator on either port.

The entire legal cell Hilbert space has dimension

    30 * 2^2 * 4^2 = 1920.

Bell syndrome labels at the left and right are conserved by the
isolated cell. Their zero/nonzero types give four 120-dimensional
blocks with multiplicities 1,3,3,9. This exhausts the legal Hilbert
space; it is not a low-energy truncation. The (0,0) block is the clock
with the right-c input pin. Every other block is at least 3/4 above
the bare hopping ground energy. Thus a lower bound on the complete
cell gap is min(omega,g_in,3/4). Section 3 separately controls all
illegal clock sectors.

The selected cell has no energy shift from its boundary or input
terms. Every other length has nonnegative additional energy. Thus
Section 5 selects b=4 on every compatible ring with unchanged mu.
For n=1 the two Bell checks share one partner and cannot both vanish;
this only raises that competing cell's energy. For n>2 extra unused
bulk variables remain. They do not compete with the selected length,
but they are one reason this small program is not an unbounded-target
solution.

[Proved nontrivial second order] Work in the gauged clock frame and
fix the logical input bit i=0,1. Let H_c=H_path-E_1. On a right-port
syndrome the clock operator is H_R,i=H_c+diag(v_R,i), and similarly
on the left. For i=0 both c bits are always zero, so
v_R,0=v_L,0=5/4. For i=1 the occupied c bit moves between the sites.
The corresponding potential is 3/4 while the bit is at that end,
and 5/4 otherwise. Each end is occupied for a nonempty time interval.

Define f_(p,i)(tau)=<phi|exp(-tau H_(p,i))|phi>. The diagonal logical
bond through sum_a sigma_a tensor sigma_a is

    K2_(ij,ij) = -3 integral_0^infinity
                          f_(R,i)(tau) f_(L,j)(tau) d tau.

There are no off-diagonal logical entries, because all gauged checks
and potentials are diagonal in the initial c bits. Its ZZ coefficient
is

    beta_ZZ = -(3/4) integral_0^infinity
             [f_(R,0)-f_(R,1)] [f_(L,0)-f_(L,1)] d tau < 0.

To establish the strict sign, H_(p,1) has a pointwise smaller diagonal
potential than H_(p,0), with a strict reduction somewhere. Their
negative off-diagonal entries are identical. Exponentials of their
negatives are entrywise positive for positive time. The Duhamel
identity expresses exp(-tau H_(p,1))-exp(-tau H_(p,0)) as an integral
of two such exponentials around the nonnegative, nonzero diagonal
potential difference. Its expectation in phi is strictly positive.
Both bracketed differences are therefore negative for tau>0. The
integral converges since every syndrome potential is at least 3/4.
This proves a nonzero interaction without relying on its measured
coefficient.

[Numerical] With the coefficients above, the full legal-cell gap is
0.000941163372332468 and the code residual is 4.81e-16. In units in
which the effective bond is lambda^2 K2,

    K2 = -1.358800177739892 I
         +0.0764274635157273 Z tensor I
         +0.0937497178206302 I tensor Z
         -0.0113770035964675 Z tensor Z.

The four diagonal entries are
(-1.2,-1.364745428448327,-1.330100919838521,-1.540354362672717).
As an analytic normalization check, input 00 leaves both source
potentials at 5/4, so its entry is -3/(5/4+5/4)=-1.2.
As a separate check, exact two-cell eigenvalues were computed in
the invariant equal-Bell-syndrome sector. For each pair of logical
bits this sector has dimension 4*30^2=3600, and includes the entire
virtual orbit reached from that code state. At lambda=.004,.002,.001,
the errors to lambda^2 K2 are respectively
1.03844e-7, 1.29672e-8, 1.62007e-9. Error/lambda^3 approaches 1.62.
Eigenpair residuals are below 6.1e-15. These are exact invariant-sector
checks, not a diagonalization of the entire two-cell Hilbert space.

The model provides a concrete proof that an own-clock selector, an
exact code and a nontrivial route-A interaction are compatible once
the boundary term can see the extra c factor. Its effective model
here is diagonal. No claim about c>0 or arbitrarily many modes is
made for this example.


## 9. A target determined by b, and the compiler that remains

[Proved auxiliary specification] The supplied script target_spec.py is
a classical target generator, not a physical history machine. It takes
only the integer b>=3 and defines

    n=b-2,
    r=floor(log2(n+1)/10),
    N=2^r.

The r=0 case is a trivial target. The slow growth is a conservative
choice, not a proved compiler resource bound. It demonstrates that the
target need not have a separate externally chosen internal length.
No runtime, tape-size or precision guarantee for a finite-state
implementation is inferred from this Python function.

For r>=1 let gamma_(2f-1),gamma_(2f) be the ordinary Jordan-Wigner
Majoranas on r logical qubits, let P=product Z_f and eta_a=i gamma_a P.
Use the operator lists

    A=(eta_1,...,eta_(2r),Z_1,...,Z_r,I),
    B=(gamma_1,...,gamma_(2r),Z_1,...,Z_r,I).

There are m=3r+1 slots. For hopping magnitude one, phase theta_f and
chemical potential mu_f, the only nonzero entries of the real m-by-m
coefficient matrix are, with o=2f-1 and e=2f,

    K_oo=K_ee=-sin(theta_f)/2,
    K_oe=-cos(theta_f)/2, K_eo=cos(theta_f)/2,
    K_(2r+f,m)=K_(m,2r+f)=mu_f/4.

Then sum K_ij A_i tensor B_j is

    -sum_f [ exp(i theta_f) c_(j,f)^dagger c_(j+1,f) + h.c. ]
    -sum_f (mu_f/2)[ n_(j,f)+n_(j+1,f)-1 ].

Here the physical ordering of fermions is by cell, then by flavor;
the cell parity string is included in eta. Direct expansion gives
the displayed coefficient signs. The script independently constructs
all Jordan-Wigner matrices on two logical cells to verify them.

One explicit computable detuning is obtained from the first 2r primes
p_1,...,p_(2r), with Q=10r p_(2r):

    theta_f=pi sqrt(p_(2f-1))/Q,
    kF_f=pi/3+pi sqrt(p_(2f))/Q,
    mu_f=-2 cos(kF_f).

The two Fermi points are -theta_f plus or minus kF_f. There is no
nonzero integer linear combination of these 2r numbers equal to an
integer multiple of 2pi. Indeed divide such a relation by pi and
work in the multiquadratic field generated by these square roots.
Changing the sign of sqrt(p_(2f-1)) forces the sum of the two branch
coefficients for flavor f to vanish. Changing the sign of
sqrt(p_(2f)) forces their difference to vanish. All coefficients
are then zero. This is a kinematic nonresonance statement only.
The supplied stability analysis remains an independent physical input.

In C9 choose k=m+2=3r+3 energies. Its general dilation then fits in

    D_r=4 ceil(max(16N,(6k+1)N)/4)
       =4 ceil((18r+19)2^r/4).

[Assumed input] C9's positive weights and dilation produce the stated
abstract second-order bond. This count concerns an entire cell,
not d0. Its arbitrary completion also does not by itself preserve a
specified microscopic grading. The script checks the operator-slot
count, dimension inequality and target signs, not the missing
spatial implementation.

[Assumed from the live C19 result, preferred target] C19 now supplies
a stronger recipe for the same task: n+1 consecutive spectral energies,
Jacobi-filter coefficients, an explicit reflection dilation, exact
parity, and fields generated by the exchange itself. Its reflection
completion has whole-cell dimension 16(r+2)2^r. It also describes an
additive two-end alternative of dimension 16(r+2)^2 2^r, so no target
term needs both distant ports. Its coefficients and square roots have
uniform formulas. These are better inputs to the compiler than the
dense C9 completion. They were not re-derived or numerically replayed
here. The small target_spec.py calculation is retained as an auxiliary
checked formula from before that handoff, not as a competing target
synthesis project.

[Located obstruction, not an impossibility theorem] A completed local
machine must still perform all of the following within the fixed work
alphabet of Section 2:

1. Count n and derive r and every register boundary. Select exactly
   the intended logical input space while initializing all remaining
   work factors. Any later checks need a bound on their common kernel
   and their hopping gap.
2. Run the finite arithmetic and reversible routing for C19's
   coefficients and flags, with explicit time and space bounds below
   the clock's 2^n complete sweeps. Extra physical tracks or controller
   lanes must be included in w.
3. Turn the formulas into a fixed local energy menu on bulk bonds and
   the port/partner bonds. A positive bare Bell check cannot do this
   while retaining the timewise exact code, by Section 7. The extra
   factor used in Section 8 is a possible interface, but it is not a
   compiler for C19's syndrome matrices and flags.
4. Account for the sine weights of every event, or prove a different
   synthesis that uses them directly. An energy program correct for
   a uniform Laplacian history need not be correct here.
5. Preserve microscopic parity and the physical source resolvents.
   If the code is approximate, control its splitting and first-order
   port compression below the intended second-order exchange.

A universal gate set is insufficient to meet these conditions. A
finite collision table for a reversible work machine, such as C8's
interface, determines a circuit, but does not by itself supply this
energy program. Conversely the unique classical clock path already
ensures that arbitrary work configurations cannot create longer paths.
Once the actual program and checks are supplied, there is no remaining
freedom to improve the kinetic energy by changing a tape split.

## 10. Alphabet accounting and numerical record

[Proved counts] The following counts include the entire local direct
sum for the stated object. They are not upper bounds on a completed
unbounded-c construction.

| Object | Site dimension | What has actually been supplied |
| --- | ---: | --- |
| C8 binary clock with two idle qubit ports | 9w+4 | The same wrapper with the inherited nine-letter rules |
| Full-sweep clock with two idle qubit ports | 15w+4 | All rules, defects, selector and gap estimates above |
| Bare full-sweep clock and ports | 19 | w=1, no quantum work payload |
| Full-sweep clock with one bulk qubit | 34 | w=2, no complete route-A program |
| Section 8 partner and auxiliary qubits | 64 | w=4, complete two-state diagonal route-A example |
| Unbounded route-A target | not established | No complete work alphabet or energy transition table |

The number 45 cannot be inferred by replacing w with an estimated
universal-machine size. The clock and work alphabet here have a
different architecture from BCO. This work neither proves that 45 is
possible nor supplies a lower bound excluding it. The d0=64 example
is also not a bound for unbounded c.

[Numerical] The four scripts read no external data. They use Python,
NumPy and SciPy and write only their JSON outputs in this folder.
The final commands are

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clock_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 interface_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 routeA_toy.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 target_spec.py

Each computation takes seconds in this environment. No target chain
with exponentially many physical sites is allocated. In particular,
target_spec.py's large b values are only integer inputs to its formula.

clock_results.json records the following:

- All 6,069 rule/context combinations over the 17 clock/port letters
  preserve the complete local defect pattern.
- Full open-segment graphs at n=1,2,3,4 contain respectively
  15,225,3375,50625 words. Their only legal components have
  9,30,79,192 vertices. At n=4 there are 16,477 components, including
  3,855 cyclic illegal components. No cycles were discarded.
- Every illegal-component spectrum was diagonalized through n=3.
  With kappa=8,J=1 their minimum energies are 8,
  5.381966011250105 and 4.763932022500208. All obey the analytic bound.
  Illegal spectra at n=4 were not diagonalized.
- Full port/clock rings at L_ring=3,4 contain 4,913 and 83,521 words.
  Every transition conserves every violation; all legal marker-free
  words have no head. These tests enumerate classical clock words,
  not the tensor-product quantum work Hilbert space.
- Selection and both exact midpoint margins are checked for selected
  n=2,...,40 against all candidate n=1,...,100. The proof, not this
  finite cutoff, excludes still longer cells.
- Endpoint-pin eigenvalues are checked through C=5124. At the largest
  C, the measured lift is 7.33301e-11 and C^3 times the lift is
  9.86528, approaching pi^2. The first run exposed loss of significance
  from direct cosine subtraction in the clock gap at n=40; the final
  code uses a sine product. This was a numerical evaluation defect,
  not a vanishing analytic gap.

Selected midpoint values, for J=1, are

| n | b | C | gamma_len | Delta_density | omega |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 4 | 30 | 3.94531e-2 | 7.89062e-3 | 3.06788e-2 |
| 4 | 6 | 192 | 5.30376e-4 | 7.57680e-5 | 7.94801e-4 |
| 8 | 10 | 5124 | 5.89491e-7 | 5.35900e-8 | 1.12728e-6 |
| 12 | 14 | 114696 | 1.07140e-9 | 7.14266e-11 | 2.25070e-9 |

interface_results.json records a 288-dimensional two-port bare-check
example. Its resolvent nonscalar residual is below 2.7e-15 and its
second-order nonscalar residual is 1.10e-14, agreeing with Section 7.
It also records the three-time positive-frame failure and a separate
40-dimensional richer-check test. In the latter, a nontrivial target
resolvent is approached with errors 8.20e-5,8.34e-6,8.36e-7 at
J=100,1000,10000; all are below the proved Schur bounds.

routeA_toy_results.json contains the complete legal-cell block spectra
summary and the independent finite-lambda checks in Section 8.
target_results.json contains r=1,2,3 target checks at two-cell logical
dimensions 4,16,64, with Jordan-Wigner residuals at most 3.62e-16.
The abstract dimension bounds are 76,220,584. None is a microscopic
site dimension.

## 11. Requirement status and ranked next steps

[Proved] Task 1 is met for the explicit clock with a fixed arbitrary
unitary work decoration and correctly supplied product initialization
checks. There is one complete legal clock path at every n, all
illegal clock components pay conserved penalties, and all rejected
product inputs have a finite hopping gap. This includes cycles and
truncated histories. It does not assert that an unspecified route-A
interpreter has been initialized or compiled.

[Partial] Task 2 is met by the literal small two-state example and by
a target-only formula depending on b. The requested unbounded route-A
history program is not supplied. C19 now gives the stronger abstract
uniform target; the remaining task is its spatial energy compiler.

[Proved and numerical] Task 3 is met for the actual K(b) of this clock,
including the exact midpoint excesses and input-gap distinction.
The generic substitution of an arbitrary growing K(b) is false
without an additional convexity check.

[Partial, exact counts for completed objects] Task 4 is met for the
clock wrapper and the d0=64 example, with exhaustive small clock tests
and finite-cell spectral tests. It is not met for unbounded c because
w for a complete program is not known.

The next steps, in order, are:

1. Compile C19's additive two-end target using a richer boundary menu
   that annihilates the code. The Section 8 term demonstrates a
   two-site interface with a partner and an auxiliary factor on one
   bulk site. A useful next deliverable is a finite menu spanning the
   required syndrome operators, with every work factor counted.
2. Implement that menu in C8's sweep-machine interface. Supply the
   initialization and computed logical-region checks, arithmetic
   precision and runtime bounds, and weighted event schedule. Retain
   the fixed clock graph so those choices cannot change its length.
3. Prove the finite-cell source-resolvent error and compare all errors
   with both gamma_len and the smaller input gap. Approximate codes
   also need their first-order leakage below the exchange, not merely
   below the clock gap.
4. Only after those steps give the microscopic d0. Mobile-cell
   transport and finite-coupling infrared control are subsequent
   problems and were not treated here. The frozen construction still
   has translated ground sectors on compatible rings.

The merged clock is a viable selector. The remaining obstruction is
now a specific local energy-program and alphabet problem, and the
bare static Bell check does not provide an exact-code shortcut.
