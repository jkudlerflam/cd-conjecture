# Independent audit of the all-string positive completion

Status: **accepted as a finite local construction and full rooted-space
reduction, with the strengthened heat-trace constant specified below**.
This is not acceptance of the remaining massive influence-action or
multichannel infrared theorem.

Author: OpenAI GPT-6-family Codex agent `exact_escape`, 2026-10-28.
Independent from the author of the sibling completion proposal. No human
review. Sources checked directly: C8 `selector_wrapper.md`,
`arrow_wrapper.md`, `helix_equalizer.md`; C32 `PROOFS.md` Section 1;
C39 `PROOFS.md` Sections 1 and 5; C44 `PROOFS.md` Sections 1-3 and 6;
C44's literal 46-letter rules in `core_literal.py`; and C45 `PROOF.md`
Sections 2-3. All claims below concern those stated rules or the explicitly
defined positive modification, not an unspecified decoder.

## 1. The positive modification is legitimate and useful

Use selector arrows

    # : R -> L,   L : L -> L,   X : L -> R,   R : R -> R.

The old endpoint-helix pin term on a physical bond is

    (|X#><X#|+|R#><R#|) tensor (a_B,left+a_A,right),

with the positive one-qubit pin a=I-|+><+|, times its positive scale J_h.
Replace its selector projector by I_left tensor |#><#|_right. The
difference is exactly

    (|L#><L#|+|##><##|) tensor (a_B,left+a_A,right)>=0.     (1)

The pairs L# and ## are forbidden, while X# and R# are the only legal
pairs entering a marker. Therefore (1) vanishes on the entire legal
selector sector. It does not merely vanish after code projection. The
head and marker rewrites preserve selector violations, and ports/readout
do not change selector letters, so all dynamics within that legal sector
are unchanged at every coupling.

On an arbitrary rooted cell, the new A-track pin acts at its first
marker independently of the previous cell. The B-track pin acts at its
last site independently of the next cell's internal state. Thus these
terms are truly on-cell in rooted coordinates. Under the old conditional
pinning, a first-site pin could depend on the previous cell's final
selector letter, which was an intercell dependence in malformed sectors.

At consecutive markers a length-one cell receives both pins on its sole
site: A from the incoming bond and B from the outgoing bond. There is no
missing pin or duplicate energy allocation. The two pins act on distinct
qubit tracks. The ## selector penalty belongs to the cell on its left.

The change is a positive nearest-neighbor term on the same 4464-state
alphabet. It defines a slightly modified candidate bond. Its legal-sector
code, source matrix elements, second-order retained Hamiltonian, and
finite-coupling legal-sector dynamics agree exactly with the old ones.

## 2. A lower bound on each illegal cell

Let ell be the rooted cell length, counting its first # and excluding
the next #. Include the closing bond's selector penalty in this cell.
Let z be the number of forbidden selector pairs among its internal
bonds and that closing bond. Each cell contains exactly one marker, so
its marker chemical energy is mu_sel, not a fraction of it.

Suppose z>=1. Remove the internal forbidden bonds. There are at most
z+1 nonempty legal fragments. A fragment has no interior marker and at
most one X head. Its head hopping is the adjacency operator of a path,
with lower bound -2t_sel, including a fragment of one site. Its wrapper
Hamiltonian is an open wrapper chain with lower bound -4U from C8/C44.
All helix and new pin terms are positive and may be dropped. Hence, with

    C=2t_sel+4U,       K_sel=4t_sel+9U,

the full cell satisfies

    H_cell >= mu_sel+K_sel z-C(z+1)
           = mu_sel+U+(z-1)(2t_sel+5U)
           >=mu_sel+U.                                  (2)

If the closing bond is forbidden, it contributes to z without making
an additional internal fragment, so z+1 is still a valid upper bound.
The formula remains valid when forbidden bonds overlap.

For ell=1, the only selector pair is ##. There is no internal head or
wrapper bond; the exact lower bound is mu_sel+K_sel. This is stronger
than (2). Its positive pin terms are again harmless.

After subtracting the coexistence energy e ell, with e<0,

    A_cell >=mu_sel+U+|e|ell                               (3)

on every selector-illegal cell. This is an on-cell operator inequality
on the full wrapper/core Hilbert space. It uses the complete open-wrapper
bound, including illegal wrapper arrows and arbitrary core tokens.

On selector-legal cells of lengths outside b,b+1, the existing C44 bound
is

    A_cell>=gamma(1+|ell-b|),   gamma=t_sel eta_b-U>0.     (4)

At the two selected lengths, the recorded complete code and positive
excitation gaps are unchanged by (1).

**Necessary constant adjustment.** Inequality (3) need not dominate
(4) with the old gamma, because |e| and gamma are different parameters.
A uniform all-string replacement is

    gamma_all=min(gamma, U/b, |e|)>0.                     (5)

For ell<=b, 1+|ell-b|<=b, so U dominates the desired bound. For ell>=b,
1+|ell-b|<=ell, so |e|ell dominates it. Therefore every wrong-length
cell, including selector-illegal ones, satisfies

    A_cell>=gamma_all(1+|ell-b|).                         (6)

One must replace gamma by gamma_all in an all-string heat-trace estimate,
or retain separate legal and illegal geometric tails. The earlier C46
constant cannot simply be quoted unchanged. The all-string on-cell gap
can be chosen as the minimum of gamma_all and the complete selected-cell
gaps; it remains positive at each fixed n and fixed microscopic schedule.

The sibling's finalized note uses an equivalent, more conservative
constant

    kappa=min(Delta,gamma,U,|e|)/(b+2),
    A_ell>=kappa(ell+1)Q_ell.                             (6a)

This is valid: at selected lengths ell+1<=b+2; for wrong legal lengths,
1+|ell-b|>=(ell+1)/(b+2); and (3) bounds illegal lengths. It yields

    Tr_Q exp(-beta A)
      <=d exp(-2beta kappa)/(1-d exp(-beta kappa)),
    beta kappa>log d.                                   (6b)

I independently checked this geometric sum and accept the sibling's
choice as satisfying the needed constant adjustment.

## 3. The unbounded readout has an invariant positive location

Define N_token on a rooted cell by counting its occupied core letters
of types d, x, A_s, B_s with nonempty token label s, including the vacuum
token v and both excited buses. Bare L_a,R_a,A_empty,B_empty and all
wrapper holes count zero. This is a count on the complete alphabet,
including malformed words.

The literal core rules conserve this count:

* A rail hop moves one token between neighboring occupied sites.
* The A/data and express/B links move a token without creating one.
* A data/express junction keeps the same token.
* A charged/neutral conversion changes its species and data bit only.
* Diagonal penalties and endpoint masses do not change any token.

Wrapper hopping moves a complete occupied label past a hole. Ports
exchange v with an excited label at the same endpoint. Crucially, the
marker move transfers H_0 to H_2, or the reverse, and never transfers
an occupied label. Thus every rooted cell's N_token is preserved by
the base, ports, marker motion, and readout. The root-seam move has the
same property in the ordered-cell convention. This direct count is
enough; no claim about every possible malformed grammar flag is needed.

The readout R_j is a sum of diagonal projectors onto a vacuum token at
a data site with the counted data bit. Thus, on the full cell space,

    0<=R_j<=N_token,j.                                    (7)

Let Q_bad,j=1_(N_token,j>=2), and split

    R_j=R_good,j+R_bad,j,
    R_bad,j=Q_bad,j R_j,
    0<=R_good,j<=(1-Q_bad,j),  ||R_good,j||<=1.             (8)

The coefficient c of this readout is positive in C39/C42/C44:
c=(lambda^2/s)J0 mu S/q^2, with every displayed factor positive.
Include c R_bad,j in A_j exactly and retain c R_good,j in the bounded
perturbation. The new on-cell operator is still nonnegative and has
the same code. It retains (3), (6), and the selected gap. This remains
true if R_bad does not commute with the rest of A: positivity and the
min-max principle suffice. Each fixed-length block is finite dimensional,
so the same bounds control its heat trace.

## 4. Compatibility with the local Schrieffer-Wolff construction

The complete code has exactly one token per cell. A local all-code
projector P_X therefore has no support on Q_bad,j for j in X. Every
bounded perturbation and the on-cell resolvent preserve the individual
token counts. Consequently each off-diagonal local generator built from

    Q_X A_X^(-1) V_X P_X - its adjoint

has zero action on a support containing any bad-token cell. The same is
true of the higher generators in C45, because their nested bounded
commutators preserve those token sectors. Thus

    [S_X, R_bad,j]=0                                    (9)

for j inside X, by two-sided annihilation, and for j outside X by disjoint
support. The unbounded positive readout is untouched by the full local
unitary and never enters a bounded-interaction norm estimate.

All remaining readout pieces have norm at most their fixed coefficient
per cell. Ports and marker transfers have fixed finite incoming/outgoing
degree and bounded compensators; marker chemical adjustments are bounded
per cell. Therefore the decomposition A+lambda V+lambda^2 F can be made
on the full all-string rooted space of the completed bond, with bounded
V,F and a positive uniform cell gap.

The dependence of A on lambda through the absorbed positive term causes
no denominator problem: generators connected to the code lie in the
one-token sectors where that term vanishes. Uniform lower bounds on A
also only improve with the added positive term.

## 5. Verdict and remaining obligations

The completion repairs the explicit on-cell pinning problem and supplies
the full-space positive base needed by the length-fugacity thermal
comparison. The token-count split removes the particular unbounded
readout obstacle without discarding malformed states. I independently
read and reran the sibling's `ensemble/illegal_checks.py`: all 32
selector-neighbor checks, 3280 fragment cases, 2080 marker reanchoring
cases, 176 core pair rewrites, and 12 endpoint-port cases passed. The
symbolic proofs above cover arbitrary length beyond those finite tests.

This acceptance requires using the completed microscopic bond and the
revised gamma_all, or an equivalent split-tail heat-trace estimate. It
does not retroactively assert that the old conditional pins were on-cell
on all malformed strings. It does not prove convergence of the auxiliary
influence action, the marginal RG flow, or a microscopic T^2 coefficient.
Those infrared obligations remain after this finite construction audit.
