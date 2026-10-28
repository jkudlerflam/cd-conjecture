# Long-cell matching: a false two-sided bound and an exact exclusion

Status: **[proved analytic statements, with independent finite numerical checks]**.
The literal two-sided matching bound in the infrared referee response cannot
hold with `k` equal to the number of wrong cells or interfaces when wrong-cell
lengths are unbounded. Its failure has the **favorable sign**: the selected
crystal gains an extensive negative charge-relaxation energy. This does not
refute the one-sided defect exclusion inequality, the length-weighted treatment
already suggested in the infrared report, or the proposed c=24 phase.

A separate nonperturbative bound excludes a legal one-marker cell occupying an
entire sufficiently large ring, without a small-J expansion or a lower bound on
its vanishing walker level spacing. General islands of wrong cells and arbitrary
illegal patterns are not settled by these statements.

Author: GPT-6-family Codex matching-attack subagent, 2026-10-28. The exact
commutator was derived analytically and independently checked by sparse ED.
The parallel selector-numerical worker independently found the same favorable
extensive residual in the actual lambda=0 small-b quantum-selector model.
No human review is claimed.

## 1. Notation and the exact assertion being tested

There are N occupied material sites and physical ring length L=2N. Write

    H = H_gap + J sum_j A_j h_j,       A_j=1[g_j=0].

Let E_s^23(N) be the exact frozen selected-sector energy, and E_w^23(N) the
frozen energy of another material sector. Let alpha_L be the free-sea contact
probability. Define the charge-dressing residual in sector a by

    R_a(J,N) = E_a(J,N)-E_gap(N)-J alpha_L E_a^23(N).

The literal second part of `infrared/REFEREE_RESPONSE.md` asks for an absolute
bound of the form

    |R_w-R_s| <= J k epsilon(J/t) + boundary rounding,

uniformly in allowed defect patterns. If k counts wrong cells, a ring with one
cell of length N has k=1 for every N. The argument below proves that its left
side instead grows at least as const*J^2*N/t-O(J)-O(t/N), at every sufficiently
small fixed J>0. Subtracting an identical additional bulk term from both
sectors does not change R_w-R_s.

This is a claim about a two-sided *remainder estimate*. The sufficient sector
inequality `E_w-E_s >= J alpha gamma k - J k epsilon` is one-sided. Our growing
correction has positive sign and therefore helps that inequality.

## 2. A nonperturbative extensive contact-relaxation bound

The total material contact operator has the exact physical-coordinate form

    C=sum_j A_j=sum_x n_x n_(x+1).

Thus the uniform scalar contact problem is the canonical hard-core hopping
chain

    H_c(v)=t T+v C,
    T=-sum_x (b_x^dagger b_(x+1)+h.c.).

Let |0> be its free half-filled ground state and define the local Hermitian
operator G=i[T,C]. Put

    B_L=<0|[C,[T,C]]|0>
       =2 <0|C (T-E_0/t) C|0> >=0.

Use the normalized trial state exp(i theta G)|0>. Its energy derivative at
zero is -v B_L. Each summand of G has support on four sites and norm at most
one: it is the hopping current on a bond multiplied by
`n_(x-1)-n_(x+2)`. A two-site summand of T or C overlaps at most five such
summands. The resulting commutator has support on at most five consecutive
sites and overlaps at most eight further summands of G; using nine is a safe
upper bound. Therefore, uniformly in theta and system size,

    |d^2 E(theta)/d theta^2| <= 180 L (t+|v|).

Taylor's theorem and optimization over theta give the exact variational bound

    E_c(v) <= E_c(0)+v <C>
                -v^2 B_L^2 / [360 L(t+|v|)].                 (1)

This is not a perturbation series with a size-dependent radius of convergence.

The coefficient is strictly extensive. Wick's theorem gives

    lim_(L->infinity) B_L/L = 1/pi+16/(9 pi^3)
                            =0.37564594739836754... .        (2)

For completeness, at half filling `f_r=<c_0^dagger c_r>` has
`f_0=1/2`, `f_1=1/pi`, `f_2=0`, and `f_3=-1/(3pi)`. A local double commutator
is

    (n_(x-1)-n_(x+2))^2 (c_x^dagger c_(x+1)+h.c.).

The three-bilinear Wick determinant is

    <n_(-1) n_2 c_0^dagger c_1>=1/(4pi)-4/(9pi^3),

which immediately gives (2). Finite rings use
`f_r=L^-1 sum_(filled k) exp(i k r)` and have the same positive limit.
In particular, B_L/L>=1/4 for all sufficiently large L. For |v|<=t, (1)
then lowers the energy by at least

    v^2 L/(11520 t) = v^2 N/(5760 t).                       (3)

## 3. Apply the scalar trial to the exact selected internal state

Take an exact selected-sector ground state of H_23, choosing it also to
diagonalize cyclic material translation in the union of translated selected
sectors, and a translation-invariant charge trial state for (1). The product
expectation is exactly

    <H>=<H_gap> + J (<C>/N) E_s^23(N),

regardless of the detailed periodic profile or noncommutation of the h_j.
Set

    e_s(N)=E_s^23(N)/N,             v=J e_s(N).

The original empty-data product-code trial has zero port expectation, so
`e_s(N)<=e_b<0`, including nonzero lambda and g. All internal parameters are
fixed, hence |e_s(N)| is uniformly bounded by a finite interaction norm M.
For J M<=t, equations (1)-(3) imply

    R_s(J,N) <= -J^2 e_b^2 N/(5760 t) + O(t/N).             (4)

The O(t/N) is the possible cyclic-root/charge-flux gluing cost. It cannot be
discarded on a finite ring, but does not affect this argument. One can spread
a bounded twist over all physical bonds; a real zero-current charge trial
changes its hopping expectation by at most O(t/L). The uniform contact term
is unchanged. Equivalently, compatible internal cyclic momentum is supplied
by the usual charge-sea twist. This argument uses no finite-size charge gap.

## 4. One giant cell has only a bounded internal energy

A legal ring with exactly one marker has exactly one token. There cannot be
two simultaneous endpoint tokens, so every port term annihilates this sector.
The nonnegative marker, data, express, excitation and illegal-pair diagonal
terms remain nonnegative after multiplying their bond contributions by A_j.
The remaining internal hopping/conversion matrix has a row-sum bound
independent of N: a conservative bound is

    D=3+|t_d|+g_*.

The gap-coordinate projectors have norm one and do not increase this bound.
Consequently the full doped one-cell sector obeys the exact operator lower
bound

    H_giant >= H_gap-J D,
    E_giant(J,N) >= E_gap(N)-J D.                           (5)

This deliberately avoids claiming that the excited-token gap remains uniform
under charge fluctuations. Its negative kinetic norm alone suffices.

The frozen giant-cell energy is exactly

    E_giant^23(N)=mu_sel-2 cos[pi/(N+1)],

and is bounded independently of N. Thus (5) yields `R_giant>=-J D'` for a
fixed finite D'. Combined with (4),

    R_giant-R_s >= J^2 e_b^2 N/(5760 t)-J D'-O(t/N).         (6)

This diverges with N at fixed J, contradicting the literal absolute per-cell
bound with k=1. The growth is not an infrared instability: it is the ordinary
extensive relaxation energy of the selected bulk which is missing from the
one-token giant cell.

There is also a useful exact one-sided consequence. The selected product-code
and free-charge trial already give

    E_s(J,N) <= E_gap(N)+J alpha_L N e_b+O(t/N).

Hence

    E_giant-E_s >= J[alpha_L |e_b|N-D]-O(t/N).              (7)

For every fixed J>0 this is strictly positive at all sufficiently large N.
For the candidate, `e_b=-0.00126401272424...` and alpha tends to
`0.2973576327...`. The sufficient length from this deliberately coarse bound
is much larger than the sharper infrared cutoff 73. It is nevertheless a
nonperturbative exclusion of an actual unbounded-cell sector, not an appeal to
its shrinking gap.

## 5. The correct weight and what remains

An all-pattern two-sided matching estimate can use a weight such as

    W = number of short wrong cells + number of illegal interfaces
        + total material length in long wrong cells/fragments.

Alternatively, measure a giant replacement by the number of selected cells
removed. The long-cell bare margin is linear in precisely this weight, so the
repair is compatible with sector exclusion. If only a one-sided lower bound
is sought, the favorable extensive term in (6) may instead be retained.

Equation (7) does not prove the corresponding statement for a giant-cell
island surrounded by an infinite selected fluid. In that problem the common
selected bulk and its interface energy still have to be compared. Nor does
it settle arbitrary collections of short defects or illegal fragments.
Those are the remaining matching problems. No claim of a verified c=24 limit
follows from this folder.

## 6. Exact finite-cell cancellation relevant to that remaining problem

At lambda=0, split every microscopic bond into its contributions to individual
cells, retaining both distributed onsite halves. A cell c of length ell has
contributions h_(c,j) on a set of at most ell+1 consecutive charge bonds, with

    sum_j h_(c,j)=H_c.

A boundary charge bond may carry contributions from both neighboring cells;
these act on separate internal cell spaces. Let m be the number of charge
bonds used, define

    f_j=h_(c,j)-H_c/m,
    B_j=sum_(r<=j) f_r,        B_0=B_m=0.

Discrete summation by parts gives the exact operator identity

    sum_j A_j h_(c,j)
       = Abar_c H_c + sum_(j=1)^(m-1) (A_j-A_(j+1)) B_j.   (8)

The first term is diagonal in the isolated cell energy basis. Every
nontrivial cell transition therefore couples to a discrete contact
difference. Equivalently, for different eigenstates e,g,

    sum_j <e|h_(c,j)|g>=0,
    <e|sum_j exp(i q j)h_(c,j)|g>=O(q).

The smooth density part of A has dimension one. Its off-diagonal machinery
coupling gains a spatial derivative and has local spectral density
proportional to omega^3, instead of a generic Ohmic omega term. Dimension-two
contact descendants are also super-Ohmic; the first harmonic has dimension
four. The exact q=0 microscopic contact response need not vanish: (8) is not
such a claim.

This cancellation is stronger than the generic noncommuting spin-impurity
toy. A vacuum excursion must leave and reenter the vacuum through the
contact-difference terms, so it has at least two gradient endpoints. The
parallel uniform-matching worker independently derived this endpoint
statement and is examining its all-order implications. Equation (8) alone
does not establish convergence of the required microscopic linked expansion.

## Reproduction

Run from the repository root:

    PYTHONPATH=/private/tmp/cd-campaign-python \
      /Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
      contrib/counterexample_audit_20261028/matching_attack/check_uniform_contact.py

The script uses numpy/scipy, runs in about 0.2 seconds on this host, and writes
`uniform_contact_results.json`. It builds hard-core hopping independently in
the fixed-particle physical basis, computes the commutator spectral moment,
checks the finite-ring Wick determinant, and checks (1) for v=+/-0.02,+/-0.04.
For L=6,8,10,12,14,16, B_L/L decreases from 0.4444444 to 0.3827658 toward (2).
All 24 variational-bound comparisons pass; the largest Wick discrepancy is
3.2e-14. These tests validate the analytic contact lemma, not the full d=24
phase.
