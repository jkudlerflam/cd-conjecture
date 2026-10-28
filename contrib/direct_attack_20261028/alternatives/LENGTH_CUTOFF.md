# Finite-order length-cutoff stabilization for the mobile wrapper

## Result and scope

**[Proved conditional on the literal C44/C45 operator definitions audited below.]** Let \(k=n+2\) be the conserved occupied count of each ordered cell and \(b=k+1\). The retained code has physical cell lengths \(b,b+1\), equivalently one or two holes per cell. Write

\[
H(\lambda)=A+\lambda V+\lambda^2(F_0+\mu T),
\]

where \(T\) is unit marker transport and \(A,V,F_0\) preserve every cell's length. All internal wrapper parameters and \(\mu\) are fixed before expanding in \(\lambda\).

Define the cutoff by **retaining every complete cell sector of length \(\ell\le K\)**, with \(K\ge b+1\), and compress every operator to that space. Do not add a boundary counterterm or change the retained on-site matrices. Put

\[
C(r)=\sum_{j=1}^r\left\lceil\frac j4\right\rceil,
\qquad r\ge1.
\]

Then every projected finite-volume Feshbach coefficient, and every coefficient in a fixed algebraic Kato/Schrieffer–Wolff convention built from the same inputs, is identical before and after this cutoff at every degree

\[
\boxed{p<4C(K-b)\quad\text{in }\lambda.}
\tag{1}
\]

More generally, in a multivariable expansion, a projected coefficient involving \(q\) elementary marker moves is cutoff independent whenever

\[
\boxed{q<2C(K-b).}
\tag{2}
\]

These are sufficient bounds. They give the first **possible** cutoff-dependent order, not a claim that the coefficient there is nonzero. The argument is uniform in the number of cells, but supplies no uniform bound on the sizes or sum of the perturbative coefficients.

The simpler bound, requiring only unit length changes and not the limited hole supply, is \(q<2(K-b)\), or \(p<4(K-b)\). The quadratic improvement in (1) uses exact occupied-count conservation and the fact that both endpoint states have at most two holes per cell.

## Audit of which operators change length

The following facts follow from the displayed microscopic rules, rather than an effective low-energy analogy.

| Operator or ingredient | Cell-length action | Occupied-count action | Source |
|---|---|---|---|
| Positive core, occupation filter, endpoint helices, selector head, penalties | Preserves every ordered length | Preserves every ordered occupied count | [C44 §§1–3](../../../codex/work_C44/PROOFS.md) |
| Wrapper vacancy swaps | Moves a hole within a fixed marker interval | Moves the complete occupied label within that interval | C44 §1: \(H_0L\leftrightarrow LH_1\), \(H_1B\leftrightarrow BH_1\), \(H_1R\leftrightarrow RH_2\) |
| Compensated port \(V\) | Preserves lengths | Acts on occupied endpoint tags; does not turn a hole into an occupied label | C44 §§1,4,6 |
| Conditional readout and marker chemical potential | Preserves lengths | Diagonal in the relevant occupation/marker labels | C44 §§1,6 |
| Marker transport \(T\) | \(\ell_i\mapsto\ell_i+1\), \(\ell_{i+1}\mapsto\ell_{i+1}-1\), or the reverse | Transfers only a typed hole, so each cell's occupied count stays fixed | [C8 (11)](../../../codex/work_C8/selector_wrapper.md), C44 §4 |
| Code projections and unperturbed resolvents | Preserves each joint length sector | Preserves each joint occupied-count sector | [C45 §§2–3](../../../codex/work_C45/PROOF.md) |

In particular, marker transport is the physical \(R\#\leftrightarrow\#L\) move with the accompanying \(H_0\leftrightarrow H_2\) conversion. Its helix rotations change internal states, not lengths or occupied counts. The ordered occupied labels stay in their cells. The selector head motion \(LX\leftrightarrow XR\) moves no marker and therefore changes no interval length.

The unperturbed operator has the direct-sum form

\[
A=\sum_i A_i,\qquad
\mathcal K_i=\bigoplus_\ell\mathcal K_{i,\ell}.
\]

Thus every joint length projector commutes with \(A\). The code projection is itself a direct sum of its length-\(b\) and length-\(b+1\) parts, so it also commutes with length projectors. The global reduced resolvent, its powers, and the local reduced inverses \(Q_X A_X^{-1}Q_X\) used in C45 consequently preserve lengths. They may mix many internal states inside a length sector, which does not affect the argument.

Starting and ending in the selected code fixes the occupied number in every cell to \(k\). All intermediate states reached by these operators retain that value. Hence the only unbounded coordinate relevant to this lemma is the hole count \(h_i=\ell_i-k\ge0\). No assertion that marker motion preserves the full cell energy or the code space is made.

The ordered-cell cover can be used on a ring. Cyclic relabeling and the physical root quotient preserve the condition that every length is at most \(K\), so the cutoff and the coefficient identity descend to that quotient. No assertion that only even port powers occur on odd rings is needed.

## The elementary path lemma

Resolve a product of perturbing operators and length-preserving factors into intermediate joint length sectors. Every nonzero summand is a path of hole-count configurations. An elementary marker factor transfers one unit between adjacent cells; all other factors leave the configuration unchanged.

Let \(B=b+1\) be the maximum selected length. If a projected word starts and ends in selected lengths and reaches length \(B+r\) at some intermediate time, that cell needs at least \(r\) increases before the peak and at least \(r\) decreases afterwards. Therefore the word has at least \(2r\) marker factors.

Equivalently, a projected word with \(q\) marker factors can never visit a cell longer than

\[
B+\lfloor q/2\rfloor.
\tag{3}
\]

For a word with only a selected initial endpoint, the corresponding elementary bound is \(B+q\). The return requirement is essential to the factor of two. It must not be applied to an arbitrary full-space remainder or an unprojected generator.

## A stronger path lemma from the finite supply of holes

**[Proved.]** Suppose both endpoints have \(h_i\in\{1,2\}\). A projected word visiting a cell with \(h_0=2+r\) requires at least \(2C(r)\) marker moves.

To see this, label holes for bookkeeping along the discrete path. Labels are only a counting device, not extra physical degrees of freedom. At most two holes at the peak cell can have started there, so at least \(r\) of its holes came from other cells. A hole starting at graph distance \(d\) must have participated in at least \(d\) marker transfers before reaching cell zero. On a line or ring there are at most two other cells at each positive distance, and each held at most two holes initially. Thus at most four of the contributing holes can have initial distance one, at most four more can have distance two, and so forth. The number of marker moves before the peak is at least

\[
1+1+1+1+2+2+2+2+\cdots
\]

truncated after \(r\) terms, which is \(C(r)\). Reverse the path from its final selected configuration to the same peak. The identical argument costs another \(C(r)\). The two path segments are disjoint in time, proving the claim.

Small rings and open boundaries only reduce the number of available nearby source cells. They therefore strengthen this lower bound. Internal charge motion or an arbitrary length-preserving resolvent cannot evade it: those operations do not transport a hole between ordered cells.

For \(r=4a+t\), \(0\le t<4\),

\[
C(r)=2a(a+1)+t(a+1)
     =\frac{r^2}{8}+\frac r2+\frac{t(4-t)}8.
\tag{4}
\]

The threshold is consequently quadratic in the excess length. In the unrestricted nearest-neighbor hole-transfer graph the bound is attained when enough donor cells are available: start with two holes per cell, collect the nearest holes one at a time, and reverse the path. This establishes sharpness of the **counting bound in that larger graph**, not nonvanishing of the corresponding microscopic wrapper matrix element. Typed-hole availability, internal overlaps, projections, and cancellations may delay the actual coefficient.

## Exact Feshbach coefficient equality

Fix a finite number of cells, with periodic or open boundary conditions and no exchange with an external hole reservoir. Let \(P\) be the selected code and \(Q=1-P\). Assume the stated unperturbed gap \(QAQ\ge\Delta Q\). For \(|z|<\Delta\), define

\[
R_0(z)=Q(QAQ-zQ)^{-1}Q.
\]

The cutoff projection \(\Pi_K\) commutes with \(A,P,Q,R_0(z)\), and \(P\Pi_K=P\). Its restricted unperturbed resolvent is exactly the compression of \(R_0(z)\); deleting entire length sectors does not alter the operator in any retained sector.

Put \(D=\lambda V+\lambda^2(F_0+\mu T)\). The Feshbach operator is

\[
\mathcal F(z)=PDP-zP-PDQ\,[Q(A+D-z)Q]^{-1}QDP.
\tag{5}
\]

Expand the inverse as a formal series in \(D\). Each coefficient is a finite sum of words made from perturbations, \(P,Q\), and length-preserving resolvents. The corresponding cutoff word is obtained by inserting \(\Pi_K\) between the factors. Therefore the difference consists only of length paths visiting some discarded sector, whose first possible length is \(K+1\).

The excess above the selected maximum is

\[
r=(K+1)-(b+1)=K-b.
\]

The stronger path lemma requires at least \(2C(K-b)\) marker factors. Each marker factor has degree two in \(\lambda\), proving (1) and (2). Powers of the resolvent arising from energy derivatives remain length preserving, so the statement also holds coefficientwise as an analytic identity in \(z\) throughout the unperturbed resolvent domain.

This is an exact coefficient identity, not an estimate based on a large but finite cutoff. For a fixed finite physical-length sector, every matrix is finite-dimensional. At fixed cell number, the sectors connected to the selected code have total physical length at most \((b+1)M\), so even their all-length representation contains only finitely many accessible length configurations. These observations make finite-volume perturbative coefficients well defined without a thermodynamic convergence assumption.

## Kato and Schrieffer–Wolff conventions

The same coefficient identity holds for any **fixed** perturbative prescription whose coefficients are built algebraically from the original perturbations, selected projections, and unperturbed resolvents or local reduced inverses. Expand the formula into elementary words before counting marker moves. Intermediate selected projections only split a word into smaller excursions, each of which satisfies the same return bound.

This includes contour-resolvent Kato constructions and the local finite-order SW prescription written in C45 §3. Its length-preserving port generator has order \(\lambda\); its length-changing generator components originate from the marker term at order \(\lambda^2\). Commutators, normalizations, and products cannot create length changes with fewer elementary marker factors.

One must use the same gauge/normalization prescription at every cutoff. Effective Hamiltonians chosen with arbitrary cutoff-dependent unitary rotations inside \(P\) need not have identical matrices, even when their spectra agree. No cutoff independence is claimed for full-space operator norms or for the entire off-diagonal SW generator.

## The cutoff convention matters

For a cutoff retaining \(\ell\le K\), merely touching length \(K\) produces no difference: that entire block was retained unchanged. The first omitted block is \(K+1\).

If instead the symbol \(K\) denotes the **first excluded** length, the same result is written with \(r=K-b-1\). The elementary round-trip bound then reads \(2(K-b-1)\); the stronger bound reads \(2C(K-b-1)\). Likewise, adding a counterterm at the last retained length can make merely touching that boundary matter, and requires its own convention.

For the convention used in this note:

| Cutoff retains through | First excluded excess \(r=K-b\) | Minimum marker factors | First possible degree in \(\lambda\) |
|---|---:|---:|---:|
| \(b+1\) | 1 | 2 | 4 |
| \(b+2\) | 2 | 4 | 8 |
| \(b+3\) | 3 | 6 | 12 |
| \(b+4\) | 4 | 8 | 16 |
| \(b+5\) | 5 | 12 | 24 |
| \(b+8\) | 8 | 24 | 48 |
| \(b+12\) | 12 | 48 | 96 |

## What does not follow

At fixed finite volume the first missing coefficient gives a Taylor-order statement inside that volume's convergence neighborhood. The elementary Neumann/Kato radius can shrink with the cell count because the perturbation's full operator norm is extensive. Neither coefficient stabilization nor the transport count supplies a volume-independent connected-kernel bound.

The [finite-cell massive-matching theorem](../infrared/FINITE_CELL_MATCHING.md) also has cutoff-dependent coefficient-conversion constants. Its dedicated-mode embedding grows with the number of excited states. The normal-monomial norm of an auxiliary vacuum projector is \((1+h^2)^m\). Stabilization through an order growing quadratically with the length cutoff does not control those constants, establish a common analytic radius, or justify interchanging cutoff, perturbative-sum, thermodynamic, and zero-temperature limits.

Thus this lemma supplies exact finite-order data and a sharper target for a length-weighted resummation. A further cutoff-uniform connected estimate, a suitably controlled remainder, or an alternative representation is still needed for all-order matching. No new central-charge claim follows.

## Reproduction and provenance

Run

```sh
python3 contrib/direct_attack_20261028/alternatives/check_length_cutoff.py
```

The standard-library script checks the formula for \(C(r)\), 870 finite-ring transport lower bounds, and 24 explicit abstract round trips attaining the infinite-line counting bound. Output is [length_cutoff_checks.json](length_cutoff_checks.json). Runtime is below one second. These are exact integer checks of the combinatorics, not a microscopic matrix-element simulation.

Written 2026-10-28 by a GPT-6-based Codex agent. The local-rule audit uses the displayed C44/C45 definitions; the hole-transport counting and coefficient-stabilization proof were developed in this run. No independent human review or priority claim is made.
