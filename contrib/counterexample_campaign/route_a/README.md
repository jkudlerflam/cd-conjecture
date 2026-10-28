# Global sector control for the 23-letter chain

## Result and scope

[proved] The complete frozen 23-letter bond of
[`explicit_candidate`](../../explicit_candidate/README.md) selects its intended
cell length on every compatible ring at a fixed **nonzero** port coupling.
Neither illegal words nor other marker patterns undercut the selected sector.
This fills the global-sector qualification in that candidate's Section 2.3.
The proof compares exact interacting cell chains; it does not compare an
extensive perturbation bound with a finite defect gap.

This is not by itself a counterexample. The C33 and C37 results give a
certified nonresonant 23-flavour quadratic reference for this d=23 bond.
The C33 matrices and certificate were independently rebuilt and audited
here; see [the audit](C33_AUDIT.md). The interacting infrared transfer and the interpretation of frozen
translation order remain separate questions. A vacancy-based mobile
completion is being investigated by the coordinating worker.

Author: Codex GPT-6 family, route-A worker, 2026-10-28. The proof and checks
are AI-assisted and have not been human-refereed. Reproduction:

```
python3 contrib/counterexample_campaign/route_a/check_global_sectors.py
```

## 1. Parameters

[definition] Use the exact alphabet and term list in the linked candidate.
Set t=t_x=t_sw=1, delta=0, and D_x,D_M>=0. Rail and loop phases are arbitrary.
Let

    g_* = max(|g_+|,|g_-|),   kappa = Delta-2-|t_d|-g_* > 0,
    epsilon >= 0,            mu_sel >= 0.

The selected length is b>=2. Put

    E_l = mu_sel - 2 cos(pi/(l+1)),
    e_b = E_b/b < 0,
    gamma_b = min_(l>=1, l!=b) (E_l-l e_b) > 0.

The standard midpoint selector supplies these strict inequalities. Write
B=sqrt(3)|lambda| and assume

    gamma_b > 5 B,       K > 4 + 7 B.                     (1)

No coupling depends on the total ring length L. The theorem below concerns
L=mb with m>=2. It establishes a sector gap, not a spectral gap inside the
selected critical sector.

## 2. Cell and fragment lower bounds

[proved] At lambda=0, every complete legal cell of length l has minimum
energy E_l. In its vacuum sector the walker is a length-l path. The
epsilon term is nonnegative and vanishes on the empty data word, so its
lowest eigenvalue remains E_l.

The excited token graph, before conversions, is the folded path with an
extra edge of weight |t_d|. Its unit-edge adjacency norm is at most 2,
and the extra edge has norm |t_d|. This assertion holds also for truncated
end fragments: deleting vertices/edges cannot increase this absolute-row-sum
bound. Each conversion block consists of disjoint two-state transitions,
so its norm is g_*. Extra express/marker energies and the data field are
nonnegative. Therefore every excited token, before the nonnegative marker
potential is added, has energy at least kappa>0. The vacuum ground energy
before the marker potential is at most zero. The complete cell minimum
is consequently attained in the vacuum sector, for **every** l, including
arbitrarily long competing cells. At l=1 the cell has only a marker and
the assertion is immediate.

[proved] A legal open fragment between a cut and the nearest marker, or
between two cuts when there is no marker, has at most one token. A vacuum
token contributes energy at least -2; an excited token contributes at least
kappa. Thus either end fragment is bounded below by -2, independently of
its length. A legal open interval contains complete cells and at most two
such end fragments. This follows directly from the arrow grammar: a bulk
token changes L to R, and a marker is required to change R back to L.
The junction-switch term has the identity on the preceding marker, so it
does not couple the internal Hilbert spaces of consecutive cells.

[proved] The boundary port operator has norm sqrt(3). In the active
four-species tensor product it is a star linking |v,v> to the three neutral
excited pairs, together with disjoint unit-amplitude exchanges. Hence its
largest singular value is sqrt(3). Restriction to complete cells or
fragments cannot increase this norm.

## 3. A boundary-energy lemma

[proved] Consider the homogeneous chain of selected b-cells, with the full
cell Hamiltonian as its onsite term and lambda V as its bond. Let A_m be
its open-chain ground energy and u the thermodynamic ground energy per
cell. Then

    |A_m - m u| <= B.                                    (2)

Indeed concatenating q copies of an m-cell interval inserts q-1 bonds,
whose total norm is at most (q-1)B. The variational principle, in both
directions, gives

    q A_m-(q-1)B <= A_(qm) <= q A_m+(q-1)B.

Divide by qm and take q to infinity. Existence of u also follows from
these almost-additivity inequalities. Adding the closing bond shows that
the selected periodic-sector ground energy S_m satisfies

    S_m <= m u + 2B.                                     (3)

Moreover u<=E_b: the product of empty-data vacuum code states has zero
port expectation. Thus u/b<=e_b<0. No effective Hamiltonian approximation
has entered (2) or (3).

## 4. Wrong lengths at nonzero lambda

[proved] In a legal ring with k>=1 wrong-length cells, remove every bond
touching such a cell. There are at most 2k removed bonds. The remaining
selected cells form at most k open runs. Their exact ground energies are
bounded using (2). Wrong cells contribute at least E_l>=l e_b+gamma_b.
Since u/b<=e_b, the ground energy in this entire marker sector obeys

    E_wrong >= (u/b)L + k gamma_b - 3k B.

Subtracting (3) for the compatible selected ring gives

    E_wrong-S_m >= k(gamma_b-3B)-2B
                 >= k(gamma_b-5B) > 0.                  (4)

This also covers a ring with no selected cells, with a looser bond count.
The marker-free legal sectors consist entirely of L or entirely of R.
Their energies are nonnegative and cannot undercut the selected sector,
whose product-code trial energy is L e_b<0.

## 5. Illegal words at nonzero lambda

[proved] Every forbidden-pair projector commutes with the full Hamiltonian.
Every non-diagonal bond term acts only on a matched pair and preserves its
external arrows. A sector with r>=1 forbidden bonds therefore splits
exactly into r legal intervals after these bonds are cut; it pays rK.

Let k be the number of complete wrong-length cells in these intervals.
There are at most 2r end fragments, at most k+r selected-cell runs, and at
most 2k+4r port bonds touching wrong cells or end fragments. Apply (2),
the complete-cell lower bounds, and the -2 bound per end fragment. The
result is

    E_illegal >= (u/b)L + k(gamma_b-3B) + r(K-4-5B).

Here replacing the end-fragment baseline zero by (u/b) times its length
only weakens the lower bound, since u<0. Subtracting (3), and using r>=1,
gives

    E_illegal-S_m >= k(gamma_b-5B)+r(K-4-7B) > 0.         (5)

Equations (4)-(5) establish the claimed global sector exclusion. They do
not assume that the selected many-body sector has a spectral gap.

## 6. The certified n=23 parameter point

[source-reported; independently applicable bounds] C33 uses b=24,
Delta=4-2cos(pi/25), D_x=D_M=1000, g=2^-365 and a loop magnitude
t_d=beta*2^-136 with 0<beta<1. Its data field is positive. Thus kappa>0.
For the midpoint selector, the checking script encloses gamma_24 and the
selector parameters using rational alternating-series bounds on pi and
cosines. It proves, using rational inequalities and sqrt(3)<2, that every

    0 < |lambda| <= 2^-20,       K=10                    (6)

satisfies (1). Infrared control may require a vastly smaller lambda;
(6) is only a sufficient global-sector interval. Reducing lambda preserves
all conclusions here. The finite band certificate itself is not rerun by
this script; it is in `codex/work_C33/root_certificate_exact_n23.json`.

The calculation also verifies the port norm and the complete-cell lower
bound for small legal cells with generic fluxes. Numerical checks are
diagnostics for the proof, not substitutes for the all-length bounds.
