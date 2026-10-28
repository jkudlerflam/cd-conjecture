# Exact quantum-selector stress test

Status: **[numerical] no wrong-cell sector undercuts the selected sector in
264 exact diagonalizations, through physical length 20.** The calculations
keep the complete vacancy bath and every vacuum-walker position in every
cell. They expose a growing, favorable giant-cell remainder and a substantial
quantum correction absent from a fixed scalar contact source. They do not
establish a bound uniform in length, prove the n=23 matching hypothesis, or
test nonzero port coupling.

Independent implementation and analysis by a GPT-6-family Codex subagent,
2026-10-28, coordinated with the `uniform_matching_proof` and `matching_attack`
workers. AI-assisted; no human review. Only this folder was written.

## Exact model and scope

The five occupied letters are precisely

    L0, R0, Xd_v0, M_empty, M_v.

Together with the physical vacancy they give six local states. This is an
invariant restriction of the literal h23 bond at lambda=0 and data bits zero.
No replacement impurity or XX payload is introduced. The occupied kinetic
moves have amplitude -1:

    (Xd_v0,R0) <-> (L0,Xd_v0),
    (Xd_v0,M_empty) <-> (L0,M_v).

The marker potential is distributed as mu_sel/2 at each endpoint of an
occupied-occupied pair, including bonds between cells. In particular, neither
marker onsite term is moved outside its contact projector. All material words
in a tested sector obey the original arrow grammar, so their illegal-pair
penalties vanish. The physical occupied-number potential is constant at fixed
N and is omitted.

For a cyclic cell-length pattern ell=(ell_1,...,ell_m), let N=sum ell_i and
L=2N. Cell i has an ell_i-state walker, with every position retained. Write
g_j for the vacancies after material site j, A_j=1[g_j=0], and H_gap for the
unit-amplitude phase-model hopping with sum g_j=L-N. The exact restriction is

    H = H_gap
        - J sum_cells sum_(r=0)^(ell-2)
            A_(start+r) (|r+1><r| + |r><r+1|)
        + (J mu_sel/2) sum_markers (A_(marker-1)+A_marker).

Charge and walker operators are diagonalized together. They are not replaced
by their expectations. Selected length b=2 or 3 uses the original midpoint
formula

    a_l=2[1-cos(pi/(l+1))],
    delta_l=(l+1)a_l-l a_(l+1),
    mu_sel=2-(delta_b+delta_(b-1))/2.

The frozen energy is E_frozen(ell)=sum_i[mu_sel-2cos(pi/(ell_i+1))]. Both
midpoints select their intended length among the tested legal sectors.
For b=2, mu_sel is negative and the selected path has only one kinetic edge;
its internal bond operators commute. This is a useful calibration, while
**b=3 supplies the noncommuting quantum-selector test** and has positive
mu_sel. These toy parameters replace the n=23 value b=24; no conclusion is
transferred uniformly in b.

## Physical ring, rather than a rooted covering

The code implements the physical momentum-zero ring quotient. This is enough
for the exact ground energy of each tested sector: its offdiagonal matrix
elements are nonpositive, so a Perron ground vector can be chosen positive
and invariant under physical translation.

A root is placed at the start of a cell-length word. Physical translation
has phase one in this subspace. If the length word repeats r times, there
remain r equivalent marker roots: rotate the gaps by one primitive word and
rotate the walker coordinates by the corresponding number of cells. The code
projects onto normalized orbits of this simultaneous rotation, retaining
shorter orbits and their square-root normalization factors. A nonrepeating
word has a unique distinguished marker root. Thus neither an extra root
position nor unphysical cyclic sectors are counted.

An independent direct construction enumerates physical six-letter words and
acts with the local physical bond around the ring. The direct and quotient
ground energies agree in five checks, including repeated and nonrepeated
patterns, with maximum error 1.25e-14. The quotient dimension need not equal
the covering dimension divided by r because of stabilizers; for (2,2) at
L=8 it is 73 rather than 70.

`validate()` additionally compares the independently written 25-by-25 occupied
bond with the corresponding restriction of `explicit_candidate/bond23.py`:
entry error and leakage are both zero. That source is imported only for this
comparison. The Hamiltonian construction imports no earlier squeezed-model,
selector, scalar-impurity, or quantum-impurity implementation. Its coordinate
convention was checked against the ring-quotient derivation in the vacancy
audit; this provenance is conceptual rather than reused code.

## Baseline and measured defect correction

The exact finite-ring charge ground energy and mean contact are

    E_charge = -2 sin(pi N/L)/sin(pi/L),
    alpha_L = N/L - sin(pi N/L)^2/[LN sin(pi/L)^2].

The program independently solves the charge bath and verifies its energy.
It also verifies that the product of its ground state with the frozen sine
walkers has first-order energy alpha_L E_frozen in every sector.

For every wrong pattern, the reported common-charge-subtracted remainder is

    R = E_exact(wrong)-E_exact(selected)
        - J alpha_L [E_frozen(wrong)-E_frozen(selected)].

This definition also subtracts the exact selected-sector energy, so its bulk
quantum correction is retained in the comparison. It does not subtract a
separately fitted bulk energy for each wrong phase.

The small suite exhausts **all cyclic positive-integer cell-length
compositions** for (b,N)=(2,4),(2,6),(3,6), 31 sector cases. Additional targeted
patterns use (b,N)=(2,8),(3,9),(2,10), including neighboring and separated
wrong cells and an entire ring occupied by one long cell. Every case is run
at J=0.1,0.05,0.02,0.01,0.005,0.002, with vacancy hopping t=1.

There are 44 sector cases and 264 quantum ground-state diagonalizations.
The largest quotient has 2,217,072 states. All wrong-sector gaps are positive;
the smallest measured gap is 0.000123222677 at b=3,N=6,J=0.002 for (2,4).
The largest ground-eigenpair residual is 2.15e-11. These are floating-point
numerical checks, not interval-certified eigenvalue enclosures.

Local comparisons at J=0.002:

| Selected b | N | Wrong pattern | R/J^2 |
|---:|---:|:---|---:|
| 2 | 6 | (1,2,3) | 0.03348750 |
| 2 | 8 | (1,2,2,3) | 0.03238342 |
| 2 | 10 | (1,3,2,2,2) | 0.03142824 |
| 2 | 8 | (1,2,3,2) | 0.03213141 |
| 2 | 10 | (1,2,3,2,2) | 0.03145273 |
| 3 | 6 | (2,4) | 0.01843732 |
| 3 | 9 | (2,3,4) | 0.01827743 |
| 3 | 9 | (6,3) | 0.03435480 |

The sampled local coefficients show no growth with system length or defect
separation. They do not prove a uniform limit. Both signs of R occur in the
exhaustive suite: b=3,N=6,pattern (1,1,4) has R/J^2=-0.00391070 at J=0.1 and
-0.00344767 at J=0.002. Its much larger positive order-J bare gap remains.
Thus these data do not suggest that every dressing correction has a
favorable sign.

## A giant cell is not a defect of bounded spatial support

For a single cell of length N, k_wrong=1 if k counts wrong cells rather than
replaced material length. At J=0.002:

| b | N | R/J^2, one giant cell |
|---:|---:|---:|
| 2 | 4 | 0.02467550 |
| 2 | 6 | 0.07333795 |
| 2 | 8 | 0.13117362 |
| 2 | 10 | 0.19388207 |
| 3 | 6 | 0.03573171 |
| 3 | 9 | 0.11858746 |

This growth is consistent with removal of the selected phase's extensive
quantum/contact relaxation. It favors the selected sector and accompanies
an order-J bare cost growing with N. Finite data alone cannot refute an
arbitrary uniform bound, but they support the companion analytic warning:
an **absolute** O(J^2 k_wrong) remainder is the wrong target if a whole
macroscopic wrong cell counts as one bounded defect. Long cells should use
their extensive selector margin, a one-sided comparison, or a defect measure
that counts displaced material length. These data provide no adversarial
undercutting from long cells.

## The quantum effect is substantial

A separate variational calculation freezes every walker to its sine ground
state, replacing its kinetic matrix on each contact by a scalar expectation,
while still diagonalizing the full charge bath. The exact quantum energy
must lie below this constrained energy. The check finds

| b | N | Pattern | (E_quantum-E_fixed_walker)/J^2 at J=0.002 |
|---:|---:|:---|---:|
| 2 | 8 | (2,2,2,2) | 0 within numerical error |
| 3 | 6 | (3,3) | -0.09553550 |
| 3 | 9 | (3,3,3) | -0.14396964 |
| 3 | 9 | (2,3,4) | -0.12568076 |
| 3 | 9 | (9) | -0.13898225 |

The b=2 equality follows from its single commuting edge. The b=3 results
show quantum relaxation consistent with extensivity, missing from a scalar
source test. For b=3 the independently checked
transition form factors also satisfy sum_j <m|h_j|0>=0 to 1.7e-16, as expected
from the cell eigenvalue equation. This sum rule is not by itself a finite-J
uniform defect proof.

## Reproduction

From the repository root, using numpy and scipy:

```sh
PYTHONPATH=/private/tmp/cd-campaign-python \
/Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 -B \
contrib/counterexample_audit_20261028/selector_numerics/quantum_selector.py --suite small
```

Repeat with `--suite large` and `--suite extended`, then run
`contrib/counterexample_audit_20261028/selector_numerics/summarize.py` using the
same interpreter. The suites take approximately 0.5, 19, and 45 seconds on
this host; the summary's scalar comparison takes about one second. Each suite
runs the independent physical-ring and literal-bond validation first.
`--suite validation` runs just those checks.

The three `*_results.json` files retain every sector, energy, residual,
first-order check and defect subtraction. `summary.json` includes the
fixed-walker comparison. Run logs retain solver timings. No CFT central
charge is measured, no excited charged/express sector is included, and the
nonzero-lambda payload hierarchy remains outside this diagnostic.
