# A 22 state candidate for unbounded central charge

Date: 2026-11-16. AI assistance: OpenAI Codex, GPT-6 family. No independent
agent or human review has been performed.

## Claim and status

**[New explicit candidate; conditional infrared conclusion.]** There is a
single, explicitly specified nearest-neighbour Hermitian bond with 20
material states and two vacancy states. Its alphabet is independent of a
chosen integer n. The same bond acts on every physical link, with couplings
fixed before the thermodynamic limit. Its proposed low-energy theory has n
payload bosons and one vacancy boson, hence

\[
 d=22,\qquad c_{\rm proposed}=n+1.
\]

Here c is the total usual nonchiral Virasoro central charge, allowing
different positive velocities. A compact boson counts as one. The proposed
family would therefore give unbounded c at fixed d=22.

**What is established here:** the literal bond; exact preservation of its
legal-word space; a new all-length selector with a complete finite-cell
code; the two-cell quadratic coefficient and its all-n nonresonant limiting
band construction. These are finite algebraic statements, supplemented by
the checks below. They are not a theorem about the thermodynamic central
charge of the 22-state Hamiltonian.

**Still assumed:** stability of the selected material pattern at nonzero
port coupling; finite-coupling matching with the gapless vacancy fluid and
a positive dressed defect cost; stability and conformal identification of
the resulting multichannel liquid. In particular, this cheaper construction
uses the vacancy route and inherits its extra matching obligations. It does
not replace the 4464-state wrapper at the level of already completed
many-body estimates.

The reduction is nevertheless concrete: it replaces the large selector and
mobile wrapper by a biased vacuum path, a terminal penalty, and two uniquely
determined vacancy colors. It does not infer unboundedness from a finite
list of large-c numerical estimates.

## Dimension comparison and scope

As of this contribution on 2026-11-16, d=22 is the smallest explicit
**mobile candidate in this project with a written arbitrary-n leading-band
argument**. Its frozen precursor already has d=20. This is a comparison
of the constructions inspected here, not a minimality theorem or a
literature-wide priority claim. The new construction has not been
independently reviewed and has additional vacancy-matching obligations
compared with the more extensively audited large wrapper.

| Dimension | Proposed central charge | Scope of the claim |
| --- | --- | --- |
| 2 | approximately 2 | Existing numerical candidate; no unbounded-c construction |
| 13 | 14 | Existing conditional finite member; no completed arbitrary-n band argument |
| 20 | n for arbitrary n | New frozen precursor; its selected ground pattern breaks translation |
| 22 | n+1 for arbitrary n | New mobile candidate with the finite construction described here |
| 4464 | n+1 for arbitrary n | Earlier mobile wrapper with substantially more inherited audits |

The lower-dimensional finite examples do not establish unbounded c, and
the new d=22 proposal does not prove a bound impossible at d<22. See the
[qubit candidate](../d2_updated_conjecture_20261116/README.md),
[d=13 member](../../codex/work_C88/RESULT.md), and
[large-wrapper review](../referee_C44_C45/README.md).

The preceding [compactness analysis](COMPACTNESS.md) records the exact
finite-cover reduction, a sufficient finite-scale entropy criterion, and
why neither proves boundedness near singular Hamiltonians without another
dynamical input.

## The mechanism

Each material cell contains n payload bits in n+1 occupied physical sites.
An explicit vacuum head moves through the bits. Its biased path has a small
positive lowest energy E_n that decreases geometrically with n. A marker
chemical potential selects one prescribed cell length. Two independent
charged virtual channels then produce the two Gram matrices used in the
existing all-n argument of [C39](../../codex/work_C39/PROOFS.md).

The new ingredients are compatible with a small alphabet:

| Function | States |
| --- | ---: |
| Left and right bit backgrounds L_0,L_1,R_0,R_1 | 4 |
| Vacuum head V | 1 |
| Data and express rails, two buses, two charges | 8 |
| Empty markers E_0,E_1 carrying the displaced last bit | 2 |
| Vacuum terminal marker M_v | 1 |
| Excited terminal markers, two buses, two charges | 4 |
| Vacancy colors | 2 |
| Total | **22** |

Keeping the vacuum head and terminal marker distinct avoids the finite-K
code splitting of the smaller merged alphabets. Each letter has one incoming
and one outgoing arrow. Every move preserves the two outer arrows and acts
between allowed pairs. Thus every forbidden-bond projector commutes with
the entire Hamiltonian. There are **no virtual excursions across the grammar
constraint**, even at finite penalty K.

## An exact selector for arbitrary length

**[Proved here.]** The vacuum path for k data sites has positions 0,...,k and
Hamiltonian

\[
 Q_k=\sum_{j=0}^{k-1}(q|j\rangle-|j+1\rangle)
 (q\langle j|-\langle j+1|)+q^2|k\rangle\langle k|,
 \qquad 0<q\le\tfrac14.
\]

The last term is essential: without it the biased path is frustration free
at every length and cannot select a finite cell. Write E_k=λ_min(Q_k).
Factoring Q_k as an upper-bidiagonal matrix squared gives, uniformly in k,

\[
 (1-q^2)^2q^{2k+2}\le E_k\le q^{2k+2},\qquad
 \frac{E_{k+1}}{E_k}\le\frac{q^2}{(1-q^2)^2}\le\frac{16}{225}.
\]

For any target n, set the marker energy to

\[
 \mu_{\rm sel}=-(n+2)E_n.
\]

Then (μ_sel+E_k)/(k+1) is uniquely minimized at k=n, with bulk energy
per material site −E_n. A wrong cell has bulk-subtracted excess

\[
 \delta_k=E_k+(k-n-1)E_n>0\quad(k\ne n),\qquad
 \min_{k\ne n}\delta_k=E_{n+1}.
\]

The minimum is over **all** nonnegative lengths, not a finite numerical
cutoff. The gap above the vacuum path ground state is at least 1−2q≥1/2.
Every payload bit string has the same path Hamiltonian, so the selected
cell code is exactly 2^n-dimensional at zero port and zero readout.

The proof, full transition table and the qualification for illegal sectors
are in [PROOF.md](PROOF.md). The exact selector addresses the length-selector
failure of a naive vacancy lift of C39's bare 46-state core.

## Arbitrarily many nonresonant leading bands

**[Derived here using C39's all-n pencil argument.]** Let

\[
 A(w)=4I-w\,\operatorname{adj}(P_n),\qquad
 f_j(t;w)=\int_0^t[e^{-sA(w)}]_{0j}
 [e^{-(t-s)A(w)}]_{j0}\,ds,
 \qquad \Gamma_{ij}(w)=\int_0^\infty f_i(t;w)f_j(t;w)\,dt.
\]

The positive Gram matrix Γ(w) has rank n. With different bus hoppings w
and wv and a relative hopping phase i, the leading logical hopping is

\[
 T=-W\,[\Gamma(w)+i\eta\Gamma(wv)],\qquad
 W=\lambda^2g^2|v_0v_n|^2D^{-2n-2}>0.
\]

This expression is the large-D, leading-g limit of the exact second-order
port coefficient, not the full finite-coupling Hamiltonian. Finite D gives
additional quadratic terms, and higher orders give interactions.

The readout needed for the all-n proof is supplied here by the one-site
projector onto L_1. If v_j are the normalized
vacuum amplitudes, its compression is

\[
 G^\dagger\sum_x |L_1\rangle\langle L_1|_xG
 =\sum_{i=0}^{n-1}\gamma_i n_i,
 \qquad \gamma_i=\sum_{j=i+1}^{n}v_j^2
 =q^{2i+2}[1+O(q^2)].
\]

Thus q=(wu)^2 supplies an independently tunable geometric profile. After
the same positive diagonal congruence as C39, the Bloch matrix is

\[
 -2C(w)\cos k+2\eta D_v C(wv)D_v\sin k+\nu E(q),
\]

where C(w)=D_w^{-1}Γ(w)D_w^{-1}, D_w=diag(w^{2i}),
D_v=diag(v^{2i}), and E(q)_ii=γ_i/(q²w^{4i})∼u^{4i}.

For each fixed n, the generalized eigenvalues of the two Gram matrices
are κ_jv^{4j}[1+o(1)] with κ_j>0. Their sea centers are consequently
separated on geometric scales. First choose u small to separate the
readout slopes, then v small to separate the centers, then ν small and
nonzero to separate all distinct unordered pair sums. This gives 2n simple
Fermi points with positive bilinear and quartic momentum margins for every
n. The existing proof applies because the new readout has the same leading
geometric asymptotics; its exact exponential form is unnecessary.

These margins survive sufficiently large finite D and sufficiently small
finite g and λ **at the finite-order matching level**. This continuity
statement does not by itself prove that the interacting infrared theory
has n modes.

## Checks performed

`checks.py` and `checks.json` retain the calculations and their precision.

- Literal 20- and 22-state bonds at n=1,2,3: Hermiticity, charge conservation,
  marker-number conservation, occupied-number conservation, and all local
  arrow-preservation conditions pass.
- An independent grammar enumeration reproduces the complete legal cell
  space (9n+5)2^n through n=4. Every allowed material pair admits exactly
  one coloring of each intervening vacancy string of lengths 1,2,3;
  forbidden pairs admit none.
- The selector inequalities were evaluated at 160 decimal digits for
  k=0,...,24, three q values and targets n=1,...,12. Every minimum occurs
  at the predicted length. The analytical proof covers the infinite tail.
- Complete cell spectra through n=4 reproduce the 2^n code multiplicity,
  the auxiliary gap, and the projected L_1 readout. Independently assembled
  legal vacuum rings of lengths 2,3,4 reproduce the predicted ground energy
  and translated-code multiplicity.
- Literal full-payload pair calculations at n=1,2 use excited resolvents
  of dimensions 288 and 3200. The λ²g² coefficient is Gaussian, has no
  nonscalar onsite term, and agrees with an independently constructed
  scalar-tree resolvent. Relative errors are below 10^-9, with the actual
  tighter values in the saved results.
- The large-D Gram limit was checked at n=1,...,4 and D=100,1000,10000.
  The check also fixes the factor of two: 2D^{2n+2}C_n→Γ, where C_n is
  the scalar exchange coefficient before the two charge paths are summed.
- The limiting quadratic pencils at n=1,...,6 have respectively
  2,4,6,8,10,12 distinct roots and nonzero separations of all unordered
  pair sums. Calculations use 145 through 370 decimal digits. These are
  high-precision numerical checks, **not interval certificates**, and
  are not measurements of the microscopic central charge.

At the deliberately conservative displayed schedule, the pair-sum margin
decreases from approximately 5.12×10^-62 at n=1 to 1.17×10^-205 at n=6.
This is an existence construction with extremely small energy scales.

## Reproducing

Dependencies: Python, NumPy, SciPy and mpmath, listed in
[requirements.txt](requirements.txt). Recorded versions are in the JSON
result. Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python3 -B contrib/smaller_d_unbounded_20261116/checks.py
```

The run writes `checks.json` beside the script. It took about five seconds
in the recorded environment. The 22-state matrix can be constructed with
`model.mobile_bond(model.parameters(n), ...)`. Those default couplings are
for resolved finite-matrix diagnostics, not a certified many-body example
with c=n+1. The parameter-order existence argument is in the proof.

## What would establish the example

1. Establish the selected material ground sector at nonzero port coupling
   with a bulk-subtracted defect comparison.
2. Prove or independently justify the fast-vacancy matching and preserve
   the positive material-defect cost in the mobile model at fixed J/t>0.
3. Transfer the leading band construction to the interacting multichannel
   conformal limit, including the entanglement or spectral identification
   of c=n+1.

The frozen d=20 family already has the arbitrary-n finite-cell and band
ingredients. The two-color lift is needed for the intended mobile,
one-site-translation-invariant ground-state interpretation on arbitrary
large physical rings. Neither frozen translation-breaking vacua nor color
choices are counted as additional bosons.

## Relation to earlier work

- [C39](../../codex/work_C39/PROOFS.md) supplies the two-Gram-matrix all-n
  nonresonance argument. Its 46-state initialized core lacked its own
  length selector; the audited mobile wrapper has d=4464.
- [C75](../../codex/work_C75/PROOF.md) supplies the two-charge replacement
  of the neutral/positive/negative triplet and the Gaussian CAR identity.
- [C80](../../codex/work_C80/PROOF.md) supplies displaced-bit storage with
  bitless tokens. Here the vacuum token and terminal marker are retained
  as separate states to preserve exact grammar invariance.
- [C88](../../codex/work_C88/PROOF.md) and the
  [vacancy-route audit](../route_audit_20261109/README.md) explain the
  mobile matching and dressed-defect obligations. The d=13 finite member
  there does not yet have an all-n nonresonance proof.
- [Latorre and Sierra](https://arxiv.org/abs/2403.17242) state the original
  c–d conjecture. No unconditional refutation is claimed here.
