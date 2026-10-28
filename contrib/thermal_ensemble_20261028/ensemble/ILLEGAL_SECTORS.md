# Positive completion of the full rooted cell space

Status: **[proved local repair and full-space reduction]**, conditional only
on the stated finite C44 cell-gap and wrapper bounds. This modifies the
candidate bond on illegal selector words. Its alphabet remains 4464 and
its complete selector-legal Hamiltonian is unchanged. It is not a proof
of the final infrared/CFT step.

Provenance: OpenAI GPT-6-family Codex ensemble-transfer agent, 2026-10-28.
The parent agent and this agent independently identified the missing
conditional pin and its positive completion. A separate agent checked
the local rules, full-cell bounds, readout invariance and SW algebra in
the [independent audit](../exact/COMPLETION_AUDIT.md). No human review is
claimed.

## 1. The literal obstruction and the repaired bond

C8 `helix_equalizer.md`, Section 1, and C44 `PROOFS.md`, Section 1, place
the endpoint pins on the selector projector

    P_cut=|X#><X#|+|R#><R#|.

With `a_B` the positive pin on track B of the left physical site and
`a_A` that on track A of the right physical site, the pin bond is

    P_cut tensor (a_B tensor I+I tensor a_A).           (1)

The second term conditions the first-site pin of a rooted cell on the
last selector letter of its predecessor. Outside the selector-legal
space it is not a sum of independent one-cell operators. Thus the
full tensor-product onsite representation does not follow merely by
including all strings in the cell Hilbert space.

**[proved repair]** Replace (1), and only (1), by

    P_into_marker tensor (a_B tensor I+I tensor a_A),
    P_into_marker=I_selector tensor |#><#|.             (2)

Equivalently, add the repeated positive two-site term

    delta h=(|L#><L#|+|##><##|)
                 tensor (a_B tensor I+I tensor a_A).   (3)

All omitted physical factors carry identities. The original positive
pin coefficients, including any `J_h`, are part of `a_A,a_B`.

Equation (3) is Hermitian, positive, nearest neighbor, one-site
translation invariant, and conserves marker number and payload charge.
It changes no physical dimension. It vanishes identically throughout
every selector-legal sector, not just after projection onto the code.
All original unperturbed sector lower bounds remain valid because the
added operator is positive. Ports, readout, marker motion and all
coefficients of the intended selector-legal dynamics are unchanged.

## 2. Exact onsite decomposition, including consecutive markers

Use the complete cell space from `PINNED_MARKER.md`: a cell starts with
`#` and contains no further `#`. Its marker site retains all 1116
internal states. Assign to the cell every internal physical bond, the
full marker chemical potential `mu_sel`, and the selector penalty on
the bond from its final letter to the next `#`. This last penalty
depends only on the cell's final selector letter, since the incoming
selector arrow of `#` is fixed.

The wrapper and helix hopping terms vanish on a bond entering `#`, so
they lie inside individual cells. Selector head hopping `LX<->XR`
also lies inside a cell. In (2), assign the B pin to the cell on its
left and the A pin to the cell on its right. Both are now unconditional
endpoint pins. The unperturbed repaired Hamiltonian is consequently
exactly a sum of one-cell operators `A_raw,ell`, before the cyclic
quotient is imposed.

For `ell=1`, the cell is the single marker site. It receives both pins,
the marker chemical potential, and the positive `##` selector penalty.
No special unpinned or zero-cost length-one sector is omitted.

Subtract C44's coexistence energy `e<0` per physical site:

    A_ell=A_raw,ell-e ell I.                            (4)

The subtracted constant becomes `e L` at fixed physical length and can
be restored after the pressure calculation.

## 3. A linear lower bound on every illegal selector cell

Let `z` count selector violations on the cell's internal bonds and its
outgoing bond into the next marker. Selector arrows are

    # : 1->0,  L : 0->0,  X : 0->1,  R : 1->1.

All violation projectors commute with the unperturbed terms. Cut the
cell at its internal violations. There are at most `z+1` resulting
legal open strings. Every string has at most one selector head, whose
hopping energy is bounded below by `-2 t_sel`. Its open wrapper
Hamiltonian is bounded below by `-4U`, including arbitrary malformed
core words, by C44 Section 3 / C8 `arrow_wrapper.md` Section 3. The
extra endpoint helices and completion are positive. The chemical
potential `mu_sel` occurs once and is positive.

Writing `C=2 t_sel+4U` and using C44's
`K_sel=4 t_sel+9U`, for `z>=1` we therefore have

    A_raw,ell >= mu_sel + K_sel z - C(z+1)
              = mu_sel + U +(K_sel-C)(z-1)
              >= U.                                    (5)

This also covers violations on the final bond and consecutive markers;
counting a possible extra fragment only weakens the bound. From (4),

    A_ell >= U+|e| ell        on every z>=1 cell.         (6)

This is a full operator bound in each invariant violation sector, not
a bound on classical diagonal configurations. It provides the required
linear cost for long illegal cells.

For `z=0`, the repaired cell is exactly the original C44 open cell.
Writing `gamma=t_sel eta_b-U>0`, C44 gives

    A_ell >= gamma(1+|ell-b|), ell not in {b,b+1},
    ker A_b and ker A_(b+1): the complete 2^n payload codes,
    A_ell >= Delta Q_ell, ell in {b,b+1},               (7)

where `Delta>0` is C44's complete selected-length gap, including its
charge tunneling splitting. Equations (6)-(7) establish a complete
finite code `P=P_b+P_(b+1)` in the full direct-sum cell Hilbert space.

In particular, with

    Delta_full=min(Delta,gamma,U)>0,
    kappa=min(Delta,gamma,U,|e|)/(b+2)>0,

we have simultaneously

    A >= Delta_full Q,
    A_ell >= kappa(ell+1) Q_ell.                        (8)

For the wrong legal lengths the elementary inequality
`1+|ell-b| >= (ell+1)/(b+2)` suffices. At the selected lengths the
definition of `kappa` pays their finite size.

The number of cell states at length `ell` is at most `d^ell`, `d=4464`.
Consequently, if `beta kappa>log d`,

    Tr_Q exp(-beta A)
      <= sum_(ell>=1) d^ell exp[-beta kappa(ell+1)]
       = d exp(-2 beta kappa)/(1-d exp(-beta kappa)).    (9)

Thus the complete all-string auxiliary heat trace is finite in a
low-temperature interval. This includes the illegal selector sectors
that the uncompleted tensor decomposition did not justify.

## 4. The unbounded readout is a positive invariant onsite term

Define `N_tok,j` as the number of occupied core letters in cell `j`
carrying any token, vacuum or excited. The literal 46-letter alphabet
counts tokens on `d`, `x`, and occupied `A_s,B_s`, and counts zero on
`L,R,A_empty,B_empty`. Wrapper holes carry no core token.

The actual rules preserve this number for every ordered cell:

- Every core hopping rewrite in C39 `literal.py` moves one token inside
  the occupied word. The junction `d<->x` preserves it, and the onsite
  charge conversion changes only the species/data label.
- A source `F_(r,s)` changes `v` into an excited token or conversely at
  a fixed occupied endpoint. It never creates or deletes a token.
- Wrapper vacancy hopping moves the complete occupied label inside
  its cell, so its token count does not change.
- Marker hopping transfers only a wrapper hole `H_0<->H_2` between
  neighboring cells. It transfers no occupied label or token.
- Selector head hopping, filter and helix operations do not act on
  core token occupation.

These statements apply to malformed words as well. Hence

    R_j=1_{N_tok,j>=2}                                 (10)

is an invariant one-cell projector for every unperturbed term, port,
marker move and readout. It annihilates the complete code. This is a
smaller and more directly checked invariant than a vague classification
of all illegal core-arrow sectors.

Let `F_read,j` be the unit-coefficient sum of the physical readout
projectors. It counts a data bit only at a vacuum-token site. Therefore

    0 <= F_read,j <= N_tok,j,
    0 <= (1-R_j) F_read,j <= I.                         (11)

C45 (33) gives

    f_read=J0 mu S_q/(s q^2)>0,

since C42 chooses `mu>0` and every other factor is positive. Split the
physical readout exactly into the bounded one-cell term
`f_read(1-R_j)F_read,j` and the positive, possibly unbounded term

    C_j=f_read R_j F_read,j R_j >=0.                    (12)

Place `lambda^2 C_j` in the onsite operator, not in the bounded
perturbation norm. It leaves `P` unchanged and only improves (8)-(9).
Noncommutation of `C_j` with `A_j` causes no difficulty: finite-length
min-max inequalities, followed by the direct sum, give the heat-trace
comparison from `A_j+lambda^2 C_j >= A_j`.

All remaining intercell perturbations are uniformly bounded in the
full direct-sum space: the ports are fixed endpoint matrices with
fixed compensator norms; marker moves are bounded partial isometries;
the weak chemical adjustment is one constant per cell; and (11)
bounds the retained readout. No factor proportional to cell length
appears in those bounded norms.

## 5. Why this absorption preserves the local SW calculation

Here is the relevant general algebra. Suppose `R_j<=Q_j` commutes
with the onsite `A_j` and every interaction meeting cell `j`. For a
support `X` containing `j`, let `P_X` be its all-code projector. Then

    R_j P_X=0,
    R_j V_X P_X=V_X R_j P_X=0.                          (13)

The reduced resolvent of `A_X` also commutes with `R_j`. Therefore
C45's first generator

    S1_X=(Q_X A_X^-1 Q_X) V_X P_X
            -P_X V_X (Q_X A_X^-1 Q_X)

annihilates `R_j` from both sides. The same property is retained by
the off-diagonal part of the second-order interaction and its generator
`S2_X`. Consequently every generator commutes with `C_j`: for `j` in
its support both products vanish, and for `j` outside it they commute
by disjoint support.

Let `S=lambda S1+lambda^2 S2`. If the bounded calculation gives the
exact identity

    exp(S) (A+lambda V+lambda^2 F_good) exp(-S)=D+Remainder,

then it gives, without any new remainder,

    exp(S) (A+lambda^2 sum_j C_j+lambda V+lambda^2 F_good)
                  exp(-S)
      = D+lambda^2 sum_j C_j+Remainder.                 (14)

The reduced resolvents constructed with the enlarged onsite operator
agree on the code-connected matrix elements, because `C_j` vanishes
there. Thus the original projected second-order Hamiltonian and its
local remainder estimates are unchanged by the readout absorption.

The positive pin completion is a separate change of the physical bond.
It vanishes in the entire selector-legal invariant sector. All paths
starting in the code remain in that sector, since the actual selector
head and marker moves preserve selector violations and ports leave
selector letters unchanged. Therefore the projected coefficients and
the candidate's reference bands also remain unchanged by that repair.

## 6. What has and has not been closed

For the explicitly repaired bond (2), the full all-string rooted
representation has a positive product onsite reference with finite
code, uniform gap, exponential heat-trace control, and bounded residual
interactions. No illegal selector or multi-token sector needs to be
silently omitted to satisfy those structural hypotheses. The pressure
transfer in `PINNED_MARKER.md` applies to the repaired physical bond.

The claimed closure is this full-space microscopic reduction. The
operator-valued massive integration theorem, its precise coefficient
norms, the multichannel RG/thermal comparison, and the final CFT
identification must still be completed and checked. The positive
completion alone does not solve the c-d conjecture.
