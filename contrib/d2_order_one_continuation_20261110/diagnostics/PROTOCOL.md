# Saved-state diagnostics, registered before postprocessing

OpenAI Codex, GPT-6 family, 2026-11-10. No new variational optimization.
Read the frozen first-search records/checkpoints only; write all new files here.

Compare the R=28 and R=14 bridges, the known ordered flag
`bridge_r1.5_u1.5_k3`, the simple R=8/3 seed, and XX / two-XX controls.
Use the lowest variational-energy completed state in each recorded parity;
also select the lower-energy parity for cross-size summaries. This selection
does not certify an exact ground state.

1. Resolve Schmidt probabilities at every cut into the two left-block parity
   sectors. Decompose S = H(p_even,p_odd) + sum_q p_q S_q. Compare the binary
   contribution and conditional contribution in both original entropy windows.
   A developing two-branch parity cat can contribute at most ln2 to a cut;
   equal parity weights alone are not a cat diagnostic (critical XX also has
   them). Retain the leading central-cut eigenvalues by parity and their pairing.
2. Use the stored physical XX/YY/ZZ correlators, originating at L/4, to report
   RMS and mean absolute transverse correlations over separations L/8..L/4
   and L/4..L/2. These are envelopes, not fixed-point exponent determinations.
3. Report central density mean/range and windowed Fourier maxima, plus
   oscillation amplitudes after a constant/linear trend is removed. Do not
   interpret open-end oscillations alone as phase separation or long-range order.
4. Recompute the original profile fit windows, parity splittings and entropy
   slopes. Examine whether the within-parity entropy growth explains the high
   total slope. Any additional diagnostic selected afterward will be labeled.

The known ordered flag and exact-control labels calibrate interpretation.
No gap, central charge or number of modes is inferred from a single scalar.

## Adaptive follow-up, registered before its computation

The frozen R=28 transverse correlator changes sign on a long wavelength, so
testing only uniform or staggered order would miss an incommensurate component.
For the lower-energy stored parity at each size, evaluate central-half-window
`<|sum_j exp(i q j) O_j|^2>/W^2` for physical X at q=pi and Y at q=0, and
for physical Y at a q selected by the maximum of the already stored signed
Y-correlator Fourier magnitude (2001 points in [0,pi]). The latter is adaptive
to the old data and is a diagnostic, not a fitted asymptotic order parameter.
Do this with exact MPS transfer contractions, one single-thread worker and no
variational updates. Check the contraction against explicit tiny-chain state
vectors and the physical X/Y identity at q=pi/0 on the symmetric bridge.
Report the q and both normalizations `<Q*Q>/W` and `<Q*Q>/W^2`. A small finite
size trend does not establish long-range order or a number of critical modes.

## Preferred attractive candidates, subsequent read-only extension

Apply the same cut-parity decomposition and long-correlation RMS envelopes to
completed `pref_ref` and `pref_rat` records from the neighboring
`attractive_scan/results` folder, including both parities when finished. For
each physical X, Y, Z operator measure central-half-window Q†Q at q=0, q=pi,
and the maximum selected from its stored signed long connected correlator.
Keep both the raw value and the connected value subtracting |<Q>|², which
matters for density. These are fixed-state measurements, not fresh DMRG.
The selection of these candidates was made by the attractive-route agent;
this is independent numerical postprocessing, not independent candidate selection.
Check the checkpoint/record history consistency and hash the inputs; skip any
record whose checkpoint changes while it is being measured. Incremental reruns
may add finished sizes or replace a row if its input hash changed.

The parent subsequently requested the identical checks for the existing
broader Peotta reference: old `reference_peotta_exploratory` L48/96 records
and new `attr_ref` L192. They are grouped as `attr_ref` only after checking
the same physical couplings; original names/source paths are retained.

## New intermediate bridge lead, same diagnostics

The bridge-scan agent selected `(Jx,Jy,Jz,D,hz)=(1/2,-1/2,5/8,1/10,-7/8)`
after its24/48 screen. Apply exactly the same entropy decomposition, density
Fourier diagnostics and X/Y/Z whole-window structure factors to completed
records `bridge_scan/dmrg/b_r5_2_d1_10_q0_hm7_8_L*.json`, including96when
available. `analyze_attractive.py --family bridge` writes a separate
`bridge_diagnostics.json`. No fresh optimization is done by this task.

## Adaptive Rényi cross-check

To test whether the entropy excess depends strongly on the small Schmidt
probabilities, measure Rényi entropies alpha=2,3 on the same completed states.
Use the existing two windows; compare slopes with factors8and9 respectively
(the open-chain CFT expectation is `(c/12)(1+1/alpha) log L`). Agreement is
supporting evidence only: finite-size oscillations, boundary effects and
ordered crossovers can affect each Rényi index differently. No threshold is
tuned after output. Preserve snapshots keyed by input hash, and skip changing
inputs. Store this separately in `renyi_diagnostics.json`.

The parent subsequently added `pref_zrat`, the zero-field reserve with
`(Jx,Jy,D,Jz)=(7,-3,6,-9)/12`, to the preferred attractive queue. Apply
the identical completed-state diagnostics without changing any windows.

## Exact small-ring product-state control

An algebraic observation suggests a useful nearby control for Jx=-Jy=1/2:
at `hz=-2 sqrt(Jz²-1/4)` and Jz>1/2, a uniform product state with
`<Z>=sqrt((Jz-1/2)/(Jz+1/2))` and transverse moment along Y is a ring
eigenstate; the DM action telescopes on a uniform product. Verify this
directly for L=4,6,8, and compare with the exact lowest ring energy for the
R28 and new R8.75 couplings. This checks an eigenstate/control only. Do not
assume it is the ground state, and do not infer an ordered neighborhood or a
no-go statement from factorization. No new DMRG is authorized for this check.

At the parent's request, test the finite-ring perturbative structure at these
exact D=0 factorization-line fields: ground-space dimension, gap above it,
and the compression of unit DM and total Z to that subspace, at L4/6/8.
If the DM compression vanishes, report the second-order compressed operator
using the full excited-space resolvent. This is finite-ring degenerate
perturbation theory only; a size-dependent small denominator prevents any
uniform thermodynamic inference or numerical field-refinement recommendation.

## Exact weighted-Dicke reference on the D=0 line

The rank-one parent suggests a sharper non-CFT control. After odd-site X
rotation, fix total occupation N=L/2. Let `r=2Jz-sqrt(4Jz²-1)` and choose
site amplitudes `w_even=-r,w_odd=1`; the state sums all N-particle basis
strings with amplitude `prod_(occupied j) w_j`. Verify directly on L4/6/8
that the corresponding original-frame state is a ground eigenstate at D=0
and the exact factorization-line field. Its Schmidt probabilities follow
from coefficients of `prod_j(1+|w_j|² x)` in each side of the cut. Compute
the two existing von Neumann entropy windows for L24/48/96/192 and their
OBC-normalized slopes. This is an exact degenerate ground-state reference,
not a central-charge or mode-count claim for nonzero D. Compare its parity
decomposition to demonstrate whether a non-CFT entropy excess can survive
the binary-weight check. No new Rényi continuation or DMRG is performed.

## Post-hoc ground-energy consistency check

At the parent's request, once both bridge L192 parities are complete, compute
`D(L)=E(4L)-3E(2L)+2E(L)` for L24 and48, selecting the lowest recorded
variational parity at each size. This cancels terms linear and constant inL.
For an energy expansion `E(L)=eL+es+a/L^p+...`, their ratio yields
`p_eff=-log2(D(48)/D(24))`, when the ratio is positive. A conformal open-chain
Casimir term predicts p=1 under the applicable boundary assumptions. This
is a post-hoc consistency indicator, not a dynamical exponent or velocity
measurement. Compare with the exact zero-field XX chain on the same lengths,
using single-particle energies `-2 cos(pi m/(L+1))` at half filling.
If a denominator is tiny or the signs differ, report that directly instead
of forcing a power fit. Hash the input records and preserve numerical chi
differences as the relevant energy-accuracy check.
