# Independent saved-state diagnostics

OpenAI Codex, GPT-6 family, 2026-11-10. All results in this folder are
postprocessing of explicitly identified MPS checkpoints. No DMRG optimization,
old-data modification, commit or push was performed by this diagnostic task.

The strongest candidate examined here is the new R=8.75 bridge,
`(Jx,Jy,Jz,D,hz)=(1/2,-1/2,5/8,1/10,-7/8)`. Its even-branch entropy size
slopes through L192 stabilize near 1.45, agree with the L192 profile fits,
and persist after removing the changing binary parity contribution.
Transverse correlations and central density modulation both decrease with
size. The registered L192 refinement from chi192 to256 changes the last
fractional-window slope by only 0.00214, to 1.45303. This is promising finite-size evidence,
not a verified central charge or proof of a unitary CFT.

The old R=28 bridge remains less clean. Its entropy growth is not primarily
a changing binary parity contribution either, but persistent density
modulation and a nearby exact non-CFT weighted-Dicke control require care.
The preferred-ratio attractive reference instead weakens at L192 and shows
increasing evidence of an ordered crossover.

## Inputs and reproduction

`analyze_saved.py` reads completed records from
`contrib/d2_order_one_20261110/dmrg` and the checkpoint paths they contain.
It chooses the lowest recorded variational energy within each size/parity and
then the lower variational-energy parity for cross-size tables. This is not a
certification of the true ground-state parity. Input record/checkpoint SHA256
hashes and exact source paths are saved in `saved_diagnostics.json`.
The checkpoint coupling/history metadata are checked. Reconstructed entropy
profiles agree with the frozen JSON profiles to better than 1e-9.

From the repository root, with NumPy and SciPy in the existing task venv:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/cd-conjecture-venv/bin/python contrib/d2_order_one_continuation_20261110/diagnostics/analyze_saved.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /private/tmp/cd-conjecture-venv/bin/python contrib/d2_order_one_continuation_20261110/diagnostics/window_order.py
/private/tmp/cd-conjecture-venv/bin/python contrib/d2_order_one_continuation_20261110/diagnostics/check_mappings.py
```

Run each command sequentially, using one worker. The scripts set the thread
limits themselves before importing NumPy. `PROTOCOL.md` separates the initial
registered checks from the adaptive whole-window extension.

Names below: R28 = `(Jx,Jy,Jz,D,hz)=(1/2,-1/2,1,1/16,-7/4)`;
R14 = `(1/2,-1/2,3/4,1/12,-7/6)`; ordered A =
`(1/2,-1/2,3/8,1/4,-1/4)`; near-one B =
`(1/2,-1/2,1/2,3/16,-1/2)`. X/Y/Z always mean physical Pauli
operators, undoing the solver's site-dependent gauge.

## Entropy resolved by left-block parity

For Schmidt probabilities `p_(q,a)`, let `w_q=sum_a p_(q,a)`.
The exact identity is

`S = H(w_even,w_odd) + sum_q w_q H(p_(q,a)/w_q)`.

The first term is at most ln(2). The following are the original fractional
window slopes `6 Delta S / ln(L2/L1)`, decomposed term by term:

| State | Sizes | Total | Binary parity | Conditional |
|---|---:|---:|---:|---:|
| R28 | 24→48 | 3.02421 | 0.14402 | 2.88019 |
| R28 | 48→96 | 2.68340 | 0.04675 | 2.63665 |
| R14 | 24→48 | 6.86759 | -0.00485 | 6.87244 |
| Ordered A | 24→48 | 0.13471 | 0.06202 | 0.07269 |
| Near-one B | 24→48 | 1.00182 | 0.24307 | 0.75875 |
| Near-one B | 48→96 | 1.00371 | 0.24412 | 0.75959 |
| XX control | 24→48 | 0.96905 | 0.20089 | 0.76816 |
| Two-XX control | 24→48 | 1.71086 | 0.02568 | 1.68518 |

Only about 1.7% of R28's last total slope comes from the binary contribution.
This rejects the narrow explanation that its growth comes mainly from the
development of a two-sector parity entropy. It does not exclude translation
cats, domain-wall superpositions, crossover, or order coexisting with gapless
correlations. The conditional slope is **not** a central charge. Critical XX
also has a nonzero binary contribution, so equal parity weights alone are not
evidence of a cat. Both original entropy windows are retained in
`entropy_decomposition_slopes.json`.

At the R28 central cut, L48 odd parity has equal sector weights and almost
identical sector spectra. L96 even parity instead has weights 0.53346/0.46654
and visibly different leading spectra. This distinction is not itself a
phase diagnostic: the physical Hamiltonian has antiunitary symmetry RK,
where R reflects the chain and K complex-conjugates in the Z basis. R and K
each reverse the DM term. For an RK-symmetric odd-parity pure state at the
middle cut, reflection plus the Schmidt correspondence pairs the two parity
blocks. Thus the odd-parity pairing need not indicate topological edge states
or an ordered cat.

## Transverse envelope and whole-window order

The stored correlator starts at L/4. Its RMS over distances L/4 through L/2
avoids treating a single oscillatory endpoint as the envelope.

| State | L24 | L48 | L96 |
|---|---:|---:|---:|
| R28, X or Y far RMS | 0.08861 | 0.17623 | 0.11946 |
| Ordered A, X or Y far RMS | 0.38319 | 0.39629 | — |
| XX control, X or Y far RMS | 0.19700 | 0.13623 | — |

For R28 the last-size decrease is present within both fixed parities: the
effective two-size envelope exponents are 0.615 (even) and 0.544 (odd).
These are not asymptotic exponent estimates. R14's RMS falls even more
dramatically, concurrently with a strong profile change.

`window_order.py` computes `Q=sum_(j in central half) exp(i q j) O_j` and
`m²=<Q†Q>/W²`, where W=L/2, by exact transfer recurrences. Its six-site
complex-MPS check against a full state vector has error 2.7e-15. For the
adaptive channel q is selected from the previously stored signed Y correlator;
it is not optimized on a newly measured full correlation matrix.

| State | L | q for Y | m² at q | m² at uniform Y / staggered X |
|---|---:|---:|---:|---:|
| R28 | 24 | 0.32044 | 0.14960 | 0.19437 |
| R28 | 48 | 0.18692 | 0.14313 | 0.09432 |
| R28 | 96 | 0.19321 | 0.10833 | 0.03511 |
| Ordered A | 24 | 0 | 0.44321 | 0.44321 |
| Ordered A | 48 | 0 | 0.42662 | 0.42662 |
| XX control | 24 | 0 | 0.39915 | — |
| XX control | 48 | 0 | 0.29141 | — |

R28 has a slowly decaying incommensurate component. A plateau has not been
demonstrated. Its structure factor `<Q†Q>/W` at the selected q grows from
3.435 to 5.200 over 48→96; growth of this quantity is compatible with an
algebraic correlation function and is not sufficient for order.

The equality of X and Y envelopes on the symmetric bridge is also not two
independent modes. For Jx=-Jy, the antiunitary symmetry combining K with
staggered z rotations by ±pi/2 maps X_j to (-1)^j Y_j. A symmetric state
therefore has `<X_i X_j>=(-1)^(i+j)<Y_i Y_j>`.

## Density modulation and finite-size warnings

Central-half density Fourier peaks after removing a constant and linear trend:

| State | L | q_density | Fourier amplitude | Detrended RMS |
|---|---:|---:|---:|---:|
| R28 | 24 | 2.76303 | 0.29664 | 0.19060 |
| R28 | 48 | 2.78659 | 0.15932 | 0.11982 |
| R28 | 96 | 2.77717 | 0.15785 | 0.10865 |
| R14 | 24 | 2.56040 | 0.01690 | 0.01645 |
| R14 | 48 | 2.95938 | 0.11560 | 0.07439 |
| Near-one B | 96 | 2.77246 | 0.20622 | 0.14237 |

R28's density amplitude barely decreases over the last doubling. This could
be slowly decaying Friedel oscillations or pinned density order; the present
data do not distinguish them. The similar wavevector in near-one B warns
against identifying a mode count from this modulation. R14 changes its
modulation markedly between 24 and 48, making its enormous entropy slope
especially vulnerable to a shell/branch/crossover interpretation.

For R28, `(pi-q_density)/2` is 0.17750 at L48 and 0.18221 at L96, close to
the selected Y wavevectors 0.18692 and 0.19321 at the available window
resolution. These can be harmonics of one incommensurate scale; their
coexistence does not demonstrate several independent soft sectors.

The frozen R28 L96 chi192 state has profile estimates 2.65 and 2.09 in the
two original windows, both with RMS residual about 0.023. Its last-sweep
energy/entropy changes pass the prior sweep gate, but larger-chi convergence
was not established by that dataset. The parent's continuation addresses it.

## Exact alternate representation

Rotating odd sites by pi about X maps this Jx=-Jy family to

`sum_j [ Jx(XX+YY) - Jz ZZ - (-1)^j D(XY+YX) ] + hz sum_j (-1)^j Z_j`.

At D=0 the transformed model has U(1) symmetry, corresponding to conservation
of staggered Z in the original frame. At R28 the transformed XXZ anisotropy
is -2 and the staggered field is -7/4, with a small alternating pairing term.
This is an exact change of basis, not an enlargement of the physical unit
cell or a no-go theorem. `check_mappings.py` verifies the mapping and the two
antiunitary symmetries on four sites to 1.8e-15.

## Discriminating follow-ups

1. Establish chi and seed stability at the largest sizes and compare fixed
   parities as well as the lower-energy envelope. A parity splitting alone is
   not a bulk gap; a small energy variance alone does not certify its sign.
2. Add intermediate lengths near the last sizes, rather than only doublings.
   A smooth entropy trend and stable q are stronger than a favorable shell
   sequence. Retain the same fitting windows and report oscillatory residuals.
3. Measure central-window order at the incommensurate q and the connected
   density structure factor across sizes. Endpoint correlations and density
   homogeneity alone can both miss cat order.
4. Use small end-Z fields to change density pinning without breaking the
   solver's exact parity. Compare extrapolated central modulation and entropy.
   An end-X field would require a solver allowing parity mixing.
5. For coexistence, use small coupling offsets, opposite initial seeds and
   energy derivatives/order observables. Abrupt branch switches and hysteresis
   are qualitatively different evidence from a stable critical interval.
6. A claim of several modes needs an independent bulk diagnostic, such as
   low-energy finite-size spectra and velocities in multiple sectors. Entropy
   plus several Fourier peaks is insufficient.

An honest *promising candidate* has reproducible above-one entropy growth in
both windows after numerical convergence, persists through an additional
size step, and survives these elementary order/shell checks. A verified
c>1 unitary CFT requires substantially more than this screen.

## Preferred attractive route extension

`analyze_attractive.py` incrementally reads completed `pref_ref`, `pref_rat`,
`pref_zrat`, and `attr_ref` records in the neighboring attractive scan, plus
the archived smaller-size `attr_ref` records. It saves the same
entropy decomposition plus X/Y/Z order at q=0, pi and adaptive q to
`attractive_diagnostics.json`. It records hashes and discards measurements
if either input changes during the contraction, because the owning agent may
later resume a checkpoint at larger chi. A row is replaced when the source
record hash changes. This extension is independent postprocessing of the
other agent's selected points, not an independent candidate selection.

For context, Peotta et al. report first-order Luttinger-liquid–ferromagnet
boundaries in this model family and leave the meeting region unresolved.
That motivates order and energy-derivative controls without assigning the
new phi=1.2 point a phase from a phi=1 figure.
Source: [Peotta et al., sections 3.1–3.2](https://arxiv.org/html/1403.4568v3).

First completed preferred-reference comparison (chi128 at L48 and L96,
lower variational branch even at both sizes): the fractional slope 1.42149
decomposes into binary 0.43650 plus conditional 0.98499. The central-three-cut
version is 1.36483 = 0.29018 + 1.07465. About 31% of the fractional growth is
changing parity weights, substantially more than for R28, although this
alone cannot assign a physical central charge. The central-half staggered-X
order diagnostic falls 0.32640→0.28295, even while the original single-origin
far-envelope RMS grows 0.17276→0.18967. Thus the endpoint evidence does not
establish a plateau; larger sizes are needed to distinguish onset of order
from algebraic behavior. The odd branch gives m² 0.37794→0.30963.

The rational preferred point reproduces this pattern at L48/96: fractional
1.427594 = 0.446552 binary + 0.981042 conditional, with even-branch
staggered-X m² decreasing 0.337799→0.298782.

The preferred reference's completed even L192 state weakens the lead:
the 96→192 fractional slope falls to 1.142452 = 0.319231 binary + 0.823221
conditional (central-three-cut 1.091919). Whole-window staggered-X m² falls
only 0.282950→0.267526, a 5.45% drop versus 13.3% over the previous pair;
the single-origin endpoint grows. Binary entropy reaches 0.686223, leaving
at most about 0.060 in its slope contribution over one more doubling.
Combined with the original profile estimates dropping below one, this is
substantially more consistent with an ordered crossover than with a stable
above-one conformal regime. It remains a finite-size diagnosis, not a proof
of a gap.

The broader attractive reference now has a completed even L192 state. Its
96→192 fractional slope is 1.310805 = 0.349998 binary + 0.960807 conditional;
the central-three-cut version is 1.290066 = 0.233011 + 1.057055. The
staggered-X whole-window m² falls 0.244717→0.214373 over that pair, even
though the single-origin endpoint grows. The decay is slowing, so both
order onset and algebraic behavior remain plausible at these sizes.

The independent archive audit is in `ARCHIVED_SHORTLIST.md`, with its
recomputable underlying rows in `archived_lead_audit.json`.

## New intermediate bridge and independent Rényi check

For `(Jx,Jy,Jz,D,hz)=(1/2,-1/2,5/8,1/10,-7/8)`, R=8.75, the completed
24→48 even-branch fractional slope 1.703756 decomposes as 0.011200 binary plus
1.692556 conditional; the odd branch gives 1.711698 with binary 0.000194.
The evolving binary contribution is negligible. Uniform-Y/staggered-X
whole-window m² decreases 0.447424→0.362263, with similar odd-branch values.
The central density detrended RMS decreases 0.009104→0.006496 and its Fourier
amplitude 0.013700→0.009227.

The completed L96 chi192 states strengthen these diagnostics. On the even
branch the 48→96 fractional slope is 1.436219 = 0.001735 binary + 1.434485
conditional; the odd branch gives 1.455205 with binary contribution 0.000019.
The whole-window transverse m² falls further to 0.289044 (even) / 0.289093
(odd). The even-branch far XX/Y RMS is 0.203217, down from 0.276539 at L48.
The density detrended RMS falls to 0.001596 and Fourier amplitude to 0.002330.
Thus this lead survives the simple evolving-parity-weight and large-density-
modulation explanations better than R28 or the attractive candidates. The
entropy slope still trends down from the first size pair, and a plateau in
the large transverse correlation amplitude has not been excluded at later
sizes. These diagnostics do not certify c=3/2. Input hashes remained stable
during measurement; read-only results are in `bridge_diagnostics.json`.

The even L96 refinement from chi192 to256 changes fractional entropy by
only 8.11e-6, transverse m² by 6.35e-7, and density detrended RMS by 1.61e-8.
The earlier chi192 diagnostic row is preserved in `superseded_rows`, with
its matching frozen record path, and the latest row is chi256. Before that
refinement, L times the parity splitting is 0.0066055,0.0060714,0.0064209
at L24,48,96. This is compatible with inverse-size splitting over those
sizes; it is not a velocity estimate, and the absolute splitting is small.

At L192 chi192, the even-branch 96→192 slopes are 1.450888 (fractional) and
1.445203 (central-three-cut), using the refined L96 chi256 state. The binary
contribution is slightly negative, -0.004711, leaving conditional 1.455599.
The two L192 profile estimates are 1.45145/1.45986, with residual RMS about
0.0020/0.0018. Whole-window transverse m² falls further to 0.223056 and far
XX/Y RMS to 0.142527. The latter has a two-size effective exponent about 0.51;
this is not an asymptotic exponent determination. Density detrended RMS is
0.0006196 and Fourier amplitude 0.0009106. Thus the last two size pairs agree
far better than the first pair, while the order diagnostics continue to
decrease. All input hashes remained stable during the L192 measurement.

The completed odd L192 branch gives fractional slope 1.510321 over 96→192,
with binary contribution -2.1e-6. Its transverse window m² is 0.222972,
essentially the even-branch value 0.223056. At L192, L times the parity
splitting is 0.00714585, compared with 0.00661/0.00607/0.00642 at the previous
three lengths. The selected even state is explicitly linked to its
hash-verified preserved chi192 checkpoint in the diagnostic JSON.

The final L192 chi256 refinement passes the same input-hash and entropy-profile
checks. Relative to chi192, fractional entropy changes by +0.00024702,
energy by -3.58e-7, transverse window m² by +1.16e-5, and density detrended
RMS by +1.70e-7. The final 96→192 slopes are **1.453026 fractional and
1.447459 central-three-cut**. The fractional binary contribution is
-0.004700, leaving conditional 1.457726. Final L192 profile estimates are
1.453933/1.462960, with residual RMS 0.002012/0.001810. Transverse m² is
0.223068, far XX/Y RMS is 0.142569, and density detrended RMS is 0.00061977.
Thus the registered finite-chi check preserves the positive assessment;
larger sizes and independent conformal-spectrum evidence remain necessary
to identify the limiting theory.

The registered post-hoc energy check is in `energy_consistency.py/json`.
`D(L)=E(4L)-3E(2L)+2E(L)` cancels bulk and boundary constants. For this bridge,
D(24)=-0.00176176447 and D(48)=-0.00085910337 give
`p_eff=-log2[D(48)/D(24)]=1.036117`, versus 0.949698 for exact XX on the same
sizes. The form `E=eL+es+a/L` fitted at 24/48/96 predicts E192 within 2.18e-5.
This supports an inverse-size energy correction; it is neither a dynamical
exponent nor a velocity measurement. The broader attractive reference gives
p_eff=0.389865 with the same analysis, illustrating its lack of a comparable
clean energy trend on these lengths.

Preservation limitation: the first L96 even chi192 checkpoint was overwritten
before an NPZ copy was made. Its completed diagnostics and raw JSON remain,
but that exact old state cannot be reloaded. Other preserved bridge states
are listed in `../bridge_scan/state_preservation.json`, with the original
L192 even chi192 state in `../bridge_scan/largest_size_preservation.json`.
All eight final size/parity states are archived in
`../bridge_scan/final_state_preservation.json`; the diagnostic rows link
matching record/checkpoint snapshots by verified hash. The owning agent's
`STATE_PRESERVATION.md` records this limitation explicitly.
The parent also archived 12 completed candidate states, including the
pre-refinement L192 even chi192 state and independent seed/size controls,
in durable workspace storage `../checkpoints_local/`. Exact checkpoint
hashes and paths are recorded in `../local_checkpoint_manifest.json`.
These binary files are Git-ignored; the manifest can be matched to the
diagnostic input hashes if the temporary snapshot directory is removed.

`renyi_saved.py` and `renyi_diagnostics.json` are an independent replication
of a check also performed by the parent. The canonical continuation outputs
are `../renyi.py` and `../renyi_exact_controls.json`, including exact controls
through 192 sites. No duplicate continuation of the Rényi calculation is
planned here. On the new bridge 24→48, the normalized fractional-window
slopes for alpha 1,2,3 are 1.70376,1.43205,1.30143. Their mismatch is a
finite-size warning, not an exclusion: the short two-XX control also has
large Rényi-index dependence. The coefficients follow the open-boundary
replica formula in [Calabrese and Cardy, section 3.2](https://arxiv.org/html/0905.4013v2).

## Additional controls and limits

The zero-field reserve `pref_zrat`, `(Jx,Jy,D,Jz)=(7,-3,6,-9)/12`, gives a
useful caution. Its even-branch uniform-Z window m² is 0.406341 at L48 and
0.386568 at L96, while far ZZ RMS grows 0.276584→0.318660. Its entropy slope
1.265862 is still mostly conditional (1.151463, with binary 0.114399).
Z order preserves total parity, so parity decomposition cannot exclude this
ordered/cat explanation. This is why the separate transverse/density order
tests are required for the bridge lead.

`factorization_control.py` checks another exact reference algebraically on
L4/6/8 rings. For Jx=-Jy=1/2, Jz>1/2 and
`hz=-2 sqrt(Jz²-1/4)`, a uniform product state with
`<Z>=sqrt((Jz-1/2)/(Jz+1/2))` and transverse moment along Y is an exact
ring eigenstate, even with nonzero uniform D; its residual is below 8e-15.
At D=0, the shifted bond `h_bond+Jz I` is positive and rank one, with
eigenvalues 0,0,0,4Jz. Restoring the actual D values changes the ground state:
on L8 **at the factorization-line fields** (-sqrt(3) and -3/4, respectively),
the product state's energy exceeds the exact ground energy by 0.03154
for the R28 exchange/D values and 0.13300 for the R8.75 exchange/D values.
These fields differ from the actual candidate fields -7/4 and -7/8. This
check neither excludes nor establishes order at the actual candidate fields.
Raw results are in `factorization_control.json`.

The registered finite-ring degenerate-perturbation extension is
`factorization_projection.py/json`. At D=0 on the exact line, the ground
space has dimension L+1 for L4/6/8. The compression of the DM perturbation to
that space vanishes to 8e-15. Its second-order lowering at L8 is -0.03139
(R28 exchange/D values) and -0.13013 (R8.75 values), close to the exact
lowering -0.03154 and -0.13300. However, the gap above the degenerate space
decreases with L (0.536→0.247 and 1.000→0.426 over L4→8). This prevents a
uniform perturbative inference from these short rings. The calculation
suggests a low-energy competition near the factorization line, but supplies
neither a reliable thermodynamic field shift nor a count of soft modes.

`dicke_control.py/json` verifies a sharper exact reference on that D=0 line.
After odd-site X rotations, the ground state with transformed occupation
N=L/2 is a weighted Dicke state with site amplitudes
`w_even=-(2Jz-sqrt(4Jz²-1)), w_odd=1`. Polynomial coefficients give its exact
Schmidt probabilities, checked against full vectors on L4/6/8. The
state is a ground state of the periodic D=0 parent and of its open sum of
positive bond projectors, whose end fields are hz/2 rather than the full hz
used in the candidate search. This boundary distinction is intentional and
must be retained when comparing the control. Its
fractional-window OBC-normalized slopes over 24→48→96→192 are
3.02637,3.00578,3.00137 at Jz=1 and 2.93851,2.96766,2.98346 at Jz=5/8.
The binary parity contribution is negligible after L48. These degenerate
ground states are an explicit non-CFT entropy-growth control in the same
family: a slope above one plus negligible binary growth cannot by itself
establish a conformal theory. This does not identify the nonzero-D states
with the Dicke reference. Their Rényi trends are also qualitatively different
on the new bridge: a Dicke leading term S_alpha≈(1/2)lnL gives normalized
slopes 3,4,4.5 for alpha 1,2,3, whereas the bridge's 96→192 slopes at the
L192 chi192 stage were about 1.451,1.479,1.369 in the canonical parent report.
The smaller size pairs
decrease with alpha, but this monotonicity does not persist at 96→192.
That contrast is supporting diagnostic information, not a
universal exclusion of ordered or nonconformal crossover.

The degeneracy has a direct algebraic explanation: after the odd-site
rotation, each positive bond projector relates configurations differing by
an adjacent 01↔10 exchange. All configurations with fixed N are connected,
so their amplitudes are uniquely fixed up to normalization by the displayed
weights. There is one orthogonal exact ground state for every N=0,...,L.
The resulting L+1 ground-state degeneracy is incompatible with interpreting
these particular entropy slopes as those of an ordinary finite-c unitary
CFT with a fixed vacuum degeneracy.
