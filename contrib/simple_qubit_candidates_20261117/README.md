# Simpler qubit candidates: calibrated full-chain diagnostics

Date: 2026-11-17. OpenAI Codex, GPT-6 family. **Bounded campaign complete:
46 runs, both candidate parity sectors, exact matched controls, and no
new established thermodynamic c>1 point.** The simplest surviving lead
is the older V=4,kappa=1,a=1 point: its previously missing L=96 chi
comparison now passes the local gates. This makes it the preferred next
test under the user's simplicity priority, not a verified CFT.

Raw records are under `outputs/`; `analysis.json` contains every gate
and parity/size comparison, and `followup_analysis.json` contains the
legacy check, parameter probes and independent start. A run marked
`done` means its bounded schedule finished, not that every gate passed.

## Claim and model class

We investigate Hermitian, one-site translation-invariant,
nearest-neighbour chains of physical qubits in the four-coupling family

```
H = sum_j [(XX-YY)/2 + kappa/(4V)(XY-YX) + (V/4)ZZ]_(j,j+1)
    + (-V/2+a/V) sum_j Z_j.
```

Coefficients are fixed independently of L and the full onsite field acts
on every open-chain endpoint. The target is a unitary compact conformal
bulk with total Virasoro c>1. No c_eff reinterpretation is used.
[MODEL.md](MODEL.md) gives exact source conventions and infrared limits.

The older C190/C195 production point remains the evidence benchmark,
with V=40,kappa=1,hz=-19.976047250613814. Its previously audited
256-to-512 main-window reading is 1.982737, but its last entropy
increment 0.010375 remains above the registered 0.01 threshold.
The [earlier registration](../qubit_convergence_20261116/REGISTRATION_20261116.md)
and conclusions are unchanged. The epsilon=10^-3 sixth-order construction
is retained as a [mechanism study](../qubit_resonant_liquids_20261116/README.md),
not the practical candidate for this simplicity-first search.

## Simplicity comparison

There are **four independent real coupling parameters, three after
overall scale fixing**. In the original common Pauli basis there are
six nonzero monomials: XX,YY,XY,YX,ZZ,Z. The same common onsite rotation
U=exp(-i pi Z/8) at all sites gives **four monomials**:

```
(1/2+D)XY+(1/2-D)YX+(V/4)ZZ+hz Z.
```

The total field is counted once, and constants and telescoping terms
are removed. No unequal coefficients are grouped to reduce the count.
Both literal coefficient spreads are reported with Delta=1:

| Follow-up rank / role | Point | Original spread | Rotated spread | Leading t2 | Evidence and limitation |
|---|---|---:|---:|---:|---|
| 1: simplest retained lead | V=4,kappa=1,a=1 | 28 | 4 | 1/8 | Both L=96,chi=256 sectors pass local gates; strong modulation and anomalous size fits remain |
| 2: dilute benchmark comparison | V=8,kappa=1,a=24/25 | 124.16 | 8.2773 | 1/16 | L=32,64,128 calibrated campaign; largest-size chi gate still fails |
| 3: no present simplicity advantage | V=12,kappa=1,a=24/25 | 284.16 | 12.3548 | 1/24 | Similar size evidence to V=8, with a larger hierarchy and failed chi gate |
| Existing evidence benchmark | V=40, production field | 3196.17 | 40.4578 | 1/80 | Strongest existing large-size evidence, still provisional |

This is a ranking of the next research targets, not established CFTs.
V=4 moves ahead only after its new convergence check; it was discovered
through the sibling source audit after the initial V=8/V=12 selection.
Its larger entropy slope is not a reason for the promotion. None of
these points has a verified critical parameter neighbourhood, and no
mode velocity is measured here; t2 is a leading generated hopping scale.

The rotation shortens the expression but encodes the small hopping in
the difference of the two cross couplings. Invariant ratios
(Delta/2)/D=2V and Jz/D=V^2 remain. A generated small kinetic scale or
mode velocity is distinct from the microscopic coefficient spread.
The [sibling term analysis](../qubit_minimal_terms_20261117/README.md)
independently confirms the count and leaves some simpler faces open;
it does not establish that all four structures are universally necessary.

The preferred simplicity lead is the explicit rational Hamiltonian

```
H = sum_j [(9/16) X_j Y_(j+1)+(7/16) Y_j X_(j+1)+Z_j Z_(j+1)]
    -(7/4) sum_j Z_j.
```

Equivalently its two-site bond is
`h=(9/16)XY+(7/16)YX+ZZ-(7/8)(ZI+IZ)`.
For open chains add `-(7/8)(Z_0+Z_(L-1))` to the bond sum. In the
original common basis this is `(XX-YY)/2+(XY-YX)/16+ZZ`, with the
same full onsite field -7/4. Its invariant pairing/DM and ZZ/DM ratios
are still 8 and 16. Its criticality and parameter tolerance remain open.

The explicit rational V=8 comparison Hamiltonian, in the common rotated basis,
is

```
H = sum_j [(17/32) X_j Y_(j+1)+(15/32) Y_j X_(j+1)+2 Z_j Z_(j+1)]
    -(97/25) sum_j Z_j.
```

It is exactly the original family with pairing scale one, D=1/32,
V=8 and a=24/25. Its finite-size evidence must not be read as proof
that this exact rational point has a c>1 thermodynamic limit.
Equivalently, a two-site bond is
`h=(17/32)XY+(15/32)YX+2ZZ-(97/50)(ZI+IZ)`.
For open chains, summing this bond requires the additional endpoint
terms `-(97/50)(Z_0+Z_(L-1))`; all runs use the full uniform onsite
field in the displayed Hamiltonian. The field is counted once in the
simplicity comparison.

## Numerical status

The new V=4 check is at the **same L=96** as the older calculation,
with an independent seed and chi=128,192,256 stages. At the common
chi=192 its energies agree with the old runs to 8.6e-9 or better.
The final chi-step entropy changes are 0.005120 (even) and 0.005618
(odd), below 0.01; both final stages also satisfy the strict stopping
rule. Variance per site is 1.52e-7 and 1.68e-7. The even sector
has the lower variational energy, by 0.00019166, but the global variance
does not rigorously certify such a small parity splitting.

This closes the missing local chi comparison, **not the CFT question**.
The even state's central density is 0.29636 and its central density
range remains 0.30147. The current smoothed profile estimator is
2.882/3.101 in even/odd sectors; the older 48-to-96 estimates span
2.51--3.00 with their parity/window choices. These are different
estimators and neither supplies a reliable central charge. Oscillatory
pair correlations also remain: the even connected separation-two pair
values at r=8,16,32 are -0.00521,0.01421,0.00734. A new calibrated
size sequence and parameter neighbourhood are still required for V=4;
this campaign's new matched-size sequence is at V=8 and V=12.

The initial a=1 pilots at V=8,12 were substantially denser than the
older candidate. Sequential rational-field bracketing selected the
common a=24/25,kappa=1. At L=32,64 both candidate ground branches pass
the residual, sweep-stationarity, discarded-weight and chi-step gates.
The main-window 32-to-64 slopes are about 2.4; the exact two-leg control
itself reads 2.389 there. Narrow and wide windows show large finite-size
effects, so those slopes alone do not identify a CFT.

The completed V=8 extension reaches chi=256 at L=128. Both parities
pass residual, sweep-stationarity and discarded-weight gates and satisfy
the stricter stopping rule, but **both still fail the chi-step gate**:
Sbar6 changes by 0.013937 (even) and 0.013014 (odd) from chi=192,
against the 0.01 threshold. V=12 remains at chi=192, with changes
0.051429 and 0.047392 from chi=128. This bounded campaign does not
promote either largest-size estimate to a converged CFT measurement.

The estimator is `c_two=6[Sbar_W(2L)-Sbar_W(L)]/log(2)`, where Sbar
averages the 2W+1 central entropy cuts. All entries below use 64-to-128;
the candidate rows use the lower reported-energy sector at each size.

| Data | W=1 | W=6 | W=12 | Accuracy status |
|---|---:|---:|---:|---|
| V=8, chi=192/256 | 2.39955 | 2.17284 | 2.23981 | Last chi-step fails |
| V=12, chi=192/192 | 2.25971 | 2.10901 | 2.18333 | Last chi-step fails |
| Exact one-leg XX | 1.02415 | 1.00403 | 1.04507 | Exact finite-size calibration |
| Exact two-leg XX | 2.15587 | 1.85827 | 2.08369 | Exact finite-size calibration |

For V=8 all four parity pairings span 2.17284--2.48030 in the main
window. The lower-energy branch switches from even to odd between
64 and 128; its central density changes 0.19944 to 0.19677, while the
central density range remains 0.16261 at 128. The density-corrected
secondary estimate 2.21753 is an interpolation here, but does not repair
the failed chi gate. At 32-to-64 it is an extrapolation and stays flagged.
The older V=8 infinite-chain reading near 1.4 remains in tension with a
clean two-mode interpretation. Re-fitting saved chi=64--192 states gives
1.40345 and 1.40092 for two- and four-site cells; their optimizer errors
remain about 1e-5--7e-5. Including less-converged higher-chi states lowers
the estimates further. Its fitted coefficients differ from this rational
point. Neither c=2 nor a smaller infrared CFT is settled by those data
(`source_audit.json`, `followup_analysis.json`).

The one-leg DMRG control at 128 agrees with its exact entropy profile to
7.2e-9. The two-leg control remains flagged: at 128,chi=192 its energy
drift is 1.22e-6 and Sbar6 error is -0.003144, giving a size reading
1.83144 instead of the exact 1.85827. Its exact value supplies the
calibration reference; the error is not transferred as an error bar to
an interacting candidate. The two-leg control has range two and is
used only for calibration, not as an admissible nearest-neighbour model.

`analysis.json` records densities, density modulation, valley weights,
connected separation-two pair correlations, parity gaps, every parity
pairing of size estimators, and the extrapolation flag on any secondary
density-corrected estimate. Balanced valleys do not rule out pairing.
No Luttinger exponent is fitted from a few short distances. Symmetry-
allowed two-pair charge locking remains a physical concern on this line.
The raw valley convention assigns both boundary momentum bins to one
half; `valley_half_boundary` splits those bins equally and removes the
spurious roughly 0.16% imbalance. Neither convention diagnoses charge
locking. The stricter registered stopping rule and the four local
diagnostic gates are reported separately; `schedule_limit_hit` remains
visible even when the looser local gates pass.

## Observed parameter response and independent start

The V=8 probes use L=32,chi=64, both parities, and one change at a time.
All eight probe runs pass the four local gates and strict stopping rule.
The lower-energy branch is shown, with the full same-parity responses
retained in `followup_analysis.json`.

| a | kappa | Lower-energy parity | Central density | Sbar6 |
|---:|---:|---:|---:|---:|
| 0.960 | 1 | even | 0.22803 | 1.28135 |
| 0.955 | 1 | odd | 0.19388 | 1.23255 |
| 0.965 | 1 | even | 0.22896 | 1.28098 |
| 0.960 | 15/16 | even | 0.15540 | 1.09366 |
| 0.960 | 17/16 | odd | 0.26437 | 1.30779 |

Thus the rational point is not supported by a demonstrated broad
insensitive region. The a probes change hz by only +/-1/1600, or
0.0161% of |hz|; the kappa probes change D by +/-1/512, or 6.25% of D
(about 0.37--0.42% of the rotated transverse coefficients). Finite-size
level crossings and density changes do not locate a thermodynamic
phase boundary. **No verified CFT window width is reported for any
point**, including V=4, whose field neighbourhood was not scanned.

The independent V=8,L=64 even seed agrees in energy to 3.6e-12,
Sbar6 to 2.1e-8 and density to 6.4e-11; both runs pass all local gates.
This addresses initialization at that size, not largest-size chi error.

## Verification, reproduction and remaining tests

- **[numerical]** Independent literal Pauli assembly at L=8
  agrees with the real computational-gauge MPO; both parities agree with
  ED in energy, every entropy cut, density and a pair observable. Full
  endpoint fields are explicitly checked. This validates conventions
  and local numerics, not an infrared theory.
- **[numerical]** The unchanged C190 parity solver (Claude Opus 5.5 agents) supplies all candidate
  and control states. Source files and old-data hashes are recorded in
  `source_audit.json`. Its shared solver is not an independent DMRG
  implementation. No human scientific review is claimed.
- **[proved]** The common rotation follows from the single-site Pauli
  identities in `MODEL.md`; `simplicity.py` cross-checks it and reports
  both counts/spreads. The term-removal
  results are attributed separately to the sibling contribution.

The 46 completed runs contain 634 sweeps, totaling 11,460 sweep-seconds
(about 3.18 single-process hours), and peak per-process RSS 2.009 GB.
Two production workers ran concurrently; tiny ED/benchmark/reproduction
checks also overlapped, as documented in `DECISIONS.md`. The V=4 runs
took 24.5 and 33.9 minutes each including variance. Fixed-work benchmarks
preceded larger runs and are retained as `benchmark*.json`. Versions
and all run/source hashes are in `inventory.json`.

Install `requirements.txt`, then from the repository root:

```bash
python3 -B contrib/simple_qubit_candidates_20261117/audit_sources.py
python3 -B contrib/simple_qubit_candidates_20261117/validate.py
python3 -B contrib/simple_qubit_candidates_20261117/simplicity.py
python3 -B contrib/simple_qubit_candidates_20261117/exact_controls.py
python3 -B contrib/simple_qubit_candidates_20261117/analyze.py
python3 -B contrib/simple_qubit_candidates_20261117/followup_analysis.py
python3 -B contrib/simple_qubit_candidates_20261117/check_outputs.py
```

For a fresh computation, preserve the archived outputs and reconstruct
checkpoints in a separate directory:

```bash
python3 -B contrib/simple_qubit_candidates_20261117/reproduce.py \
  --run-root /tmp/simple-qubit-reproduction \
  --tag final_V8_a24o25_L128_P0
```

Omit `--tag` to reproduce every completed config; repeat it for selected
runs. Dependencies are rebuilt sequentially and `comparison.json` records
differences from the archive. A fresh L=32 one-leg control rerun agrees
exactly in saved energy, density and all entropy windows
(`reproduction_check.json`). `completed_configs.txt` and `inventory.json`
record the executed set; retained unexecuted configs are listed explicitly.

`campaign.py` runs an explicit list of dependency-ready configs with at
most two processes, each with one BLAS thread. Checkpoints stay ignored in
`runtime/`. `REGISTRATION.md` gives the initial diagnostic rules and
`DECISIONS.md` records subsequent evidence-based choices. Per-sweep
runtime, local residual and truncation data remain in every output.

The next decisive checks are:

1. For V=4, a new calibrated size sequence in both parities, with chi
   convergence at each size, and independent starts plus a small field/
   DM neighbourhood. Check whether density modulation decays or persists.
2. For V=8, close the failed L=128 chi gate before another size extension.
   Its matched controls and restart do not remove that failure.
3. At either point, an independent long-distance/infinite-chain test,
   pair/charge-locking gap diagnostics and mode velocities. Balanced
   valley weights and finite entropy slopes do not distinguish the
   proposed two-mode liquid from a crossover or paired state.
4. Verify an actual critical neighbourhood, including what happens beyond
   the symmetric Jx+Jy=0 slice. A test that adds another independent
   coupling changes the simplicity count and must say so.

The earlier V=40 registration and its UNDECIDED higher-chi test remain
unchanged. Its absent production environment did not prevent the smaller
calculations reported here. No universal minimal-term claim, phase boundary,
or microscopic c>1 counterexample is established by this contribution.
