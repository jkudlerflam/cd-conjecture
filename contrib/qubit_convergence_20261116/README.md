# Microscopic qubit convergence audit and registered continuation

**Result: UNDECIDED for the new convergence test. No chi > 2048 production
state was computed.** The saved data still support a finite-size c near 2,
but do not close the bond-dimension gap or establish a counterexample.
The locally feasible audits, exact controls, Hamiltonian checks, checkpoint
authentication, interruption/resume tests, configuration path and bounded
benchmark are complete. [assessment.json](assessment.json) is the reproducible
machine-readable verdict; [RUNBOOK.md](RUNBOOK.md) is the executable next step.

## Claim and scope

[numerical] For the full microscopic d=2 Hermitian nearest-neighbour chain
with Jx=1/2, Jy=-1/2, D=1/160, Jz=10, the saved open-chain calculation uses
**hz=-19.976047250613814 on every site**, including both endpoints, held fixed
with size. The displayed field -19.976047 is rounded. The computational
staggered onsite gauge does not alter the physical one-site translation
symmetry. Only parity is exact; particle number is not conserved.

Here c is the total nonchiral Virasoro central charge of the complete unitary
conformal bulk theory, in the usual convention that one XX liquid has c=1.
The measured quantity is the finite-size estimator c_two, not c itself.
Existence of a unitary conformal infrared theory and survival of the second
gapless mode remain assumptions to test.

The immutable [registration](REGISTRATION_20261116.md) was committed and
pushed as **a010e9d5 before new solver outputs**. SHA-256:
`671d219d6a009947a82032984b011b06731bed4a279b2af9b615cdd715761f4f`.
It selects chi=2560, both parities, numerical branch checks and matched
controls. Unlike historical C195 T2, this new test requires C4 to pass.

## Ground-branch arithmetic from raw saved data

`audit.py` reads original JSON, hashes each source, retains canonical final
and last-sweep entropies separately, compares parity sectors by energy and
retains all parity pairings. C4 compares consecutive stages of the **same
initialization/continuation**, not different random/product starts.

| Sizes | Ground parities | Highest chi | c_two, W=6 | c_two, W=12 |
|---|---|---|---:|---:|
| 64 -> 128 | even -> odd | 384 -> 512 | 2.268064 | 2.356424 |
| 128 -> 256 | odd -> odd | 512 -> 1536 | 2.025339 | 2.058874 |
| 256 -> 512 | odd -> odd | 1536 -> 2048 | 1.982737 | 1.993043 |

These reproduce the later C190/C195 evidence. The 2026-11-10 audit's lower
chi verdict is superseded, while its rule to select parity by energy remains
essential. The final canonical Sbar6 values are 1.831004633 at L=128,
2.064980917 at L=256 and 2.294035715 at L=512. The last-sweep values differ
slightly; both are stored in [audit.json](audit.json).

| L, parity, chi | T1 residual/stationarity/truncation | Last chi step in Sbar6 | C4 <=0.01 |
|---|---|---:|---|
| 128, even, 512 | pass | 0.003358 | pass |
| 128, odd, 512 | pass | 0.003601 | pass |
| 256, even, 1024 | pass | 0.006232 | pass |
| 256, odd, 1536 | pass | 0.003076 | pass |
| 512, even, 1536 | pass | 0.031267 | **fail** |
| 512, odd, 1536 | pass | 0.031285 | **fail** |
| 512, odd, 2048 | pass | 0.010375 | **fail** |

At L=512, odd, chi=2048: maximum local residual 9.9994e-10, last |dE|
1.4250e-9, largest last-sweep |dSbar_W| 1.2935e-4, |dn_mid| 1.5754e-8,
discarded weight 4.9812e-11. Thus small residual and stationary energy are
already satisfied; **the remaining failure is the chi increment**, not
sweep stationarity. The density difference from the odd L=256 reference is
about +0.00214, within T4=0.006.

At matched chi=1536, E_even-E_odd at L=512 is +1.5682581e-6. Continuing only
odd to 2048 changes the comparison to +2.0338057e-6. The latter is not a
matched-accuracy comparison. At L=256 the corresponding splitting is
+6.2940749e-6 at matched chi=1024 and +6.4805026e-6 with odd at 1536.
Odd is the numerical ground branch throughout these comparisons, but neither
the local residual nor the variational energy supplies a rigorous bound on
the competing sector. The registration therefore also advances even.

No production MPS checkpoints exist in the local C190/C193/C195 directories.
The saved JSON are measurements, not resumable tensor states. Remote paths
in historical configurations are provenance, not verified present-day files.

## What the older headline establishes

REPORT 9.192-9.205 supports the arithmetic above. C195 explicitly allowed a
C4 failure without voiding its older T2 entropy-threshold pass. Consequently
STATUS's “passing the same test” does **not** mean that L=512 passed chi
convergence. That distinction reconciles the headline with the detailed data.

[heuristic] Linear extrapolation of the L=512 entropy in discarded weight
through chi=1536,2048, with the raw converged L=256 reference held fixed,
gives c_two=2.041091 at W=6 and 2.051359 at W=12. These are retained beside
the raw values in `assessment.json`, not used to pass the test. This differs
from extrapolating both sizes, and the choice of fit is not an error bound.

[proved, elementary] Variational optimization orders energies, not
entanglement entropies. For example H=I-|00><00| has a product ground state;
cos(theta)|00>+sin(theta)|11> has energy sin²(theta) arbitrarily near zero
and strictly greater entropy. No entropy lower-bound statement follows from
an energy upper bound. An observed increasing chi sequence likewise does not
certify monotonicity to infinite chi or bound the two-size difference. The
“lower bounds” language in REPORT 9.192/9.198/9.203 must be read as a trend,
not a theorem. C4 itself tests one increment, not a certified remainder.

[read] The open-chain entropy coefficient c/6 follows **assuming** a CFT
continuum limit; see Calabrese and Cardy, Eq. (2) and Sec. III.2 of
[Entanglement Entropy and Quantum Field Theory](https://arxiv.org/html/hep-th/0405152v3).
It does not establish that assumption for this bond.

## Controls and local validation

`controls.py` independently builds open-chain sine-wave orbitals and computes
entropies from the restricted correlation-matrix eigenvalues. For two XX
legs it sums the two leg entropies across the physical cut; it does not
mistake a spin next-neighbour hopping term for a single free JW chain.
The Gaussian correlation-matrix method is [read] from Peschel,
[Calculation of reduced density matrices from correlation functions](https://arxiv.org/html/cond-mat/0212631v1),
Eqs. (8)-(12). This formula is used only on the free controls.

| Saved control, 128 -> 256 | DMRG c_two W=6 | Exact finite-size c_two W=6 | Largest final Sbar error at 256 |
|---|---:|---:|---:|
| one XX, C190 solver family | 0.9736015 | 0.9736018 | 3.60e-8 |
| two XX, C190 solver family | 1.8983922 | 1.8990499 | 1.2841e-4 |

The control finite-size drift is real, not mainly solver error. All raw
W=1,6,12 estimates are in [controls.json](controls.json). C198's L=512 c2
control uses `fibdmrg` and N=128 (density 1/4), while the registered matched
control here has N=100 (density near 0.1957) and the C195 solver. C198 is useful
corroboration but does not meet the new matched-solver/filling gate.

New exact 256->512 calibration at the registered control fillings gives
(W=6,W=12)=(0.997509,0.997972) for one XX, and (2.025926,1.949638) for two
XX. Configurations for their C195 runs are saved in `control_configs/`;
**those large control DMRG runs have not been run**.

New local validation, actual C195 solver and driver, 145 s total on this Mac:

- Independent literal-Pauli Hamiltonian vs production MPO in the staggered
  gauge: maximum matrix discrepancy 2.0e-13 at L=8,10. Both parity ground
  energies agree with ED to 6.4e-13. Dropping endpoint half-fields gives an
  operator error of 19.97604725, so it is emphatically not an innocuous split.
- L=16, both parities: uninterrupted in-memory vs interrupted/resumed
  offload bookkeeping agrees to 2.85e-14 in energy and 2.08e-12 in entropy.
  A separate authenticated chi=32->40 continuation also passes. Wrong-parity
  checkpoint metadata are rejected. CPU emulation validates bookkeeping;
  a new GPU/hardware run still requires the registered benchmark/gate.
- One-XX and two-XX controls at L=32,64 using the C195 driver: largest
  entropy error 0.000415 and largest energy error 6.85e-6 (the two-leg
  L=64 chi=128 run). These are small-chain validation, not production controls.
- The real-gauge central-update benchmark runs at L=32,chi=64: three
  100-matvec updates take 0.171-0.188 s each. The production-sized invocation
  with a declared 48 GB budget rejects before allocation. No large random
  state or large-chi candidate state was constructed locally.

`validation.json` contains measured tolerances. `local_outputs/` preserves
the small raw run records; the large tensor files are ignored. No independent
human or second-model review has occurred.

## Correlation evidence and remaining infrared assumptions

The frozen C193 analysis was recomputed from its saved full correlators,
without changing fit ranges: [correlations.json](correlations.json). These
are **older chi<=768 states**, not the latest entropy checkpoints.
On the symmetric ground branch, K_s from the pair estimate is
2.4470, 2.0893, 1.9952 at L=64,128,256 (errors 0.315,0.282,0.204).
Its B2 inequalities pass. The structure-factor route gives
-0.6177, 3.4597, 1.7973 and triggers B3's accelerating-drift clause. The
negative small-size estimate and the previously documented failed control
calibration make that route unreliable there, but its registered B3 outcome
must remain visible. The later boundary-corrected fit is a distinct,
unregistered diagnostic, not a replacement of those frozen estimates.

The saved C195 48-site one-body window cannot cover even 0.10L=51.2 at
L=512, and its CP observable is not C193's connected Dpair. Therefore no
new 256->512 K_s/drift test is possible from these JSON alone. No new claim
that the second mode is gapless is made. C196 and C196/b's failed infinite
routes were reviewed; no retry or imposed exact U(1) was attempted.

## Resources and blocker

The read-only check found 48 GiB RAM, 18 CPU cores, no CUDA GPU and no local
production checkpoints. `ssh -o BatchMode=yes -o ConnectTimeout=5 typhon true`
timed out even with network access allowed. No remote job was submitted,
cancelled or modified. No fresh allocation/account access is assumed.

For **real float64** sym40, L=512, chi=2560, the dense stored tensors and
environments occupy approximately 51.4+128.6=180.1 GB; 25% headroom gives
225.1 GB. Recommend 256 GB host memory, one 80 GB CUDA GPU after benchmarking,
and 300 GB scratch, running sectors sequentially. Chi=3072 would need
259.0 GB before headroom (323.8 GB with it), so 2560 is the useful first step.

Scaling the archived **real sym40** chi=2048 sweep times gives 1.48-2.13 h
per 2560 sweep and 4.4-17.0 h for 3-8 sweeps per sector. Existing final
measurements took 4.36 h at 2048; cubic scaling estimates 8.5 h at 2560.
Including both sectors, even's intermediate 2048 stage and the L=256 parity
reference gives roughly **34-66 h**, or **68-132 h with a factor-two margin**,
before the new large control campaign and full C193 measurements. This is an
extrapolation, not a fresh GPU benchmark. Full C193 at chi=2560 converts the
MPS to complex128: over 100 GB before correlator stacks, so its unchanged
GPU measurement cannot fit on an 80 GB card. Use CPU with a separate budget
or validate a streaming measurement path first.

## Reproducing and provenance

Author and checker: OpenAI Codex, GPT-6 family, 2026-11-16. No subagents used.
Status: numerical audit and validation; infrared interpretation heuristic.
Dependencies: Python 3, NumPy, SciPy, threadpoolctl. Tested with Python 3.12.3,
NumPy 2.5.5, SciPy 1.18.1. Scripts use stable basic array/linear-algebra APIs.

From the repository root, with the dependencies installed:

```sh
python3 contrib/qubit_convergence_20261116/audit.py
python3 contrib/qubit_convergence_20261116/controls.py
python3 contrib/qubit_convergence_20261116/correlations.py
python3 contrib/qubit_convergence_20261116/validate.py
python3 contrib/qubit_convergence_20261116/resources.py
python3 contrib/qubit_convergence_20261116/assessment.py
```

The audit takes under a second; controls/correlations seconds to tens of
seconds; validation about 2.5 minutes on this host. For a fresh interruption
test, move only this contribution's `runtime/` aside first; completed runs
are otherwise skipped by the production driver. Benchmark commands and the
exact checkpoint continuation procedure are in [RUNBOOK.md](RUNBOOK.md).
`source_manifest.json` hashes the production solver and measurement modules;
all original files remain unchanged. This supplements REPORT 9.192-9.205 and
does not upgrade the project's overall counterexample status. Top-level
summary edits from other active chats were preserved.
