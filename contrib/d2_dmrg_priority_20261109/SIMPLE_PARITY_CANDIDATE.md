# Simple rational qubit candidate

GPT-6 Codex contribution, 2026-11-09. Status: both parity sectors at L=32 and L=64 are complete. The intermediate-size entropy signal is reproduced, but the registered local residual target of 1e-8 failed. This is not an established counterexample.

The literal translation-invariant nearest-neighbour qubit Hamiltonian is

\[
H=\sum_j\left[10Z_jZ_{j+1}+\frac{X_jX_{j+1}-Y_jY_{j+1}}2
+\frac{X_jX_{j+1}+Y_jY_{j+1}}{4000}
+\frac{X_jY_{j+1}-Y_jX_{j+1}}{160}\right]
-\frac{199759}{10000}\sum_j Z_j.
\]

Thus Jx=2001/4000, Jy=-1999/4000, Jz=10, Dz=1/160 and hz=-19.9759. There is no transverse DM term. The on-site dimension is exactly two; no auxiliary local degree of freedom or projector is used. The uniform field is applied to every site of the open finite chains, including their endpoints. The nonzero real hopping Jx+Jy=1/2000 moves the effective band away from the exactly imaginary-hopping point. Its hopping phase is arctan(25), close to the coordinator's promising C181 phase 1.53.

An equivalent integer bond, after distributing half the on-site field to each adjacent bond, adding a scalar 10 and rescaling by 10000, is

\[
10000(h_{\rm bond}+10I_4)=
\begin{pmatrix}
241&0&0&10000\\
0&0&5+125i&0\\
0&5-125i&0&0\\
10000&0&0&399759
\end{pmatrix},
\]

in the basis up-up, up-down, down-up, down-down. Independent Pauli reconstruction agrees to 3.6e-15. This is equivalent up to an overall energy scale and scalar for the periodic/infinite chain. The uniform-field open-chain calculations additionally include hz(Z_1+Z_L)/2 relative to the sum of these half-field bonds.

The Hamiltonian conserves total Z parity. Both even and odd sectors are computed, with a measured parity expectation verifying the selector. The lower physical energy at each length determines the ground branch; switching sectors between lengths is retained. Energies reported below remove the parity selector using its measured expectation.

## Validation and convergence protocol

Independent sparse Pauli-matrix diagonalization at L=10 verifies both sector energies to 2e-11, entropy profiles to 2.3e-8 and local magnetization to 2e-8. Measured total parity differs from its exact value by less than 2e-15. The optional restart path is separately ED-validated at L=10.

The C183 two-site Davidson solver is independent of the coordinator's TeNPy calculations. Production uses bond dimension 96 at L=32 and 128 at L=64 with repeated final sweeps and local residual target 1e-8. Residuals, discarded weights and last-sweep energy movement are retained. A local eigensolver hitting its iteration cap is not by itself accepted as converged. Saved MPS checkpoints permit additional sweeps without rebuilding the state.

## Results

[numerical] The lower **computed** energy switches from odd at L=32 to even at L=64:

| L | parity | physical energy | Sbar | central raw density | Var(N) |
|---:|---|---:|---:|---:|---:|
| 32 | even | -330.0058806527 | 1.31228419 | 0.18410 | 0.03673 |
| 32 | odd, lower | -330.0059530758 | 1.35324023 | 0.21896 | 0.03246 |
| 64 | even, lower | -650.0418240611 | 1.59917937 | 0.22416 | 0.08762 |
| 64 | odd | -650.0417845183 | 1.58153798 | 0.20696 | 0.09429 |

The parity energy differences are 7.2423e-5 at L=32 and 3.9543e-5 at L=64. Direct physical-energy contractions agree with the reported energies within 2.3e-11, and measured total parity is within 1.4e-12 of its intended value.

The primary, central-three-cut entropy slope is **c_eff=2.12889** from L=32 to L=64. The L=64 even-sector nested profile slopes are 2.374, 2.051, 2.198 and 2.136 when excluding the outer eighth, quarter, third and three eighths. The L=32 narrow-window fits remain unstable. This reproduces the coordinator's intermediate-size signal with a different solver and exact rational coefficients.

Before the L=64 entropy output, a secondary fixed-fraction average was declared: integrate the piecewise-linear entropy profile over [3L/8,5L/8] and divide by L/4. It gives c_eff=2.24750; the analogous middle-half average gives 2.23144. Neither replaces the primary estimator. Exact low-density controls at the same 32-to-64 sizes give quarter-window slopes 1.1493 and 2.8654 for c=1 and c=2, and half-window slopes 0.9959 and 1.9631. Shell effects therefore still bias short-size estimates even after smoothing. All values are retained in `rational_summary.json`.

**The registered 1e-8 local residual target did not pass.** Final maximum residuals are about 1.1e-6 at L=32, 1.18e-5 for even L=64 and 2.16e-5 for odd L=64. The two L=64 runs end with positive last-sweep energy changes 6.77e-7 and 3.74e-7, and discarded weights 8.43e-9 and 6.39e-9. Each used two final chi128 sweeps and took approximately 24–25 minutes. No entropy convergence history is available for those final two sweeps, so strict entropy convergence at L=64 is not certified.

Two additional chi96 sweeps at L=32 changed odd Sbar by 5.6e-6 and even Sbar by -3.9e-6, with energy changes at most 2.2e-9. Nevertheless, 26/62 odd-sector and 28/62 even-sector local solves still hit the 400-iteration cap. A fixed central two-site-problem audit directly verifies the residual. Increasing its Davidson subspace 24→48→80 reduces the residual 1.10e-6→5.03e-7→3.19e-7 at the same iteration cap, while local eigenvalues move less than 2e-12. Small energy motion alone does not certify the local solves.

The raw particle counts of the selected branches are 9.00485 and 16.01359, with small nonzero number variance despite exact parity conservation. The four-site-smoothed density range over the central half falls from 37.0% to 22.4%; the central-half mean density changes from 0.25495 to 0.21823. The table's central density instead averages the central quarter, so both definitions are explicit. These observables include virtual dressing and do not independently identify quasiparticle density.

All scheduled runs and continuations are finished. No further optimization or parameter scan is included in this contribution.

## Limits

The primary entropy diagnostic is c_eff=6[Sbar(64)-Sbar(32)]/log(2), where Sbar averages the three central cuts. Nested-window profile fits are also reported. These finite sizes do not establish a thermodynamic central charge: low-density shell changes, boundary oscillations, density drift and a long crossover can bias the estimate. A positive result would motivate larger lengths and additional infrared diagnostics, not establish c>1 as a theorem.

Coordinator-reported independent TeNPy data at the nearby **unrounded** phase 1.53 point already weaken the c≈2 interpretation: L=64 has Sbar≈1.60106, whereas the completed L=128 odd branch has Sbar≈1.69192, giving a size slope≈0.7865, and its central profile slope is≈1.076. The L=128 calculations themselves are not strictly converged. The even branch has Sbar≈1.68632 and is higher in reported energy by 9.0e-6, so it does not restore a c≈2 reading either. These are different coefficients from the exact rational point above and are not a proof about it, but they are an important warning against interpreting a reproduced 32-to-64 excess as stable infrared c>1. Crossover versus incomplete optimization remains unresolved. The coordinator's raw outputs are [odd](../d2_tenpy_20261109/results/c181153_L128_odd.json) and [even](../d2_tenpy_20261109/results/c181153_L128_even.json).

The earlier transverse-DM retuning in this contribution did not survive its matched-size test: the final adjusted point gives c_eff=0.9713 from L=24 to 48. The rational parity-only candidate is a distinct, simpler direction.

## Reproduction

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 contrib/d2_dmrg_priority_20261109/rational_campaign.py --stage validate --workers 2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 contrib/d2_dmrg_priority_20261109/rational_campaign.py --stage production --workers 2
python3 contrib/d2_dmrg_priority_20261109/rational_analyse.py
```
