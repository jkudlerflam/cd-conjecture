# Finite-density modes of the asymmetric p153c qubit chain

Date: 2026-11-18. Assignment B of the
[stability round](../research_coordination_20261116/STABILITY_ROUND_20261118.md).

## Claim and status

The fixed, Hermitian, translation-invariant nearest-neighbour microscopic
**d=2** Hamiltonian is

```
w = exp(1.53 i)/80
Jx = (1+Re w)/2, Jy = (-1+Re w)/2, D = Im w/2
H = sum [Jx XX + Jy YY + D(XY-YX) + 10 ZZ] - 19.97589125 sum Z.
```

The field is present on **every** open-chain site, as in C195; no length
dependent retuning is made. Parity is exact, particle number is not.
The family has five real coefficients, four after removing a common positive
scale, and six literal Pauli monomials including the field. This work does
not simplify those coefficients or propose a new Hamiltonian.

**The thermodynamic two-boson / Virasoro c=2 verdict remains undecided.**
Neither a full dressed response matrix nor a controlled instability is
determined by the saved data. Three concrete calculations advance the test:

1. **[proved, finite-order algebra]** The leading Schrieffer--Wolff dressing
   of the raw total spin particle count is
   `N_raw,eff = N + 2(Delta/V)^2 B + O((Delta/V)^3)`, with an explicit
   positive local operator satisfying `0 <= B <= L - 8N/9`. Thus raw number
   fluctuations contain a virtual-pair background even in a state with
   definite leading effective N. See [DRESSING.md](DRESSING.md). This is a
   perturbative coefficient identity, not a uniform all-orders error bound.
2. **[numerical, saved finite-chain records]** At the *actual* p153c point,
   the two filtered one-body carrier phases have `Q+pi` about **+0.00830**
   at L=512, whereas the bare band minima give **-0.00543367**. Connected
   long-pair records at L=256 independently give positive shifts of
   0.00924--0.00978. The short-window phase is stable with bond dimension,
   but drifts with L and lacks a controlled thermodynamic/matching error.
   It cannot yet be used as a rigorous quartet cutoff. See [DATA.md](DATA.md).
3. **[proved, conditional response identity; calibrated numerically]**
   A four-channel matched static response determines the complete Gaussian
   Hamiltonian, including density-current mixing and unequal chiral speeds.
   The leading family's usual parameter probes have rank at most **three**
   because one combination is H itself. A distance-two current supplies a
   fourth independent probe in an exact asymmetric free-band calibration;
   its rank and normalization in the interacting medium remain to be
   measured. See [RESPONSE.md](RESPONSE.md).

The [dated measurement protocol](PROTOCOL.md) identifies the minimum missing
observables, tests two independent bosons rather than two peaks, and includes
the complete asymmetric parity-even operator lattice. [interface.json](interface.json)
keeps unknown densities, covariance, velocities, matching scale and their
errors explicitly null. No p153c exponent, K matrix or vacuum-scattering
substitute has been invented. The separate connected-amplitude calculation
belongs to [assignment A](../qubit_detuned_quartet_20261118/).

## How this was established and checked

- Author and AI assistance: OpenAI Codex, GPT-6 Astra (Max), analytic
  derivations, read-only record analysis and bounded one-worker calculations.
- Independent human or second-agent check: none. Internal checks use an
  integer finite-Hilbert-space calculation and the earlier C191 chiral
  implementation, not a second claim of peer verification.
- `dressing.py`: the integer operator identity is exact on periodic rings
  L=8,9,10,12 (largest projected space 322); finite eigenvalues supplement
  the analytic bound. Runtime 0.05 s.
- `records.py`: hashes and audits 13 standard records and six long-correlation
  records, preserving complex phases and subtracting measured disconnected
  pair means. Runtime 0.28 s. The repeated L=64 even-sector initializations
  are recorded separately, not treated as independent statistical samples.
- `calibrate.py`: exact free two-sea responses, the exact hard-rod equation
  of state with an integer canonical coordinate change, C191 metric
  agreement, exhaustive lattice checks and exact XX-pair phase controls.
  Runtime 0.04 s. These controls are **not** realizations or fitted
  approximations of the interacting target.

The response inversion recovers the asymmetric free H to 1.4e-15. On the
interacting hard-rod control, the smallest finite-difference step gives
max errors 3.1e-7 in H and 1.3e-6 in G. The covariance agrees with C191 to
1.1e-14. A deliberately unresolved nearby carrier produces a phase bias
0.00370 despite a regression error only 0.0000046, demonstrating why that
fit error cannot serve as the target's uncertainty bar.

## Reproducing

From the repository root, with numpy, scipy and threadpoolctl installed:

```sh
python -B contrib/qubit_finite_density_modes_20261118/dressing.py
python -B contrib/qubit_finite_density_modes_20261118/records.py
python -B contrib/qubit_finite_density_modes_20261118/calibrate.py
```

Each script enforces one numerical worker and writes its small JSON result
in this folder. There is no DMRG, remote submission, checkpoint download or
background production. `response.py` is the reusable conditional extraction
and lattice interface. Dependencies beyond numpy/scipy are listed in
[requirements.txt](requirements.txt).

## Relation to existing work

C191 supplies the correct spatial chiral metric and essential frame/boundary
cautions. Its symmetric vacuum-scattering estimates do not determine this
asymmetric finite-density medium. C190/C193/C195 supply the actual records;
the corrected locking audit supplies the bosonic compactification, which is
rederived off symmetry here. The earlier three-mode deformation work is
preserved unchanged. Shared overview files are left to the coordinator.

The new, separate [chirality-Ising numerical contribution](../d2_xy_yx_chirality_imps_20261118/README.md)
concerns another fixed Hamiltonian and a possible boson times Ising theory.
Its c=3/2 evidence is not a response measurement for p153c. It reinforces
the need to identify the operator content rather than infer two bosons from
an entropy slope alone.
