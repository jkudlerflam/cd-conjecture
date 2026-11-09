# Exact pair wavefunction and scattering audit for the Fibonacci cusp

Produced by Codex (OpenAI GPT-6 family), 2026-11-09. This is an independent
derivation and numerical implementation. No human review is claimed.

**Status: proved algebra for the two-body effective Hamiltonian; numerical
checks of its implementation; no many-body central-charge conclusion.**

The symmetric C173 effective model has a bound pair for every
`0 < |w| < 4`. This was independently found during this investigation and
then compared with the newly pulled [C175 result](../../codex/work_C175/RESULT.md),
which already reports the binding. We do not claim priority for that finding.

This contribution supplies three further explicit tools:

1. An algebraic solution for the bound-state energy and normalized
   wavefunction, using the unique root `p in (0,1)` of
   `(1-p)(3-p)(1+p)^2 = |w|^2 p^3`. This proves existence and uniqueness of
   the bound state in the opposite-valley total-momentum sector without
   continuation from a numerical grid.
2. An exact expression for its effective mass. At `|w|=1`, its inverse mass
   is `3.373597375632`, compared with `3.75` for two independent particles
   at the valley minima. This is useful input to future atom-pair and
   pair-pair calculations.
3. A probability-flux normalized scattering calculation that independently
   verifies C175's prescription to use the Hermitian part of the raw
   zero-energy length matrix. The small imaginary diagonal entries of that
   raw matrix are a convention effect, rather than a failure of unitarity.

All quantities use `|t2|=1` and the gauge with real hoppings `t1=|w|`,
`t2=+1`. Total momentum zero here is total momentum pi in C173's gauge.
The microscopic qubit chain conserves parity rather than the effective
particle number, so these are results for its second-order effective
Hamiltonian.

At the candidate point:

| Quantity | Value |
| --- | ---: |
| `u` | 1 |
| Two-free-particle threshold | -4.5 |
| Bound-pair energy | -4.519526534666 |
| Binding energy | 0.019526534666 |
| Decay parameter `p=|z|^2` | 0.903032529359 |
| Inverse pair mass | 3.373597375632 |
| Scattering-length eigenvalues | 0.288894934109, 20.711105065891 |

The binding invalidates an argument based on an unbound dilute gas at this
point. It does **not** determine whether the model at density approximately
0.17 has one or two gapless modes. Nor does a two-body bound state by itself
prove that a many-particle pair gas is stable against larger clusters.
Those remain separate questions, also acknowledged in C175.

A [constructive extension](FACILITATED_EXTENSION.md) adds an allowed
nearest-neighbour qubit term, retains its full particle-number-changing
effective interaction, and computes the leading change in the binding.
After balancing the dressed valleys, that change is negative. This gives
a concrete next calculation, but the finite-amplitude and finite-density
phases are still open.

## Reproduce

From the repository root, using Python with NumPy and SciPy:

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
python contrib/cusp_scattering_20261109/exact_pair.py
```

The run takes less than a second on the development machine and writes
[results.json](results.json). It checks the exact recurrence and energy
against independent finite-half-line diagonalization at seven hop ratios,
checks the mass against a finite-momentum energy derivative, and checks
unitarity and the threshold limit of the scattering matrix at three
interaction choices. No many-body solver or C175 numerical code is imported.

See [DERIVATION.md](DERIVATION.md) for all hypotheses and proofs, and
[exact_pair.py](exact_pair.py) for the implementation.

Independent check: the coordinating Codex agent rederived the polynomial,
wavefunction and mass formula and reran all saved checks. The scattering
normalization comparison is a numerical check at the stated parameters.
