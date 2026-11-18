# Connected quartet response at the fixed asymmetric qubit point

Date: 2026-11-18. Assignment A of the
[stability round](../research_coordination_20261116/STABILITY_ROUND_20261118.md).

## Claim and status

**The stability of a two-boson / Virasoro c=2 phase at p153c remains
undecided.** This contribution computes a new finite-density, connected,
four-single-leg response in the full leading interacting blockade model.
It also gives an exact source-projected bound on that response. It does
not identify a finite-ring coefficient with a thermodynamic cosine.

The physical model is the fixed Hermitian, one-site translation-invariant,
nearest-neighbour **d=2** chain

```
w = exp(1.53 i)/80
Jx = (1+Re w)/2, Jy = (-1+Re w)/2, D = Im w/2
H = sum [Jx XX + Jy YY + D(XY-YX) + 10 ZZ] - 19.97589125 sum Z.
```

The field is present on every open-chain site in the saved C195 records.
Our small rings use periodic boundaries and the same bulk coefficients.
Nothing is retuned with length. The microscopic model conserves parity,
not particle number. The five independent real coefficients become four
after removing a positive scale; there are six literal Pauli monomials
in this basis. This work does not propose a simpler Hamiltonian.

The new results are:

1. **[proved, finite-system identity]** For a nondegenerate grand ground
   state, a specified two-mode Legendre transform gives the connected
   quartet coefficient as a sixth ground-energy derivative with two
   uniform pair-source insertions and four single-boson legs. Its 90
   interacting-resolvent words have no folded-energy, Wick or anomalous
   external-leg subtraction in the chosen momentum projection. A
   source-projected Cauchy bound controls each word. See
   [MATCHING.md](MATCHING.md).
2. **[numerical]** At the fixed p153c chemical potential, all number and
   momentum sectors are searched before selecting each reference. The
   finite-density coefficients, in the precise normalization below, are:

   | L | N | Density | `g4,L/t` | Cauchy bound on `abs(g4,L)/t` |
   |---|---|---|---|---|
   | 20 | 3 | 0.150000 | -0.000309319086 | 0.001169701615 |
   | 24 | 5 | 0.208333 | -0.001233862897 | 0.002907939935 |

   Here `g4,L = L` times the static collective four-leg Legendre
   coefficient, with Fourier operators normalized by `1/sqrt(L)` and
   `t=1/80`. These are finite response values, not a converged local
   matching coefficient. Both attractions, the blockade and the new
   asymmetric separation-two pair source are retained. Dropping that
   source changes the coefficient by 1.93% at L=20 and 0.74% at L=24.
3. **[proved counterexample]** A connected response of a *composite pair*
   followed by pair-level amputation need not be a quartet interaction.
   It is nonzero in an exactly quadratic squeezed oscillator whose
   four-single-leg vertex is identically zero. The composite diagnostic
   is retained in [medium.json](medium.json), explicitly distinguished
   from the four-leg result in [four_leg.json](four_leg.json).

The substantial size and cutoff dependence prevents a thermodynamic
estimate. A hard spectral cutoff is defined and tested, but is not
mistaken for a local Wilson matching procedure. The finite leading model
also has no established uniform microscopic remainder at epsilon=1/40.
The missing dressed covariance, velocities, field and composite-operator
normalizations, matching scale, and momentum uncertainty remain explicit.

## How it was established and checked

- Author and AI assistance: OpenAI Codex, GPT-6 Astra (Max), analytic
  derivation and bounded one-worker calculations.
- Independent human or second-agent check: none. The numerical validation
  uses a different algorithm and representation within this contribution.
- `medium.py` retains the complete leading Hamiltonian, searches all
  constrained sectors, evaluates the 90-word coefficient and bound, and
  reconstructs the interference between the two source channels.
- `check_matching.py` independently builds the 2,207-state L=16 real-space
  Hamiltonian and its constrained sparse resolvent, then uses multivariate
  Rayleigh recursion. The sixth coefficient agrees with the momentum-word
  calculation to **8.5e-15 relative**; the largest resolvent residual is
  7.3e-15. The L=16 momenta are an algorithm control, not a matching point.
- A separate finite-source energy calculation validates the composite
  fourth derivative; halving the source reduces its error by about four.

These are double-precision checks, not interval-arithmetic certificates
or independent peer verification. The bound's algebra is exact; its
displayed evaluations have numerical roundoff.

## Reproducing

From the repository root, with [requirements.txt](requirements.txt):

```sh
python -B contrib/qubit_detuned_quartet_20261118/medium.py
python -B contrib/qubit_detuned_quartet_20261118/check_matching.py
```

The saved run took 7.3 seconds for the finite responses and 0.2 seconds for
the independent check. Each script enforces one numerical worker. The
largest sector searched has dimension 1,144. No MPS optimization,
production run, remote job or checkpoint transfer is performed.

## Relation to prior work

The previous [quartet contribution](../qubit_stable_two_modes_20261118/README.md)
treated finite level couplings and a special cancellation. Here the fixed
chemical potential selects a finite-density grand ground state and the
external object is a connected four-leg response. The corrected
[locking audit](../v4_locking_operators_20261118/README.md) supplies the
interacting Hamiltonian and perturbative sources. C179/C191's symmetric
or vacuum estimates are not imported as p153c medium parameters.

The companion [medium-response contribution](../qubit_finite_density_modes_20261118/README.md)
keeps the actual dressed inputs unknown. Its finite-window carriers have
a quartet mismatch near +0.0166 at L=512, opposite in sign to the bare
value -0.0108673; neither has a controlled thermodynamic uncertainty.
The separate [chirality-Ising evidence](../d2_xy_yx_chirality_imps_20261118/README.md)
concerns another Hamiltonian and possible boson-times-Ising operator
content. It supplies no normalization or c=2 verdict for this calculation.
Shared overview files remain with the coordinator.
