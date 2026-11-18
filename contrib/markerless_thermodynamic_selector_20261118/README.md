# Half-occupancy response of the original markerless mobile model

2026-11-18. The full thermodynamic selector remains open. This round
establishes new uniform control of one exact response channel, derives
the other channels as explicit determinant sums, and controls clock
retardation within a defined nonzero-J second-cumulant calculation.
It adds no physical states or terms to the original corrected d=10 bond.

## Results

- **[Proved]** At V=N, for every even N, the exact canonical uniform
  contact susceptibility obeys **0<=t gamma_N(0)<=3/4**. Its exact
  thermodynamic value is **(24-pi^2)/(18pi^3)=0.0253181633451665...**,
  derived using the XXZ energy and compared with an independent finite
  two-particle-hole calculation. No uniform charge gap is used.
- **[Proved formula; bounded numerical checks]** Every material-momentum
  contact response is an explicit positive Bethe determinant sum.
  A same-grid two-particle-hole truncation is valid at q=0 and fails at
  q!=0. Gap-coordinate resolvents agree with the determinant sums on
  complete N=4,6 spectra. The first commensurate current branch contributes
  O(N^-6), rather than an infrared divergence.
- **[Proved projected identity]** For the actual balanced and head-number
  patches inside unchanged length-10 bulk, the canonical q=0 cross term
  is retained and has a finite favorable limit. A weaker sufficient
  static condition needs the local susceptibility and only the ten
  background harmonics, rather than a supremum over all momenta.
- **[Proved for a specified dynamical kernel]** Keeping the clock mass
  J alpha delta in the exact second-cumulant denominator cannot increase
  its curvature above the J=0 source-projected bound. This holds at
  every finite volume and J>0. It does not control the higher connected
  cumulants at the physical coupling.
- **[Exact clock algebra; conditional acoustic consequence]** A new
  closed spectral-moment identity gives the first single-clock acoustic
  logarithm as **J^5 log(t/J)**, with explicit n=9 coefficients for
  all three replacements. Its curvature stays finite. This is not a
  claim about the complete microscopic remainder.

The all-momentum bound (22) from the previous contribution is still
neither established nor refuted. The bounded half-density checks show
max t gamma_N(q)=0.09627, 0.09978, 0.10138, 0.10221 at N=4,6,8,10.
The hydrodynamic small-nonzero-q prediction is 0.103893960837...,
distinct from the canonical q=0 result. These numbers are not a
certified supremum. No L-independent J/t selector allowance is inferred.

## Definitions, proofs, and scope

[SUSCEPTIBILITY.md](SUSCEPTIBILITY.md) contains the all-momentum contact
form factors, the uniform q=0 proof, the thermodynamic calculation, the
commensurate-branch test, and the exact projected comparison including
its unchanged-background cross term. It also bounds the thermodynamic
homogeneous scalar-reference curvature at nonzero source.

[RETARDATION.md](RETARDATION.md) defines the resummed second cumulant
before bounding it, proves the dipole and spectral-moment identities,
calculates its first acoustic logarithm, and states the missing
higher-connected relative remainder. It does not substitute a frozen
clock or a scalar impurity for the original quantum clocks.

The microscopic model remains Hermitian, nearest-neighbour and one-site
translation invariant at physical d=10, with L-independent coefficients.
The target is n=9, b=10, canonical V=N. The replacements are (10,10) to
(9,11), nine length-10 cells to ten length-9 cells, and eleven length-10
cells to ten length-11 cells, with arbitrary unchanged exterior bulk.
The biased all-n clock and L1 readout of the paired contribution are not
used. The prior finite certificates are preserved unchanged.

Zero-port selection, charged excitation gaps, finite ports, the logical
2^-35 allowance, and full infrared stability remain separate. The global
charge-field repair is usable only for the original conserved charge;
at the zero field considered here it contributes no new estimate. No
unbounded-c construction or unitary compact Virasoro CFT is established.

There is no new Hamiltonian, new continuous control, or improvement to
the original coefficient hierarchy. The auxiliary eta interpolation and
homogeneous scalar reference are calculation tools, not added physical
interactions. Thus the original term counts and simplicity limitations
remain those of the preceding contribution.

## Reproduction and review

Author and AI assistance: OpenAI Codex, GPT-6 Astra / Max assignment.
No human or independent agent review has been performed. No chats or
subagents were created. No DMRG, background production, or large scans
were run; all numerical-library workers are restricted to one.

```sh
python3 contrib/markerless_thermodynamic_selector_20261118/checks.py
```

Dependencies are numpy, scipy, and mpmath. The script reads the prior
contribution without editing it and writes only [checks.json](checks.json)
here. The largest free-bath matrix has 92,378 states; the nonzero-J
coefficient audit has 4,158 states. The all-q determinant audit uses
only N=4,6. The q=0 free-fermion sums run through N=80 without a many-body
matrix. The recorded check run takes approximately 0.46 seconds, with
the runtime and package versions saved in the JSON report. Numerical values are
floating-point checks; the uniform bound is the written analytic proof.

Primary sources actually opened are the
[phase-model Bethe equations](https://arxiv.org/pdf/solv-int/9612002) and
the [XXZ ground-state energy formula, equation (11)](https://arxiv.org/pdf/cond-mat/0304475).
Their conventions are converted explicitly in the proof. The new contact
determinant, bounds, and clock identities are derived in this folder.

The next decisive calculation is to bound the explicit determinant sums
for local contact response and period-10 harmonics uniformly in N,
then control the connected eta>=3 defect remainder at eta=1. The new
second-cumulant curvature result isolates, but does not solve, that task.
