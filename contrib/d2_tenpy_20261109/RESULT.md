# Full-spin DMRG: first calibrated readings

OpenAI Codex, GPT-6 family, 2026-11-09. Interim numerical report; larger
checks continue. No d=2 counterexample is established.

The literal microscopic spin chain is used with no blockade projection.
TeNPy's matrices and L=8 ground-state entropies agree with independent
Kronecker-product exact diagonalization (matrix error below 3e-14,
entropy below 2e-12). The C183 L=12 answer also matches the independently
run repository Davidson solver and sparse ED.

For the predetermined estimator c_2=6[Sbar(2L)-Sbar(L)]/ln2, where Sbar
averages three central cuts, the completed L=32 to 64 results are:

| Full microscopic point | Sbar(32) | Sbar(64) | c_2 |
|---|---:|---:|---:|
| C183 transverse DM, hz=-7.9463 | .32228975 | .44301798 | 1.04504 |
| Simple rational transverse DM, hz=-8137/1024 | .35902599 | .47706901 | 1.02180 |
| XX c=1 estimator control | .87804764 | .98741699 | .94672 |
| Two independent XX legs, c=2 estimator control | 1.48793688 | 1.71220998 | 1.94135 |

The two-leg control is a range-two calibration, not a d=2 NN candidate.
Its largest discarded weight is 1.9e-6 at chi=96, so it calibrates the
large separation of the two readings, not high-precision c. Candidate
discarded weights are much smaller (below 1e-10 in the final sweeps).
The C183 L=64 energy changes alternate by about 4e-10 and entropy changes
are below 5e-8; the simple point has changes around 3e-9 and 2e-7. These
are useful numerical tolerances, not rigorous eigenstate error bounds.

The short-chain c>1 profile fits are not stable: for C183 the central
quarter fit falls from 2.25 at L=32 to 1.15 at L=64. For the rational
point it falls from 1.99 to 1.01. The two-size readings are close to a
one-mode liquid at these sizes. This does not prove asymptotic c=1,
exclude a narrow nearby c>1 window, or distinguish an eventual gap.
The bulk densities (middle half) change from .1826 to .1716 (C183) and
.2103 to .1956 (rational), so density and boundary effects remain.

The original V=40 point quoted in the chat is much less entangled at
L=32 (Sbar=.04685, middle-half density .02916); it needs a field scan
before it can be treated as a finite-density critical target.

The full-chain C181 parity-family runs are a separate search arm.
Neither a large single-size entropy nor an unfinished sweep is a
central-charge measurement. Its independent effective-model evidence
is audited in ../d2_priority_audit_20261109/.

An initial infinite-chain attempt failed numerical convergence and
violated the all-up variational upper bound. Its correlation length
is unusable. It is being independently diagnosed; no infinite-chain
claim follows from that run.

The rational Hamiltonian and independent small-ring comparisons are
in ../d2_simplification_20261109/. A separate exact finite-ring series
shows the surviving antiunitary allows a very small uniform Y source;
it does not establish relevance or a gap. See the corresponding audit.
