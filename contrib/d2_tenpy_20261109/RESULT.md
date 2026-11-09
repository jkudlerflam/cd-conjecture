# Full-spin DMRG: transverse-DM and simpler parity-family tests

OpenAI Codex, GPT-6 family, 2026-11-09. Interim numerical report; larger
checks continue. No d=2 counterexample is established. The strongest new
positive signal is the simpler parity family below, not the transverse-DM
point initially prioritized.

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
| Parity family, phi=1.53 (lowest tested parity at each L) | 1.35323214 | 1.60105887 | 2.14523 |

The two-leg control is a range-two calibration, not a d=2 NN candidate.
Its largest discarded weight is 1.9e-6 at chi=96, so it calibrates the
large separation of the two readings, not high-precision c. Candidate
discarded weights are much smaller (below 1e-10 in the final sweeps).
The C183 L=64 energy changes alternate by about 4e-10 and entropy changes
are below 5e-8; the simple point has changes around 3e-9 and 2e-7. These
are useful numerical tolerances, not rigorous eigenstate error bounds.
The exact free-fermion controls give .946721 and 2.047422 for this size
pair. Thus the two-leg chi=96 result has a downward bias of about .106
in the estimate; the exact comparison is retained in summary.json.

The short-chain c>1 profile fits are not stable: for C183 the central
half fit falls from 2.25 at L=32 to 1.15 at L=64. For the rational
point it falls from 1.99 to 1.01. The two-size readings are close to a
one-mode liquid at these sizes. This does not prove asymptotic c=1,
exclude a narrow nearby c>1 window, or distinguish an eventual gap.
The bulk densities (middle half) change from .1826 to .1716 (C183) and
.2103 to .1956 (rational), so density and boundary effects remain.

At L=128, C183 gives Sbar=.57013151 and c(64 to 128)=1.10032; the
central-half profile fit is .98583. The largest final discarded weight
is about 3e-10, but alternating sweep energy changes remain 5.4e-8 and
mean-entropy changes about 4.8e-6, so this is a provisional finite-size
reading, not a fully certified convergence limit. A different dense
initialization at L=64, chi=64 reproduces the chi=128 energy within
6.8e-8 and central entropy within .000195. The exact corresponding
64-to-128 control estimates are .97227 and 2.02543.

The original V=40 point quoted in the chat is much less entangled at
L=32 (Sbar=.04685, middle-half density .02916); it needs a field scan
before it can be treated as a finite-density critical target.

## Encouraging simpler full-microscopic point

Set V=40, Delta=1, w=exp(1.53 i)/80, Dy=0, and hz=-19.97589125:
Jx=(1+Re w)/2, Jy=(-1+Re w)/2, Jz=10, Dz=Im w/2. This is the
ordinary parity family with a small real hop, and contains one fewer
exchange term than the transverse-DM proposal.

Both exact parity sectors were optimized separately; small-L checks of
the selector agree with independently restricted ED. At L=32 the odd
sector has lower energy, and at L=64 the even sector does. The energy
separations are 8.05e-5 and 4.91e-5, respectively. The selected states
give c(32 to 64)=2.14523, and the L=64 profile fits across the four
preselected windows are 2.3985, 2.2288, 2.0504 and 2.3616. This is new
finite-size evidence from the full microscopic chain, not substituted
effective-model data.

Important limitations: the bulk density drifts .25494 to .21826;
the selected parity changes; the L=64 last energy changes are of order
5e-7 and mean-entropy changes about 1e-5 (odd sector 8e-7 and 1e-4).
The final discarded weights are below 2e-9, but small discarded weight
alone does not certify local eigensolver convergence or the global
ground state. No extrapolated c, gapless theorem, or counterexample
claim follows. The L=128 tests and independent solver tests of a nearby
rational bond are the next checks. The independent effective-model
evidence and its assumptions are audited in ../d2_priority_audit_20261109/.

An initial infinite-chain attempt failed numerical convergence and
violated the all-up variational upper bound. Its correlation length
is unusable. It is being independently diagnosed; no infinite-chain
claim follows from that run.

The rational Hamiltonian and independent small-ring comparisons are
in ../d2_simplification_20261109/. A separate exact finite-ring series
shows the surviving antiunitary allows a very small uniform Y source;
it does not establish relevance or a gap. See the corresponding audit.
