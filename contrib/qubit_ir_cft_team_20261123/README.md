# Fixed qubit Hamiltonian: audited first infrared investigation

**Best-supported verdict: a finite-scale near-two entropy candidate,
with a conditional two-boson effective description; the asymptotic
microscopic CFT is unresolved.** Two independently gapless modes, a
compactification lattice, radii, mixing and physical velocities have not
been identified. This investigation supplies new exact selection rules,
conditional stability constraints and independently reproduced small-ring
spectral evidence. It does not assign a torus, orbifold or named rational
CFT from the entropy slope.

The fixed translation-invariant, nearest-neighbour, Hermitian Pauli chain is

\[
H=\sum_j\left[\frac{19}{32}X_jY_{j+1}
 +\frac{13}{32}Y_jX_{j+1}+Z_jZ_{j+1}\right]
 -\frac{1783}{1000}\sum_j Z_j.
\]

Its physical on-site dimension is **d=2**. The field is counted exactly
once at every site, including both open ends. No parameter was changed
and no long DMRG or background production was launched.

The user's requested [prioritized to-do list](TODO.md) gives the precise
missing observables, discriminating outcomes, minimum next measurements,
data availability and resource limitations. [CLAIMS.md](CLAIMS.md) tracks
each claim's status. Full worker notes and independent AI checks are in
`theory/`, `evidence/` and `audit/`; [PROGRESS.md](PROGRESS.md) records the
team, serial computation and synchronization checkpoints.

## What the saved evidence establishes

| Observable at this exact target | Reproduced reading | Scope |
|---|---|---|
| Infinite-chain entropy vs longest correlation length, chi=96–256 | c_slope=1.97265 +/- 0.00986 | Rounded-table regression; longest saved target xi=483.2 |
| Open-chain central-window entropy, L=96,128,192,256,384 | c_slope=1.95401 +/- 0.07500 | Last-entropy-step filter <=0.02; L=512 excluded |
| One-body fitted exponent, L=384 | 0.3569(100) | Saved estimator output; raw target function unavailable |
| Connected separation-three pair fitted exponent, L=384 | 0.9193(225) | Fitted momentum and raw complex function unavailable |
| Raw-density structure-factor estimator, L=384 | 0.675(27) | Not independently matched to conserved emergent charge |

The errors shown are regression/estimator errors, not full systematic
uncertainties. OBC entropy still changes at L=384. The iMPS and OBC
methods share the same microscopic implementation. The historical
November 18 c=3/2 claim was withdrawn on November 19; it is not reinstated.
The saved target OBC odd-even comparison covers L=64,96,128 and has
positive differences at all three sizes. Sign-changing examples belong
to other parameter points. None of those boundary splittings identifies
a bulk primary or a second gapless mode.

No target HDF5 states or raw complex measurement arrays were recoverable
in the accessible checkout, including ignored files. The other models'
raw C193 arrays are not data for this Hamiltonian. The three zero `check`
values printed by the target exporter are literal placeholders, distinct
from its separately implemented small-state contraction test.

## New audited conclusions

1. **Exact mapping, without a perturbative phase transfer.** A uniform
   rotation gives `[(XX-YY)/2+(3/32)(XY-YX)+ZZ]-1.783 Z`.
   Thus V/Delta=4 and Delta/V=1/4. The formal second-order hopping is
   1/8, while the direct imaginary hopping magnitude is 3/16: their
   ratio is **3/2**. The equal-hop condition discussed for delta=1/8
   does not hold at this fixed target. No controlled dilute or
   all-orders expansion follows from the exact rotation.

2. **A concrete operator-identification correction.** Staggered
   distance-two chirality is odd under one-site translation. Smooth
   valley-density imbalance is even. They cannot be linearly identified
   in the unbroken two-valley dictionary. Instead the staggered chirality
   has the exact symmetry characters of one Hermitian quadrature of
   the pi inter-valley pair primary; the staggered literal range-two
   `XX-YY` bond has those of the other quadrature. This permits overlap;
   symmetry does not force its coefficient to be nonzero.

3. **The locking test can include mixing.** Assume two compact bosons,
   unbroken microscopic symmetries, generic incommensurate carriers,
   and correctly matched bulk one-body and undressed pi-pair primaries.
   Write `p=eta_pair`, `d_theta=4 eta_b-p`. The symmetry-reduced canonical
   Gaussian covariance leaves one free mixing parameter t, and every
   generic permitted nonderivative harmonic has dimension

   `x_(r,s)=2[p(r-t s)^2+s^2/d_theta]`, with integer r,s.

   At the quoted central values, the primitive quartet has
   **x4=1.8386** (propagated fit error 0.0450), while every s!=0
   harmonic has **x>=2/d_theta, approximately 3.9347**, independently of t. Higher pure quartets
   have x>=7.3544. Thus the primitive quartet is the only relevant
   generic harmonic under these assumptions; special resonances still
   require their own test. More generally any genuine pi primary in
   this bosonic lattice with exponent below one has an allowed relevant
   doubled vertex, which need not be the pure quartet.

   These are conditional algebraic results, **not a measured target
   locking dimension or gap**. The saved pair momentum is missing, the
   bulk primary assignment and convergence are unproved, and the target
   continuum coupling is unknown. None of the listed symmetries forces the quartet to
   vanish. A weak nonzero coupling with this dimension grows by only
   about 1.25 between lengths 96 and 384 in linearized RG, so a long
   crossover remains plausible; its scale is not determined.
   An explicit abstract Gaussian counterexample reproduces the same
   two scalar exponents with a neutral transfer harmonic while keeping
   every generic lock irrelevant. It is not a proposed microscopic
   realization; it proves why the missing operator assignment matters.

4. **Mixing alone cannot reconcile the two charge estimators.** Under
   the same canonical Gaussian assumptions, a pair exponent p and a
   unit-normalized total-density slope obey `K_density*p>=1`. The saved
   central values instead give **0.62053**. Even symmetry-allowed
   relative phase-gradient admixture does not evade the unit-normalization
   bound. At least one of density matching, primary identification,
   finite-window accuracy or asymptotic Gaussian scaling is inadequate.
   This refutes their simultaneous literal interpretation, not the
   existence of two bosons. The reported Ks=1.97 is an equality only in
   the unmixed convention; the mixed analysis retains a conditional
   lower bound on the backscattering dimension.

5. **A shared finite-ring pole is not a mode count.** New original-frame
   PBC exact diagonalization at L=8,10,12 finds that pair, staggered
   bond and staggered chirality probes reach the same first pi pole,
   with gaps 0.193866, 0.106237, 0.057932. At L=8 the auditor independently
   reconstructed the literal Pauli matrix, all parity/momentum sector
   minima and the probe residues, agreeing within 1e-14. Bond/chirality
   cross weight at that pole is imaginary, with unit normalized squared
   coherence. Opposite antiunitary characters do not require distinct
   modes. The separation-three pair's strongest static peak moves away
   from pi on L=10,12 although pi weight remains substantial. Small-ring
   ground sectors and neutral pole weights change with L; no conformal
   tower or physical velocity is fitted to these sizes.

## What remains unidentified, including the central-charge convention

The intended convention is the Virasoro c of a unitary compact CFT, not
an assigned nonunitary/noncompact c_eff. Here `c_slope` is an entropy
fit coefficient. The usual relation to c presupposes an appropriate
critical scaling limit; the finite-chain formula has coefficient c/6
for one open-chain cut. [Calabrese and Cardy](https://arxiv.org/abs/hep-th/0405152)
provide that CFT relation, not a proof that this Hamiltonian meets it.

If two Gaussian sectors are ultimately established, each contributes
Virasoro c=1 and the entropy coefficients add to two. The exact
symmetries allow signed velocities in pairs `(+v1,-v1,+v2,-v2)`, with
v1 and v2 unequal. That is not evidence for one common-velocity
spacetime CFT. Equal-time exponents cannot determine either absolute
speed. A deliberately loose microscopic local-response light-cone
ceiling is `v <= 153 e/16`, approximately **26.0** in physical lattice
units, from a nearest-neighbour commutator bound; there is no inferred
positive lower bound or useful velocity estimate.

The proposed bosonic integer charge/winding lattice is an explicit
conditional dictionary, not an established microscopic compactification.
Its covariance is not fixed: even the two matched exponents leave
mixing free, and the full response Hessian contains additional velocity
information. A torus, enhanced rational theory, orbifold, or alternative
critical content requires independent operator/sector data. A long
crossover to one mode or a gap remains possible. In particular,
neither exact asymptotic c=2 nor asymptotic criticality is proved here.

## Reproduction, authorship and relationship to prior work

The research workers and independent auditor are AI agents, not human
referees. Their actual identifiers and model settings are recorded in
[PROGRESS.md](PROGRESS.md). The coordinator requested GPT-6 Astra Max;
the manager selected `gpt-6.2-astra / max` for the three workers instead
of the specifically requested, available `gpt-6-astra / max`. This
selection error is disclosed; work was preserved without relabeling
models or rerunning calculations merely to change metadata.

Use Python with NumPy for the bounded new scripts. All many-body checks
were serialized, one worker/thread at a time. The main L=8,10,12 run took
about 0.15 seconds; the independent L<=8 dense check about 0.30 seconds.
Hard caps are 90 and 30 seconds respectively. Saved-table regressions
require no many-body simulation. From the repository root, run these
in sequence; the first two save their small JSON outputs and the
independent audit scripts print their comparisons:

```sh
python3 contrib/qubit_ir_cft_team_20261123/evidence/read_tables.py
python3 contrib/qubit_ir_cft_team_20261123/evidence/ring_spectra.py
python3 contrib/qubit_ir_cft_team_20261123/audit/check_ring.py
python3 contrib/qubit_ir_cft_team_20261123/audit/check_fits_metric.py
```

The local system Python lacks NumPy. The successful recorded reruns used
the already installed interpreter
`/Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python`
(Python 3.12.14, NumPy 2.3.5). Substitute that executable for `python3`
on this machine, or use any prepared Python environment with NumPy.
No package installation was needed. The manager reran all four checks
sequentially; the independent audit outputs matched their saved JSON
exactly. Internal Markdown links were checked and the largest file is
under 90 KB.

See [theory derivations](theory/THEORY.md), [evidence report](evidence/EVIDENCE.md),
[audit verdict](audit/AUDIT.md), and [verification record](verification.json).

This supplements REPORT sections 9.209b and 9.212–9.214 and preserves all
earlier contributions. It sharpens their operator and stability caveats;
it supplies no parameter transfer from the V/Delta=40 or p153c models.
Only this contribution is owned by the team. Proposed shared-summary
text is provided separately for the coordinating chat.
