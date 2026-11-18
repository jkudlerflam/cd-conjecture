# Eleven-state loop: linear centre splitting needs finite express delay

2026-11-18. **[Proved analytic partial family; finite numerical checks;
all-n nonidentity still open.]** No microscopic CFT claim.

The same eleven physical states, nearest-neighbor rules and physical
C84 phases are retained. Only the already allowed closing amplitude
t_d is turned on. The inactive D_n and X_n endpoints remain in the
2n+3-vertex carrier. The previous contribution's biased clock and L1
readout are unchanged.

The exact derivative of the otherwise zero opposite-root pair sums is
derived in [DERIVATION.md](DERIVATION.md) and implemented in loop.py.
Hopping is even in t_d, while the imaginary
onsite term is odd. The first derivative vanishes if the tree hopping
is replaced by its strict symmetric Gram limit. At finite express
energy D the tree's antisymmetric delay term combines with the loop
onsite term to give a nonzero linear response in explicit families.
For eta=t_d D^(n+1), its leading size is eta/D.

The all-n question reduces to explicit rational moment coefficients L_j
in moments.py. The last is exactly -1/(16n-12). It gives two distinct
extreme centre responses for every n>=2, also separated from all the
remaining responses. This proves full distinctness at n=2,3. The exact
n=4 coefficient L_1=-5033/2227500 completes that case as well.

Thus sufficiently small nonzero loop amplitudes give **all pair-sum
margins positive for an analytic n=1,...,4 family**, with finite D and
the actual geometric physical field. This is an existence argument
using nonzero asymptotic coefficients, not a numerical extrapolation.
It does not certify a new large-n member or establish arbitrary n.
For n>=5, proving L_j !=0 for every 1<=j<=n-3 would complete this
sufficient all-n argument. Their observed negative signs have not been
proved generally. A zero would require examining higher coefficients;
it would not by itself prove an exact resonance.

The literal modification adds no states and no new rule type beyond
restoring the two existing loop transitions. One real coefficient
controls both; their phases remain fixed. [ERRORS_AND_SCOPE.md](ERRORS_AND_SCOPE.md)
distinguishes absolute root separation, non-Cooper comparisons and
the new Cooper margin, and budgets finite-g, finite-K and port errors
against the latter. In particular the lambda=0 finite-K splitting is
an absolute physical error. The hierarchy remains severe: even the
small diagnostic example has D/t_d=10^22 before the other small
couplings and large penalty are specified.

Nonzero-port ground-sector selection, volume-uniform matching,
interacting multichannel stability and conformal-observable control
remain separate assumptions. No total nonchiral Virasoro c of a
unitary compact microscopic infrared CFT is proved or measured here.

## Verification and reproduction

Author and AI assistance: OpenAI Codex, GPT-6 family, for derivation,
code and checks. No independent human or separate-agent review is
claimed. The original bond-to-scalar Fock-space match is retained from
the previous contribution; no new logical operator replaces that bond.

Run from the repository root with requirements.txt installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python contrib/eleven_state_loop_splitting_20261118/checks.py
```

[checks.json](checks.json) records six check groups: exact loop parity,
finite-D delay convergence, endpoint shortest walks, exact rational
moment coefficients at n=2,3,4, their eigenvalue asymptotics, and an
independent 75-digit full-loop root calculation at n=2. Runtime is
about one second, one worker. The largest pair solve has dimension 81.

At n=2, D=10^4,w=1/2,q=1/16,nu=10^-6, the derivatives with respect
to eta are approximately +0.18736054 and -0.18736054. At eta=10^-6,
the full pair-sum margin is approximately 3.74721e-7. Halving eta
reduces the derivative-approximation error by four, as the odd Taylor
expansion predicts. These numerical values are validation, not rational
interval certificates or a many-body phase measurement.

## Relation to the repository

The source search belongs to the coordinating chat. The GitHub clue
has not identified a ten-state arbitrary-c construction; this folder
does not repeat that search. All edits are confined to this new folder.
The original completed contributions and top-level status files are
preserved; the coordinating chat maintains the status index.
