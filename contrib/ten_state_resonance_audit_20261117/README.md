# A nine-state resonant seed: the pair-sum certificate is not necessary

## Claim and status

**[Proved finite algebra; numerical matching; perturbative RG interpretation.]**
A literal Hermitian one-site-invariant nearest-neighbour bond on nine
physical states has a frozen two-flavour code whose leading Bloch symbol
has four simple Fermi points. Their two opposite-momentum pairs have
exactly the same total momentum, so C84's strict pair-sum certificate
fails. The actual generated quartic vertex has nonzero pair transfer,
but the leading weak-coupling RG flow is in the marginally irrelevant
basin, with `D-Y=.51894579 g²+O(g⁴)>0`. Thus reciprocity alone does
not supply an instability theorem. C86/C87 themselves correctly state
restricted certificate obstructions.

This is **not** an unbounded-c construction or a microscopic CFT theorem.
The candidate infrared theory in the selected frozen sector is c=2,
well below d-1=8. Finite-K splitting, nonzero-port ground-sector selection,
volume-uniform matching, and an arbitrary-flavour extension remain open.
The model adds explicit physical background fields and bit swaps to the
single rail; it does not assert that the unmodified single rail already
has two gapless flavours. Its hierarchy of microscopic coefficients
remains substantial.

## Derivation and scope

- [SOURCE_AUDIT.md](SOURCE_AUDIT.md): the actual C84/C86/C87 hypotheses.
- [BOND.md](BOND.md): every physical state and local rule; exact code projections.
- [SYMBOL.md](SYMBOL.md): an exact rational/algebraic four-root identity.
- [INTERACTIONS.md](INTERACTIONS.md): the generated quartic, band matching and RG signs.
- [LIMITS.md](LIMITS.md): finite-K, ground-sector, infrared and arbitrary-flavour gaps.

## Method and checks

The nine states are `L0,L1,R0,R1,DA,DB,E,MA,MB`. All three marker
states and both carrier states are counted. There is no separate
vacancy; the neutral head shares R0. The n-flavour legal cell occupies
n+2 physical sites and has code dimension 2^n. Markers are frozen.

`literal.py` builds the physical two-site rules and independently reduces
their legal cell matrices. `mechanism.py` derives the common-carrier
resolvent coefficients through fourth order in conversion. Their full
Fock-space operators agree for n=2 and 3, including the generated
quartic interaction and the physical counterterm projections.
`audit.py` checks charges, marker conservation, roots, all pair sums,
forward scattering and pair transfer. No effective interaction was
chosen independently of the literal bond.
The largest floating discrepancy in the full matching/projection checks
is 4.5e-16. `exact_symbol.py` independently certifies the rational inverse,
indefinite hopping signs and four simple roots at cosine values -7/15
and 13/15. Actual velocities and vertices use the orthonormal physical
code, not the congruence used to prove the zeros.

Author and AI assistance: OpenAI Codex, GPT-6, for derivation, code and
literature audit. No independent human or separate-agent review is claimed.
The two implementations provide an internal cross-check, not a referee.

## Reproduction

Run from the repository root with numpy and scipy, one worker:

```sh
OPENBLAS_NUM_THREADS=1 python contrib/ten_state_resonance_audit_20261117/exact_symbol.py
OPENBLAS_NUM_THREADS=1 python contrib/ten_state_resonance_audit_20261117/audit.py
```

The largest pair resolvent has dimension 6400 and is applied to 64 code
vectors. The complete audit takes a few seconds or less on this machine;
the recorded runtime is in `audit.json`. No DMRG, cluster job or broad
scan of previously excluded mechanisms is used. The complete audit took
about .1 seconds here, after its two-flavour benchmark. Saved outputs are
`exact_symbol.json` and `audit.json`, including the complete quartic
vertices, neighbouring points and finite-conversion remainder checks.
See [PROTOCOL.md](PROTOCOL.md) for the bounded scope.
