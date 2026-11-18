# Smallest-d audit and the eleven-state all-n transfer test

Started 2026-11-17; completed 2026-11-18.

## Claim and status

**[Proved here; finite numerical cross-checks; not independently reviewed.]**
An explicit nearest-neighbor, one-site translation-invariant, Hermitian,
U(1)-conserving modification of C84 uses the same **eleven physical
states**, including all markers. For every finite payload size n it has:

- A biased legal vacuum code of dimension 2^n, on n+2 physical sites,
  with a uniform auxiliary gap. The bias is an actual boundary potential.
- Full ground-length selection against all marker patterns at zero port,
  with a finite grammar penalty. The physical payload band is split at
  finite penalty; the error is retained in absolute energy units.
- A positive all-n Gram limit that retains the inactive carrier endpoint,
  and an analytic margin for 2n simple roots and all non-Cooper pair-sum
  comparisons in its leading tree symbol.

The tree has **exactly zero full pair-sum margin for n>=2** because
opposite Fermi points sum to zero in every band. The complete two-bus
nonresonance proof therefore does not transfer to this one-bus model.
This is an obstruction to that sufficient proof, not an instability
theorem or an eleven-state no-go. The genuine C84 loop can break this
reciprocity; its all-n nonidentity/margin problem remains open.

**No microscopic central charge is established here.** Physical c means
the total nonchiral Virasoro charge of a unitary compact infrared CFT.
Nonzero-port ground-sector selection, volume-uniform matching and
interacting infrared stability are still required, even for the original
certified finite member. The new couplings do not inherit its certificate.

## Sources and priorities

The reported d=10 arbitrary-c nearest-neighbor construction is **not yet
located**, not disproved. [SOURCE_LEDGER.md](SOURCE_LEDGER.md) lists the
actual d=9–14 records, every physical state, range, n dependence,
certified member, ground-sector qualification and infrared assumption.
C84 certifies n=11; the same eleven-state material family also has a
regenerated n=13 certificate in C85/C88. Neither is an all-n theorem.
The [complementary nine-state audit](../ten_state_resonance_audit_20261117/README.md)
supplies a resonant two-flavour seed and keeps the smaller-d route open;
a one-hole d=10 mobile lift is not automatic.

[PROOF.md](PROOF.md) gives the literal modification, the full finite-K
selector bound, the Gram transfer and the exact remaining resonance.
[SIMPLICITY.md](SIMPLICITY.md) separates independent couplings from
literal monomials and records the severe tuning and coefficient hierarchy.
The smaller state count does not make this a simple Hamiltonian.

## Establishment and reproduction

Author and AI assistance: OpenAI Codex, GPT-6 family, for derivation,
source audit, implementation and checks. No human or separate-agent
review is claimed. The literal bond and full Fock-space reductions
provide internal cross-checks. Old large certificates were inspected
and hashed, not rerun or newly independently certified.

From the repository root, with the dependencies in requirements.txt:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python contrib/minimal_d_unbounded_audit_20261117/checks.py
```

The script saves [checks.json](checks.json). It takes about one second
on this machine, using one worker. It checks exact recovery of the old
literal bond; symmetry and term counts; all payload clocks at n<=3;
finite-penalty forms through four data sites; actual payload splitting;
the entire two-cell Fock coefficient at n=1,2; inactive-endpoint effects;
large-D convergence; and high-precision simple-tree margins at n<=3.
The largest full cell matrix has dimension 12,288. The largest relative
onsite matching discrepancy is 2.8e-9, from cancellation; the hopping
match is at about 3e-15 and the Gaussian residual below 3e-16.

The default floating parameters in model.py are diagnostic examples.
For an all-n schedule the small clock energy and chemical potential
must be evaluated with adequate precision; double precision can lose
the selector offset. The analytic proof is independent of these defaults.

## Relation to the repository

This narrows the gap between the larger all-n construction and the
smallest fully audited finite-member alphabet without adding states.
It also explains why a failed fixed-D schedule or coarse root enclosure
does not establish a dimension obstruction. All changes in this round
are confined to this contribution; the coordinating chat handles the
top-level status documents. No new mobile extension is constructed.
