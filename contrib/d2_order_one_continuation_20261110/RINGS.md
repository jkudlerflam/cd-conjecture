# Independent periodic-ring check

Both attractive models were solved in the full microscopic qubit Hilbert space,
separately in both parity sectors. Translation orbits of returned ground
vectors were completed before selecting entropies. The lowest momentum
index in the ground manifold is the fixed selection rule. All selected
ground states here have even parity and momentum zero.

| Point | Sizes | Half-chain size slope | Profile slopes |
|---|---:|---:|---:|
| pref_rat | 14→16 | 0.951614 | 0.982112→0.975249 |
| pref_rat | 16→20 | 0.934745 | 0.975249→0.954368 |
| pref_ref | 14→16 | 0.960261 | 0.986487→0.980387 |
| pref_ref | 16→20 | 0.945799 | 0.980387→0.961500 |

The size estimator is `3 ΔS_half / log(L2/L1)` for periodic chains.
The earlier same-code L16→20 controls give approximately 0.993 (XX),
0.499 (Ising), and 1.940 (two independent XX chains, a range-two control).
These preferred attractive rings do not independently support an
above-one entropy coefficient at these lengths. They also do not rule
out longer crossover scales or a dilute additional sector. The open-chain
entropy excess must survive larger lengths and order diagnostics.

Raw Lanczos spacings remain explicitly distinguished from gaps above a
certified complete ground manifold; completeness of the low spectrum is
not established. No gap or Casimir-energy veto is applied.

The exact field is zero. The rational point has coefficients
`(Jx,Jy,D,Jz)=(14,-5,11,-12)/24`, with raw spread 2.8.
The reference is `phi=1.2, Delta=.8, V=-2, t=1` in the stated
fermion convention. It has raw spread approximately 2.656.

Reproduce with `prepare_rings.py`, the registered `ed_followup.py` command
in the continuation protocol, and `analyze_rings.py`. Per-sector residuals,
raw/resolved entropies, coupling metrics and source hashes are retained
in `preferred_rings.jsonl` and `preferred_ring_summary.json`.

## R=8.75 bridge

The same registered calculation for `(Jx,Jy,Jz,D,hz)=(1/2,-1/2,5/8,1/10,-7/8)`
also selects even parity and momentum zero at every completed length.

| Sizes | Half-chain size slope | Profile slopes |
|---:|---:|---:|
| 14→16 | 1.049799 | 1.033126→1.030991 |
| 16→20 | 1.039646 | 1.030991→1.023765 |

These small rings are close to the one-mode control. The larger open-chain
entropy excess therefore has no independent periodic-ring confirmation
yet. A longer crossover length is possible; a c>1 CFT is not established.
