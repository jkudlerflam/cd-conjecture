# Mobile toy short-ring check

[numerical] Full d=6 bond from contrib/mobile_cells/ed.py, with t=1, mu=-1, J=0.5, Jz=0, K=10, and varying marker hop u. All legal wall-number and spin sectors were checked against unrestricted ED at L=3,4,6. They agree. Gauge optimization is bracketed, not asserted solved. The lower bound uses gauge-invariant expectations in equal-marginal states, the upper bound minimizes over diagonal coboundaries.

[numerical, limitation] These short rings are not certified vacuum-compatible along the incommensurate charge fluid. Therefore Q_3 is recorded as a raw finite-ring diagnostic and no A estimate is assigned. Increasing marker mobility adds at most the one length mode in this toy. It has no parameter encoding an unbounded payload. Larger sector-resolved work is needed to infer its microscopic A/M.

| u | M interval | raw Q_3 interval |
| --: | :-- | :-- |
| 0.1 | [5, 5.5] | [0.153304099, 0.168634509] |
| 0.3 | [5, 5.5] | [0.226601938, 0.249262131] |
| 1.0 | [5, 5.75] | [0.336866524, 0.387396503] |

[proved] For any vector with identical one-site marginals, the expectation of o tensor I minus I tensor o vanishes. Two such vectors therefore certify M >= half their expectation difference. This lower certificate applies to every Hermitian gauge, even though the numerical upper search only used diagonal gauges.
