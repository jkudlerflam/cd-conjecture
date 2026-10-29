# C100 state

Completed 2026-10-29. This folder studies thermal and transport routes
to A = sum_a c_a v_a <= pi M log_2 d for nearest-neighbor chains.
The conjecture remains open. No microscopic counterexample was found
or claimed.

Read RESULT.md for the conclusions and PROOF.md for definitions,
derivations, gauge caveats, and transport-limit assumptions.

Settled within this investigation:

- With ||g|| = M in the optimal centered bond gauge, j_i = i[g_(i-1),g_i]
  has norm at most 2 M^2.
- The positive spectral moment m1(q,T) = <[h_-q,[H,h_q]]>/(2L)
  obeys 0 <= m1 <= 8 M^3 q^2. This is an unconditional lattice bound.
- Its CFT value is pi A T^2 q^2/6 + q^4 sum c_a v_a^3/(24pi).
- A positive low-energy spectral matching factor alpha at (q_*,T_*)
  implies A + q_*^2 sum c_a v_a^3/(4pi^2 T_*^2)
  <= 48 M^3/(pi alpha T_*^2).
- The CFT total-current variance is pi A T^3/3 per length. In the
  convention Re kappa = pi D_th delta + regular, D_th = pi A T/3.
  This is not an automatic statement about a lattice dc Drude peak.
- The CFT static density susceptibility is pi B T^2/3 + A q^2/(12pi),
  B = sum c_a/v_a. Its inverse-frequency moment lacks a microscopic
  bound from the norm and positive moments used here.
- Entropy matching at T_s gives A <= 3 v_max^2 log d/(pi alpha_s T_s).
  Combining with C69's conditional v_max <= 6 e M gives
  A <= 108 e^2 M^2 log d/(pi alpha_s T_s).
- The ring ground energy differs from L e0 by at most 2M. Casimir
  matching at L_* gives A <= 12 M L_*/(pi alpha_E).

All continuum-to-lattice comparisons retain a matching scale. The
analytic thermodynamic family in PROOF.md disproves an inference from
scalar thermodynamic capacity and convexity constraints alone; it is
not a spin-chain construction. The investigation does not exclude
other thermal proofs. A new bound on the static susceptibility's
q^2 coefficient is a concrete open possibility.

checks.py uses NumPy and SciPy for small normalization checks; final
results in CHECKS.md all pass. No large numerical work or external
files were modified. This local memory substitutes for a parent-state
edit, which the task's write restriction forbids.
