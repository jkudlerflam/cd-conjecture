# C80 handoff

Completed 2026-10-28. Twelve-letter bond with a certified n=12 member.
Read RESULT.md, then PROOF.md for the explicit bond and all derivations.

The tokens carry no data bit. A head at position j displaces logical bit
j to the next site; the terminal empty marker carries the last bit.
This gives fourteen letters. Merge vacuum marker M_v with empty E_0,
then vacuum token D_v with R_0, giving twelve. The legal-word map is
injective and exactly reproduces C75's constrained parent. Conversion
is now D_A R_1 <-> D_B R_0, or D_A E_1 <-> D_B E_0 at the boundary.
The express fold is unchanged. Exact U(1) and marker conservation hold.

Final schedule n=12, D=1000, Delta=4-2cos(pi/14), phi_A=pi/24,
phi_B=0, psi_A=pi, psi_B=pi/2, beta=1/d_13 with recurrence
 d_(j+1)=1004d_j-d_(j-1), d_0=1,d_1=1004.
 t_d=beta 2^-66, g=2^-216, K=2^2048, lambda=2^-480.
 Physical field=lambda^2 alpha^4 beta^2 g^2 2^-112.
 mu_sel=2-(delta_13+delta_12)/2 as defined in PROOF.md.

Certificate: 24 simple roots, delta2>0.01509618, delta4>5.5894e-7,
300 unordered pair sums with replacement, all centers distinct modulo pi.
Allowed additional dyadic Bloch norm is 2^-35. Finite-g error<4.93e-18,
finite-K error including zero-port splitting<1.16e-55. The root routine's
finite_g_power=-205 is an inherited threshold; the final schedule is -216
from finite_budget.py. No numerical output was copied from C75.

Do not claim an exactly degenerate finite-K code. Its splitting is O(1/K)
and enters the physical error budget divided by lambda^2 E_unit. The
second-order ED uses the symmetrized denominators E_q-E_a and E_q-E_b.
ED closures: 346 and 8224 states at n=1,2, with fourth-order lambda scaling.
Exact grammar audit reaches n=6. Uniform zero-port selector proof follows
the fragment and partial-isometry argument, allowing the shifted bits.

Scripts: bond_shift.py, checks.py, run_certificate.py, finite_budget.py,
selector_checks.py, grammar_audit.py. alternatives.py, merge_c75.py and
fast_carrier.py preserve the failed smaller-design tests. Scalar certificate
and parent-reference source is copied locally from C75, with attribution.

Conditional infrared prediction c=12 at a=12, or c=13 at a=13 with one
vacancy. Open: nonzero-port dressed-defect comparison, volume-uniform
infrared transfer, and vacancy matching. Twelve is the smallest certified
alphabet here, not a universal minimum. No independent referee.
Only this folder was modified. No network or cluster activity.
