# C75 handoff

Updated 2026-10-28. The smallest alphabet found is 15, with n=15 and
cell length 16. Conditional infrared prediction c=15 at d=15, or c=16
at d=16 with one vacancy. No central charge measurement, no global
minimality theorem, and no completed microscopic counterexample.

Read RESULT.md for status, PROOF.md for the explicit bond and derivations,
ALPHABET.md for reductions and failed/open alternatives.

The construction first replaces the three excited species by two of
relative charges +/-1/2, giving 18 letters and exactly half of the C33
common-loop Gaussian matrices. Then identify M_v with M_empty and
D_(v,a) with R_a. Legal words remain distinguishable. The resulting
15-letter bond has a vacuum L/R domain wall. Extra off-grammar moves
are retained at finite K.

Final schedule: n=15, D=1000, Delta=4-2cos(pi/17), phi_A=pi/30,
phi_B=0, psi_A=pi, psi_B=pi/2, port phases -pi/4,+pi/4,
t_d=beta 2^-84 with beta=1/d_16 and d_(j+1)=1004d_j-d_(j-1),
g=2^-256, K=2^1152, lambda=2^-640,
field=lambda^2 alpha^4 beta^2 g^2 2^-140.
The marker selector is the C22 midpoint at b=16.

Exact certificates: 30 simple roots, delta2>0.0145540,
delta4>1.1444e-7, all 465 unordered pair sums with replacement.
finite_budget.py transports the leading certificate to finite g,K.
The copied C33 complex product-radius formula was conservatively
corrected from a constant 100 to 200. Final n=15 outputs include it.
The root_certificate finite_g_power=-247 is the inherited scalar
routine's threshold, NOT the final 15-letter choice. Use -256 from
finite_budget.json and RESULT.md.

Uniform zero-port selection is proved: for k data sites,
mu+E_v(k)-(k-1)_+/(K-4) <= E_K(k) <= mu+E_v(k).
Split the vacuum moves into bad-pair-preserving moves and partial
isometries creating one bad pair. In a fixed bad pattern there are
at most m+1 legal fragments, each with at most one head and kinetic
row sum two. Young's inequality absorbs the remaining moves into
(K-4)N_bad. Neighboring selector density difference >2^-16,
selected density <-2^-8, so b=16 survives finite K for all lengths.
The finite-cell code gap is protected by a separate Schur estimate.

Checks: bond18/17/15 builders; n=1,2 literal legal-space projection,
U(1), Gaussian g^2 coefficient, exact second-order source resolvent;
ED dimensions 436 and 15404, errors scale as lambda^4.
Unshifted Lanczos missed degenerate-band multiplicities. Its failed
output is retained. Use shift-invert with 80 Krylov vectors at n=2.

Remaining: independent audit of the new finite-K arguments, uniform
nonzero-port dressed-defect comparison and infrared matching. Vacancy
lift adds its own existing matching assumptions. No complete 14-letter
design or replacement for the express fold was found.

Reproduce with run_certificate.py, checks.py, selector_checks.py.
All files stay here; the project STATE was read but not modified.
