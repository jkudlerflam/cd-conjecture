# C83 handoff

Completed 2026-10-29. Read RESULT.md, then PROOF.md.
All files and changes are local to work_C83.

[Proved] Two colors V_0,V_1 are necessary and sufficient within the
nearest-neighbor diagonal repair preserving every legal vacancy placement.
Dimension 14. C80 outgoing states are 0 for L and all markers, 1 for
R,D,X. Incoming sets are {0,1} for R_0,E_0, {1} for R_1,E_1 and {0}
otherwise. Vacancy V_q accepts and outputs q. Do not copy C82's incoming
sets unchanged. The literal independent bond is bond14.py.

[Proved, computer-assisted] n=13, b=14, 26 simple roots, 351 pair sums,
delta2>0.0149087301, delta4>8.70283e-8, distinct centers modulo pi.
D=1000, alpha=sqrt(2/15)sin(pi/15), Delta=4-2cos(pi/15), phi_A=pi/26,
phi_B=0, psi_A=pi, psi_B=pi/2. beta=1/d_14 from d_(j+1)=1004d_j-d_(j-1).
t_d=beta 2^-74, g=2^-256, lambda=2^-576, K>=2^2304,
field=lambda^2 alpha^4 beta^2 g^2 2^-121. Midpoint mu_sel=1.8742920463679156.
Finite errors: g<1.87e-33, K including zero-port splitting<1.12e-40,
local two-cell port remainder<2.60e-50, allowance 2^-35.
Raw root finite_g_power=-219 is an intermediate threshold only.

[Proved] F3 frozen selection transfers after assigning a left-marker
boundary penalty to the right cell, since all marker outgoing states
are identical. gamma>2^-16 and gamma>5sqrt(2)lambda. F4 finite-order
local normal form transfers with penalty JK in the lift. Neither proves
finite-J mobile matching or a mobile defect comparison.

[Numerical] Literal invariant-closure leakage tests dimensions at most
4608 recover the uncolored plateau .0256 and colored K^-2 decay. Full
nonzero-port Q=2 tests with one and two vacancies have dimensions 268,
1045. Largest residual <1.55e-11. Legal projection identity and all
38400 nonempty four-site grammar words pass. Complete frozen sector
comparison dimensions 1024 passes. SW check dimension 60 passes.

[Assumed] A positive infrared tolerance, strong-norm massive matching,
finite-J vacancy matching, mobile dressed-defect stability, half-density
tuning, retarded multichannel RG and conformal identification remain.
After choosing positive J, increase K to control nu^2/(JK) relative to
J times the logical tolerance. K=2^2304 is a frozen lower bound, not a
mobile matching window. exp(i pi Q/2) still has +/-i auxiliary eigenvalues,
so the old microscopic Z2 grading cannot be imported unchanged.

[Heuristic, conditional] d=14,c=14. No unconditional counterexample,
central-charge measurement or universal d=13 no-go claim.

Reproduce with certify_n13.py, budget_n13.py, check_lift.py,
check_frozen.py, check_sw.py, all with python3 -B in this folder.
C80 scalar source copies preserve C75/C33 attribution. No network,
cluster jobs, subagents or writes outside this folder. The project STATE
was intentionally left untouched under the explicit folder restriction.
