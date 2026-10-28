# C26 progress

2026-10-28. G5 and G6 only. All writes stay in work_C26. No network or agents.

[proved, task accounting] Read project STATE, route-A background, assembly README,
repository STATUS, C19 RESULT and cell/check sources, ir_control Sidon construction,
and research_continuation/IR_AUDIT. Read C21, C22 and C8 progress. C25, C27,
routeA_compile and routeA_spatial progress were absent; routeA_demonstrator has
source files but no progress. No local prior progress existed.

[assumed, inherited] The spatial compiler, mobile wrapper composition and
thermodynamic Gaussian identification remain separate tasks. This run treats
the complete abstract C19 reflection cell with exchange-generated fields.

[proved, derivation being recorded] With PVP=0 and R=Q H0^{-1} Q, the canonical
Hermitian fourth-order operator is
-PVRVRVRVP + {PVR^2VP,PVRVP}/2. Omitting its second term would not compute
the canonical low-energy Hamiltonian. The two-cell computation will report
both contributions, the traceless norm and the non-Gaussian residual.

Next: replace the geometric target by the Sidon points with the C19 sign
convention theta=-(p_L+p_R)/2, rerun exact filter and matrix identities,
and compute K3 and K4 at n=1,2 with explicit visibility-rescaling checks.

## G5 and two-cell G6 checks completed

[proved] The target replacement changes only t, theta and mu in the C19
weights. Positivity, reflection completion, parity and the second-order
contraction retain their proofs. The code uses theta negative of the
Fermi-sea center, as C19's dispersion requires.

[numerical, n=1,2] Full cells have dimensions 96 and 256. Source-measure,
code and direct second-order identities pass. Exact rational Jacobi tests
pass at n=1,2,4,8,12. The canonical fourth-order traceless norm divided by
rho^2 is 1.23760530 and 4.37656212. Non-Gaussian residuals divided by rho^2
are 0.00066164530 and 0.00051791245. Independent full and reducing-subspace
calculations agree below 1.7e-14. Finite-coupling canonical graph solutions
at lambda/Delta=0.08,0.04,0.02 show fifth-order remainders after K4.

[proved, derivation being recorded] Analyticity alone plus parity is not
the full rho argument. Bare logical parity conjugates H(rho) to H(-rho),
and the physical conserved parity induces that same logical action.
Together they force the canonical effective Hamiltonian to be even in
rho. At rho=0 the logical register is a spectator, so its traceless part
is divisible by rho^2 at every fixed perturbative order.

Next: include the connected adjacent-bond term at fourth order, then state
the second-order local-integration and infrared hypotheses explicitly.

## Adjacent-bond calculation and independent check completed

[numerical, n=1,2] The complete connected fourth-order calculation includes
both a repeated bond and two adjacent bonds each used twice. The latter
was computed in invariant spaces of dimensions 55296 and 1048576, with
code dimensions 8 and 64. Its traceless norm/rho^2 is 0.8622186591 and
1.6116197398. Its non-Gaussian residual/rho^2 is 0.0007054490 and
0.0003889720. Visibility multipliers 1,1/2,1/4 keep the normalized
coefficients finite. Each n=2 run took about 15 seconds.

[numerical, n=1] An independent direct tensor calculation of the complete
three-cell V agrees with the streamed K4 within 2.1e-14. Its nonlinear
canonical graph solutions at lambda/Delta=0.04 and 0.02 have residuals
below 3.9e-16 in the exact-integer energy representation. Fourth-order
truncation errors decrease by a factor 33.22 under halving lambda.

[proved] The canonical fourth-order chain has no further connected support:
a three-bond cluster needs at least five insertions to avoid a singly
acted-on outer cell. Disconnected terms cancel by factorization of the
canonical polar code map. Thus the computed bond and adjacent-bond terms
cover fourth order for the fixed-cell chain.

## Infrared reformulation and handoff

[proved, finite-order] RESULT.md gives the full derivation, including the
positive weights, reflection, parity, second-order contraction, canonical
K4 subtraction, adjacent-bond formula and the two-symmetry rho argument.
The Sidon phase convention is theta minus the sea center. The jitter is
clipped against the scaffold margins for arbitrary n; this does not change
the n=1,2 examples. Exact rational filter tests and source identities pass.

[assumed/heuristic, explicitly separated in RESULT] IR_AUDIT equations
(5)-(10) have been restated for W=lambda^2 rho^2/Delta_src. The correction
starts at lambda^3 after subtracting the complete lambda^2 target. With
a symmetry-preserving local massive-sector bound it is relatively
O(lambda/Delta_src), plus synthesis, splitting and wrapper errors. Local
convergence, dressed Gaussian stability and mobile gluing remain inputs.
No thermodynamic central charge or all-n order-one coefficient was claimed.

[proved, implementation limitation] Payload-scalar terms may depend on
cell lengths and then are not global scalars. They can appear at
lambda^2/Delta_src, without rho^2. G3 must equalize, compensate or control
them. Fixed code splitting also cannot be repaired by only shrinking
lambda. These are explicit hypotheses, not a failure of the Sidon cell.

[proved, completion] G5 and the requested G6 calculations are complete.
The full route-A construction remains conditional. All writes stayed in
work_C26, so the outside project STATE was not changed. This file is the
local handoff memory. No network, cluster computation, external messages,
subagents, commit or push were used.

## Final validation

[proved, task accounting] Read the now available C25 and C27 progress.
C25's finite-vacancy source mismatch and scalar corrections are consistent
with the wrapper qualification above. C27 supplies conditional projected
gluing; its finite-ring work was not repeated. The three Claude progress
files remain absent.

[numerical] The final Sidon script was rerun after adding explicit
assertions for every algebraic diagnostic and clipping the jitter for
arbitrarily large n. All assertions passed; n=1,2 parameters are unchanged.
The three-cell and independent suites already passed and were not repeated
for this inactive-at-tested-sizes change. All four Python sources parse,
all result JSON files parse, and the report's formula delimiters and text
constraints pass. No numerical failure was suppressed. RESULT.md is the
complete derivation and scope statement; STATE.md is a compact local index.
