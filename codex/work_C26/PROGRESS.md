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
