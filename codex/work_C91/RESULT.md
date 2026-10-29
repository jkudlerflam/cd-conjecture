# C91 result

[Assumed, unresolved main target] A verified finite-J mobile representation satisfying C48 has not been obtained. The dressed material gap, full-space grading or a replacement representation, and strong source-uniform action bounds remain open. The proposed logical action is stated conditionally in [PROOF.md](PROOF.md), section 7. This is not a completed initialization or a counterexample to the conjecture.

[Proved] The legal compression has the exact gap representation

    t H0+J sum_j |g_j=0><g_j=0| h_j,
    H0=-sum_j(E_j^+ E_(j+1)^-+h.c.),   ||E^+/-||=1.

It retains the fixed total gap and physical cyclic gluing. At finite K the legal space is not invariant, so this compression cannot replace the literal chain.

[Proved, quantitative ultraviolet progress] The collaborator's single-exponential construction admits explicit volume-independent constants for this bounded legal Hamiltonian. With M0=8192, v=64e, C_R=2v+9, R=ceil((2v+8)delta^-2), and q=20M0 R(J/t)/delta<=1/2,

    ||R_UV||_first-moment
       <=2^20 M0^2 C_R^4 t (J/t)^2 delta^-10
          +2^15 M0 C_R^2 J delta^-5 exp[-1/(4delta^2)].

An exponential rate (log 2)/(12R) is also available. The derivation includes the incidence count, filter boundary error and BCH remainder. It establishes ultraviolet locality, not a massive influence action.

[Proved, application failures] The following prevent a direct C48 substitution.

1. The filtered contact differs from its sea average by at least 0.545 times the material-term norm on some retained half-filled states for delta<=1/2. It filters energy transfer, not absolute bath energy. The frozen clock also has order-J code leakage already at lambda=0.
2. The available exponential rate tends to zero with the cutoff. It does not imply C48's fixed strong coefficient norm. Section 4 gives the stronger Hamiltonian majorant and its additional smallness condition explicitly.
3. Every diagonal product sign symmetry of the thirteen microscopic letters is blind to the stored bits. Such a grading cannot implement payload fermion parity. This does not exclude a dressed or block representation.
4. The finite-K color correction gives unequal positive payload fields with coefficients t^2 w_r/(JK), independent of lambda. Relative to the intended J alpha_c lambda^2 E_unit scale they grow as lambda decreases. Keeping this reference requires a positive lower bound on lambda at fixed J,K, as well as the matching upper bound.

[Proved, restricted parameter check] J/t=2^-5000, K=2^14000 and lambda=2^-608 satisfy the explicit ultraviolet first-moment and leading color budgets with a large margin. They do not certify the other initialization requirements. The failure to finish the proof is not a contradiction among all finite parameter choices.

[Numerical and exact finite controls] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. The 35-dimensional phase-model filter identity has residual 3.95e-15 or smaller. The length-fifteen single-contact clock leakage is 0.1507149408. All 8192 microscopic sign assignments were checked; sixteen conserve the literal moves and all are bit-blind. The color-field weights range from 0.0033428416 to 0.6864370767. These checks establish no mobile spectral gap or infrared limit.

[Proved, deliverables and limits] [PROOF.md](PROOF.md) contains the derivations, conditional action and missing input contract. [PROGRESS.md](PROGRESS.md) records the work; [checks.py](checks.py) and [checks.json](checks.json) reproduce the controls. [STATE.md](STATE.md) is the local handoff. Only work_C91 was modified. No network commands, parent-script execution, delegated agents or large computations were used.
