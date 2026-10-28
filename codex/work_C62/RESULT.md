# C62 result

[unproved target] Lemma TP remains unproved for interacting XXZ. No rigorous thermal c=1 conclusion at fixed nonzero anisotropy follows from this run. The failure is a missing uniform joint comparison, not a counterexample to the thermal formula.

[proved partial result] [PROOF.md Section 3](PROOF.md#3-a-specific-hamiltonian-reference-and-its-fermion-regulator) constructs a trace-defined spatially smoothed Luttinger reference with the prescribed dressed v,K. Its finite fermion-window Hamiltonians have Gibbs traces increasing to the full Hamiltonian trace at fixed L,beta. The proof uses form-core approximation and spectral min-max. Their exact Duhamel/Grassmann representation retains the free Gaussian trace and every scalar. Thus this regulator gives the oscillator scalar, whose thermal coefficient C51/C59 already evaluates as -pi/(6v).

[proved boundary] This identifies the total scalar for the chosen Hamiltonian-derived reference. It does not prove that its regulated effective kernels admit the uniform BM comparison in TP, or identify it with BM's different spacetime-smoothed reference. A cutoff-independent connected expansion and its common kernel norms remain unverified. The trace argument itself requires no uniform-in-T convergence because the regulator is removed at fixed L,T before the thermodynamic limit.

[proved parameter check] The Hamiltonian stiffness parameters corresponding to v,K are

    A=(v/2)(K+K^(-1)),     B=(v/2)(K^(-1)-K).

With BM's tau=lambda_infinity/(4pi v) and K=(1-tau)/(1+tau), these become

    A=v(1+tau^2)/(1-tau^2),
    B=2v tau/(1-tau^2),
    g_H=2pi B=lambda_infinity/(1-tau^2).

Identifying the raw reference couplings without this conversion gives a wrong sound speed already at second order. Matching these parameters is necessary but does not prove a full effective-action comparison.

[proved finite-cutoff obstruction] Truncated fermion currents do not obey an exact central current algebra. Their commutator is the difference of occupation operators at the two edges of the sea cutoff. Its deviation from m I has norm 2m for every cutoff. Consequently finite-cutoff oscillator substitution is invalid. The trace convergence argument above supplies a legitimate alternative for this reference.

[proved conditional scale estimate] A joint marked-tree majorant with one short-memory factor, not two independent factors, suffices to retain q_h^theta(T/q_h)^p. PROOF Section 5 gives the explicit convergent double scale sum. Thus p=3 and a sufficiently small positive theta are compatible. The obstruction is obtaining that majorant in the actual normalized expansion, not summing it afterwards.

[unproved, exact remaining estimate] [PROOF Section 7](PROOF.md#7-the-remaining-obstruction-for-xxz) requires the full lattice/reference defect and mixed history bounds

    ||X_h^0-Y_h^0|| <= C q_h^theta,
    ||(X_h^T-Y_h^T)-(X_h^0-Y_h^0)||
                             <= C q_h^theta(T/q_h)^p,

together with the compatible vacuum derivative/Hessian bounds, Gaussian factors and ultraviolet/terminal terms, uniformly in depth and regulators. The collaborator's single-shell winding, localization and conditional moment allocation do not establish this joint all-depth estimate. A global Cauchy argument also needs a matched critical interpolation; PROOF Section 6 gives an exact analytic-flow counterexample to inferring that from its two endpoint trajectories.

[numerical/exact checks] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. It tests ten finite CAR commutators, five rational parameter matches, fifteen joint scale sums, an exact critical-manifold interpolation, and finite interacting reference traces against the oscillator trace. The largest matrix has dimension 256. These are controls of the partial results, not a numerical verification of TP.

[scope] All files are in work_C62. No retarded-action extension, interacting XXZ simulation, lattice CFT assumption or author correspondence was used. The exact next task is the joint normalized GN estimate for the trace-derived reference, or that estimate for the Euclidean reference plus an o(T^2) scalar comparison between the references.
