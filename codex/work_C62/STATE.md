# C62 handoff

2026-10-28. Partial result. Lemma TP and the fixed-interaction XXZ thermal theorem remain unproved.

[proved partial] PROOF Section 3 defines a stable g2 Hamiltonian Luttinger reference with smooth compact transfer-momentum cutoff and matched v,K. Finite fermion-window compressions converge in the Gibbs trace at each fixed L,beta by form-core approximation and min-max. Their Duhamel/Grassmann normalization includes the free trace, so the total reference scalar equals the oscillator trace after removing regulators in the specified order. This does not give a uniform connected RG expansion or the common TP norm for that reference.

[proved] Hamiltonian parameters are A=v(K+1/K)/2, B=v(1/K-K)/2. In BM's tau convention, A=v(1+tau^2)/(1-tau^2), B=2v tau/(1-tau^2), g_H=2pi B=lambda_infinity/(1-tau^2). The BM Euclidean reference is a different regulated functional. Do not equate raw couplings or partition functions.

[proved] Finite-window current commutator equals lower-edge occupations minus upper-edge occupations. Its difference from m I has norm 2m at every cutoff. The infinite-sea anomaly is recovered on a finite-excitation core, not in norm. This invalidates direct finite-cutoff oscillator substitution, not TP. The trace convergence proof avoids that substitution.

[proved conditional] A joint mark scale sum with one shared short-memory factor converges with q_h^theta(T/q_h)^p for theta<min(alpha,sigma). p=3 and theta<1 are compatible. Obtaining the marked tree majorant itself is unproved. An exact analytic-flow example shows that affine interpolation between two critical bare inputs need not stay critical, so a global Cauchy proof needs its own domain control.

[remaining] PROOF Section 7 is the stopping point: prove the uniform full-kernel defect, four-history mixed comparison and complete scalar derivative/Hessian estimates in one normalized GN construction. Use the trace-derived reference above, or separately compare the BM Euclidean scalar to it with o(T^2) error. The ordinary XXZ RG construction is not being declared absent; its joint thermal vacuum extension is what was not established.

[checks] checks.py passes ten exact finite-CAR cases, five rational parameter cases, fifteen scale sums, an exact interpolation recurrence and two interacting reference trace controls. Largest matrix 256. These are not XXZ RG data. No failed checks in the final run.

[files] RESULT.md gives the verdict, PROOF.md all new arguments, SOURCE_AUDIT.md the primary-source boundary. All changes are inside work_C62. Parent STATE was left unchanged under the explicit folder restriction. No subagents, shell network commands, large computations, author messages or independent review.
