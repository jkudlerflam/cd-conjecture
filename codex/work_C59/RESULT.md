# C59 result

[proved partial results; target unproved] Estimate P has not been established here for interacting XXZ at fixed nonzero anisotropy. The requested adaptation of the Ising proof remains incomplete. [PROOF.md, Section 7](PROOF.md#7-exact-remaining-lemma) states the remaining lemma TP, including the full effective-kernel norm, actual thermal histories, vacuum scalars and reference normalization. No unconditional thermal c=1 theorem or extension to the matched retarded action is claimed.

## Literature result

[read, announcement only] The [December 17, 2024 Rutgers program, page 10](https://math.rutgers.edu/images/127_SMC_short_talks.pdf#page=10) confirms the Renzi-Giuliani-Mastropietro-Santini announcement. Searches through 2024-2026 arXiv and university results did not locate its manuscript. The abstract gives no theorem hypotheses, so its coverage of off-half-filled XXZ and of the matched retarded action remains unverified. A located 2026 Fabbri-Mastropietro-Renzi paper, [2602.03705](https://arxiv.org/abs/2602.03705), concerns massive lattice QED2, not the announced thermal theorem. The negative search is limited by failed author-page retrievals; details are in [SOURCE_AUDIT.md](SOURCE_AUDIT.md).

## Mathematical result

[proved distinction] The Ising scalar estimate cannot be applied to the XXZ interaction scalar itself. With a surviving marginal coupling, even a generated six-field reference kernel can have a nonvanishing dimensionless size O(g^2). The full interacting reference must be subtracted before the infrared gain is extracted. Convergence of the localized quartic coupling does not establish convergence of the complete effective potential to its localized part.

[proved partial lemmas] PROOF Sections 3-6 establish the following ingredients.

- A lattice-curvature covariance replacement gains one shell power, including any fixed weighted moment.
- Matched anomalous field normalizations differ by a finite amplitude times 1+O(q_h^theta), under a summable difference of their log increments.
- Smooth Gaussian shell log determinants and their first two parameter derivatives obey a thermal bound C q_h^2(T/q_h)^m for m>=3. Their four-history comparison therefore has both required gains once the state histories do.
- An explicit weighted-history criterion permits relevant amplification by two and a nonzero small marginal coupling. Its contraction factor is given in PROOF (16).
- A trace-normalized oscillator reference with momentum-dependent speed tending to v has thermal coefficient -pi/(6v). Identifying the RG reference scalar with that trace is still required.

[assumed, remaining estimate] The required all-depth comparison must establish

    ||X_h^0-Y_h^0|| <= A q_h^theta,
    ||(X_h^T-Y_h^T)-(X_h^0-Y_h^0)||
                              <= C q_h^theta(T/q_h)^p,
    ||D E_(h,T)||+||D^2 E_(h,T)|| <= C q_h^2,
    ||D(E_(h,T)-E_(h,0))|| <= C q_h^2(T/q_h)^p,
    p>2+theta,

in compatible normalized kernel/history norms, with the exact scalar expansion and terminal treatment. Lemma TP gives local/history hypotheses sufficient for these conclusions. The outstanding operation is the uniform renormalized sum retaining both a reference defect and a thermal replacement. The collaborator already proves the fixed thermal-mark allocation conditional on ordinary marked-tree routing; that result is imported, not duplicated. Its full XXZ application and the trace-normalized reference identification were not completed.

[proved conditional theorem] Under TP, with fixed physical field and the thermodynamic limit first,

    f(T)=e0-(pi/(6v))T^2+o(T^2),

where v is the dressed sound speed. The lattice/reference thermal difference is O(T^(2+theta)); an independently controlled o(T^2) reference-normalization error is also allowed. Fixed-density conclusions need estimates uniform near the selected field and the corresponding ensemble argument. No temperature-dependent physical retuning is used.

[assumed retarded extension] An instantaneous XXZ theorem would not alone establish TP for the all-degree matched retarded action. That application still needs the nonlinear endpoint estimates, rematching of the static marginal and kinetic data, and source-uniform scalar control for C58. With broken reflection the coefficient is -pi(1/v_R+1/v_L)/12. C58's physical transfer is already available conditionally and was not redone.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. The largest Poisson-identity discrepancy, including two parameter derivatives, is 3.11e-17; the test history contraction factor is 0.310149. These are checks of the partial lemmas, not an interacting XXZ calculation. C57's first-order coefficient was used without recomputation. All writes stay in work_C59; [STATE.md](STATE.md) is the handoff.
