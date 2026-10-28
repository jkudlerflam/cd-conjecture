# C56 result

[proved, verdict] C50 (7.1) remains unproved for the matched retarded action, including at n=1. The full physical counterterm map is constructed only conditional on the collaborator's explicit RG hypotheses. No all-scale Luttinger theorem or unconditional counterexample is claimed.

[proved, counterterm continuity implication] The collaborator's COUNTERTERM_INDUCTION gives a unique bounded critical trajectory under G1-G2. Together with uniform bounds over a compact set of trial Fermi points and continuity at every fixed scale, this implies a continuous infinite-scale counterterm map. The proof uses a compact graph and uniqueness, without differentiating Fermi phases or a critical self-energy. Finite-cutoff maps converge uniformly. Trajectories converge in every weaker weighted norm. See [PROOF.md sections 5-6](PROOF.md).

For the exact initialization

    a_0^phys(p')=A(p')[e(p')-mu]+u(p'),

the map is

    nu(p')=A(p')^(-1)[u(p')-a_0(p')],
    ||nu||<=C_nu epsilon.

Brouwer then fixes the roots of the original physical action with displacement at most 2 C_nu epsilon/v_*. Independent microscopic chemical-potential adjustments are unnecessary. G1-G2 and the exact initialization remain application hypotheses, so this is not an unconditional verification of C50 section 8.

[proved conditional on the imported C48 polymer construction] Keeping C48's occupied connected network strengthens its output to an exponential terminal-tree norm. With its support parameter L_0 and block time tau, safe rates are

    a_tree=alpha/8,       kappa_tree=L_0/(128 tau).

The same cluster majorant controls the stronger norm and its leakage-linear difference. This specializes the collaborator's network method to C48. It does not redo massive matching. A separate example proves that exponential diameter weights alone cannot imply the polynomial terminal-tree norm at arbitrary degree, even with exact charge conservation and one species. See PROOF section 2.

[proved, first-order map] For the compact-patch theory, the all-degree one-endpoint self-energy sum is real, small and continuous in trial Fermi points under the stated uniform moment and smooth-projection assumptions. Its truncated root equation can be solved. No uniform bound on the omitted higher-order infrared terms is proved. See PROOF section 4.

[assumed, remaining decisive estimate] The missing input is the normalized renormalized-tree comparison, including a marked-history bound that retains the extra small factor after field normalization. PROOF section 3 displays the exact estimates that would yield (7.1). The collaborator's Gaussian multivelocity algebra is used as established reference input and is not mistaken for this comparison. At n=1 there are four chiral branches and six possible density couplings; arbitrary external degree and repeated localization remain. General n changes finite constants and label counts, not the conditional continuity argument.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. The checks cover connected-network geometry, the norm counterexample, a solvable relevant recurrence and a four-branch root equation with velocity ratio 3000. The recurrence formula error is below 8.7e-19, and the root error is below 6.8e-21. These are controls of the lemmas, not an RG calculation or a simulation of C39.

[proved, scope] Files are confined to work_C56. Primary sources were consulted with the browsing tool as permitted by S2. No shell network commands, delegated agents, cluster jobs or writes outside this folder were used. The parent STATE.md is outside the allowed write scope; [STATE.md](STATE.md) is the handoff. No independent review of these new arguments has been performed.
