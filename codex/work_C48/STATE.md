# C48 local state

[proved, completed] RESULT.md and PROOF.md contain a resummed massive-matching theorem for an explicitly specified enlarged rooted model, not a proof of full microscopic S. It uses the prior operator-valued route in contrib/direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md with attribution, supplies the hard-core cluster proof, and inserts the C44/C46 heat bound.

[proved, new constants] With alpha the output spatial exponent, L=8(alpha+log32+10), C0=d0+d0^b+d0^(b+1)+b*d0^(b-1)+2*d0^(b+2), T_n=max(2L,log(2d0),log C0+2), tau=T_n/Delta and kappa=Delta*L/(16T_n). The matching constant A_dec,n=2^22*exp(L)*(2+16*(16T_n/L)^2) is independent of wrapper ratios and Delta at fixed n. Input support exponent a=L+(n+1)log128+1 must be chosen before recomputing C45's local SW constants. The explicit lambda schedule is PROOF (26),(29).

[proved, C46 mechanism] Treat the entire cell excited Hamiltonian spectrally before external coefficient expansion. The toy example then has four auxiliary levels, finitely many endpoint monomials and no expansion in u/Delta. The original raw action remains divergent. Retained source extensions off the code change and must be tracked.

[proved, limitations] Uniform beta means beta>=2T_n/Delta. The rooted heat trace may not exist at arbitrary high temperature. Invalid states of the embedding have a threshold Delta-w and exponentially small pressure contribution. This does not prove exact thermal source transfer to the physical fixed-length ensemble or dressed exclusion of physical illegal cells.

[assumed, open] Complete local fermionization and bounded interaction representation for all required physical sectors; illegal-sector comparison; fixed total length and cyclic quotient; unbounded length sources and filling derivatives; multichannel retarded RG and CFT/current/stress limit. No c=n+1 theorem is claimed.

[numerical, checks] checks.py passes. It writes checks.json, using only numpy/scipy/standard library. Largest matrix dimension is 20. Exact-embedding, pressure, noncommuting semigroup, 729 periodic polymer histories and dimensionless constant checks pass. No network, subagents or files outside work_C48 were used. Project-level STATE was not modified because the task forbids outside writes.
