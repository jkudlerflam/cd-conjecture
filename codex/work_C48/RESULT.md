# C48 result

[proved, restricted massive-matching theorem] The resummed cluster route establishes exponentially decaying, all-degree retarded kernels for the enlarged rooted cell model specified in [PROOF.md](PROOF.md), uniformly in cell count and at all sufficiently low temperatures. After C45's local transformation, the correction to the complete retained second-order Hamiltonian is linear in the remainder, rather than in the much larger interaction inside the massive sector.

[assumed, full input S remains open] This does not prove the requested statement for the complete microscopic fixed-length chain. The physical ensemble and source transfer, illegal-sector treatment, filling response and gapless multichannel RG remain separate obligations. No central charge or unconditional counterexample is claimed.

## Quantitative result

[proved conditional on the imported C42/C44/C45 rooted inputs] For an output spatial exponent alpha>=0, set

    b=n+3,  d0=4464,  m=n+1,
    C0=d0+d0^b+d0^(b+1)+b d0^(b-1)+2d0^(b+2),
    L=8(alpha+log32+10),
    T_n=max(2L,log(2d0),log C0+2),
    a=L+(n+1)log128+1.

Use the complete cell gap Delta, including the small charge tunneling splitting. PROOF (20) gives its explicit C44 lower-bound prescription. Then the kernel decay rate and temperature range can be taken as

    kappa = Delta L/(16T_n),       beta>=2T_n/Delta.

The matching and two-moment constant is

    A_dec,n = 2^22 exp(L) [2+16(16T_n/L)^2].

At fixed n and alpha it is independent of Delta, wrapper ratios and logical bandwidth. Recompute C45's v0,f0,G,C_SW at support exponent a and choose lambda by PROOF (26),(29). The zero-temperature influence action satisfies

    ||delta S_R||_LE <= A_dec,n C_SW lambda^3,
    ||delta S_R(omega)-delta S_R(0)||_LE
        <= A_dec,n C_SW lambda^3 |omega|/Delta.

Since W=c_n lambda^2 with c_n=J0/s and C42's full visibility/propagation factor

    J0=2g^2 q^(2n+2) D^(-2n-4)/S_q^2,

the relative error is at most A_dec,n C_SW lambda/c_n. A fixed positive lambda makes it arbitrarily small, independently of the ring size. Constants are deliberately very loose and are not uniform in n.

[proved, finite-temperature term] The thermal action also contains a baseline bounded, before the field-weight/moment conversion, by

    B_* exp[-beta Delta/(4T_n)],
    B_*=2^16 exp(L) Delta/T_n.

It must be retained. The infinite-cell rooted trace need not exist for all beta>0, so the theorem's temperature range cannot be replaced by an unrestricted all-beta assertion from these hypotheses. The range includes beta tending to infinity at fixed microscopic couplings.

## What was resummed

[proved] The cell's whole excited Hamiltonian is retained as an auxiliary operator. Only its finite logical matrix blocks are expanded. Their coefficient norm costs at most (8h^2)^(mr) on r cells, independently of the number of excited states. Coarse time blocks sum all repeated interval insertions before the connected expansion. A weighted activity sum treats unbounded interval supports, and a Schatten Holder bound controls thermal windings. The proof includes the hard-core tree bound and rooted-tree recursion, without invoking an unread convergence theorem.

[proved, C46 is addressed] In C46's model q_aux(Delta+u n1 n2), the raw order-lambda^2 external-degree sum still diverges for u>=Delta. The resummed representation instead uses four auxiliary eigenstates of energies Delta,Delta,Delta,Delta+u. Its order-two kernels have finitely many endpoint monomials and factors exp(-E_s t), with no expansion in u/Delta. Their operator integral is

    -lambda^2/Delta
      +lambda^2 u/[Delta(Delta+u)] n1 n2.

This changes the extension of retained variables off the code and requires the corresponding source matching. It does not assert convergence of C46's original raw action, or use a false Grassmann identity N(t)^2=N(t).

## Original model and remaining obligations

[proved] The embedding's invalid states cost at least Delta-w per invalid cell, where w is the full interaction incidence norm. Below that threshold, the enlarged and valid rooted spectra coincide. Their pressure difference is at most

    beta^(-1) log[1+C0 exp(-beta(Delta-w))].

This removes the extra states at the low-energy spectral and pressure levels. It is not an exact identity for arbitrary original thermal source functionals.

[assumed, remaining] To establish full S one still needs the complete physical bounded/local representation and sector comparison, transfer through the fixed-length constraint and cyclic quotient, control of the unbounded length source and filling derivatives, and constructive RG for the resulting arbitrary-degree retarded action. The cluster proof integrates no retained propagator and makes no claim about that last step.

## Verification and provenance

[numerical] Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py`. All checks pass. The resummed C46 kernel integral agrees with its operator resolvent to 8.7e-15. All 729 one-site periodic active/passive/vacuum histories satisfy polymer factorization to 6.8e-21, and their partition sum agrees with direct noncommuting transfer matrices to 6.7e-16. Exact spectral embedding and pressure comparisons pass on dimensions 8 and 20. These are small controls of the proof, not simulations of the microscopic wrapper or a CFT measurement.

[proved, attribution and scope] The abstract representation and coarse-time route were already developed in ../../contrib/direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md and its independent audit. C48 supplies a self-contained conservative derivation, explicit C42/C44 substitution, dimensionless uniformity and the application boundary. Only files in work_C48 were written. No network, subagents or cluster computation were used. The project STATE was outside the authorized write scope; [STATE.md](STATE.md) is the local handoff.
