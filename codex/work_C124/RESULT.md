# C124 result

2026-10-30.

[Open] The reconnection lemma remains unproved for general noncommuting stoquastic chains and for general interacting chains. No counterexample to that lemma or to A<=C M f(d) was found.

[Proved] The lemma holds for every chain with commuting adjacent bonds, including complex bonds with entangled eigenstates. An on-site basis adapted to the two commuting interaction algebras gives a nonnegative matrix B on recorded histories with

    q_(L,s)=B_ss,  q_(2L,s,t)=B_st B_ts.

Consequently

    sum_s q_(L,s)^2 <= Q_D,  G<=0

at every spacing, every N, and every L>=2. This supplies b_0=b_1=0 with that basis-selection rule. The proof applies to the positive recorded histories, including their reflected contours.

[Proved] In an arbitrary measurement basis the same commuting subclass satisfies

    R_L<=log r<=log d,
    G<=log r+2N log d,

where r is the number of central sectors of the on-site interaction algebra. Thus b_0=log d, b_1=2 log d work without a basis restriction. Under the specified conformal infrared assumption this subclass has A=0. It does not supply a new theorem for interacting critical models such as Heisenberg or SU(3).

[Proved obstruction] C123's separate sufficient boundary-overlap estimate cannot hold uniformly, even for diagonal chains. For

    k=diag(0,epsilon,1),
    h=(k tensor I+I tensor k)/2,  M=1/2,
    z(t)=1+exp(-epsilon t)+exp(-t),

its temporal transfer has lambda_S=lambda_D=1 and w_D=1. The boundary factor is

    C_L=log(||u_S||_2^2/w_D)=log z(2tau)+(2L-2)log z(tau).

Set epsilon=L^(-2) and tau=2aL. Then C_L=2L log 2+O_a(1), although the actual reconnection cost is

    G=log z(2beta)-2log z(beta)<=0.

The transfer inequality itself remains correct. Separately bounding its boundary factor discards an extensive thermal contribution that cancels in the actual ratio.

[Proved] The obstruction persists in a critical stoquastic d=4 chain, XX tensored with a slow on-site qubit. Its boundary factor diverges at least as 2L log 2-O_a(1), while its G is at most the XX value. Thus excluding chains with an empty conformal sector does not repair that sufficient condition.

[Open remaining estimate] Keep the full thermal contractions S_n,D_n at fixed tau and set

    Gamma_n=log(S_(n+1)/S_n)-log(D_(n+1)/D_n).

The exact identity is G_N=G_1+sum_(n=1)^(N-1) Gamma_n. A sufficient completion is a uniform bound on G_1 for 0<tau M/L<=a and on Gamma_n for a/2<=tau M/L<=a. Equation (16) of PROOF.md expresses Gamma_n as a ratio of Rayleigh quotients with evolved endpoint vectors. The exact target only requires controlling the sum, so the separate increment bounds are explicitly stronger sufficient conditions. Neither is proved here for noncommuting chains.

[Numerical] Independent full-matrix tests verify the commuting separator to 3.50e-14, including complex virtual bonds and multiple central sectors. Fresh XX calculations reach L=8192, doubled length 16384, and six insertions. Full-spectrum ED tests Heisenberg and XXZ at Delta=-0.9,0.5 through doubled length 12, and SU(3) on the vacuum-compatible lengths 3 and 6, through eight insertions. Direct temporal-transfer contractions agree with the ED history values of G to 1.87e-13. The extension G<=0 to noncommuting models fails in every requested family at sampled parameters. No numerical extrapolation is promoted to a theorem.

[Work record] Derivations are in [PROOF.md](PROOF.md), numerical scope and reproduction in [CHECKS.md](CHECKS.md), and the handoff in [STATE.md](STATE.md). All writes and calculations stayed in this folder. The independent parallel attempt was not read. No cluster calculation or network shell command was performed.
