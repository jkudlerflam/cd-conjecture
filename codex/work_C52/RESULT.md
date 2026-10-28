# C52 result

[proved, outcome] S has not been reduced to a published theorem. The direct exact-invariant-cell redesign is obstructed, and the unchanged C44 mobile wrapper also fails exact invariance. No successful fixed-d nearest-neighbor replacement bond was obtained. The original conditional construction is not refuted.

[proved, main theorem] Suppose a physical nearest-neighbor chain has an exactly invariant tensor product of cell codes, with distinct physical left and right endpoints in each cell. Its neighboring connected logical bonds must commute. This is C8 RESULT Section 10's theorem, whose proof is reproduced with attribution in [PROOF.md](PROOF.md), Section 2.

For a number-conserving logical bond with nearest-cell hopping T, the three-cell one-particle matrix element is

    <particle in cell 1|[K_12,K_23]|particle in cell 3> = T^2.

Consequently an invertible T is impossible in this class. Vacuum-normal-ordered quartic interactions vanish on the one-particle sector and cannot repair the contradiction. The result uses invariance of the full Hamiltonian, so internal counterterms cannot cancel the extracted leakage.

[proved, hypotheses checked for the target] C39's hopping is

    T_0 = -W[Gamma(w)+i eta Gamma(wv)],
    sigma_min(T_0) >= W lambda_min(Gamma(w)) > 0.

The positive Gram matrix proves invertibility for every n. A bond error smaller than half the displayed bound preserves it. This tolerance can be added to C42's finite all-n parameter schedule. Thus an exact product-cell realization cannot preserve that target with sufficiently small corrections, regardless of the fixed physical dimension and regardless of whether the parent consists of commuting projectors or is merely frustration free.

[proved, mobile extension] Exact invariance of the direct sum over C44's marker patterns would force the diagonal Hamiltonian block of each fixed pattern to preserve its product code. Marker-changing terms have orthogonal outputs and cannot cancel leakage in that block. The same obstruction therefore applies to C44's common payload hopping, or sufficiently close length-dependent versions.

[proved, independent wrapper obstruction] For the literal C44 marker move, write

    phi_l(i)=sqrt(2/l) sin(pi i/l),
    A_b=sum_(i=1)^(b-1) phi_b(i) phi_(b+1)(i),
    s_b=1-2 sin(pi/(b+1))^2/(b+1).

Its head output has squared norm s_b and projected amplitude A_b^2. The omitted component has squared norm

    s_b-A_b^4 >= s_b(1-s_b)>0.

The other wrapper factors give a nonzero multiplier, so nonzero physical marker hopping leaks out of the selected code even with the payload ports off. C44's exact compressed identity transport remains correct. It is not an invariant-subspace statement. Section 5 also proves that a boundary-local repair retaining the same sine-head family fails whenever two interior head positions remain untouched, because their required amplitude ratios change differently when the cell grows.

[proved, further restrictions] Direct compression through a d-dimensional physical bond gives 2 rank(T)<=d^2, or rank(T)<=floor(d^2/4) with charge-preserving encoding. Independently, an n-orbital nearest-neighbor quadratic Hamiltonian with rank-r hopping has at most 2r isolated Fermi roots. Rank reduction therefore cannot keep C39/C42's 2n simple roots. The already settled C1 bound for quasi-free nearest-neighbor chains, recorded in project STATE and STATUS, also excludes an all-n fully Gaussian physical-bus replacement.

[proved, exactly local dressing qualification] A cross-boundary dressing with physical locality radius R is still subject to the commuting-bond theorem when all cells have length greater than 4R+2. To evade this argument for b=n+3, its range must grow at least on the order of n. An arbitrary dressing confined separately inside each cell is absorbed into the code and never evades the theorem.

[assumed, open alternatives] The task allows a more general redesign. An encoding entangled across cell boundaries, a dressing whose physical range is comparable to the cell length, or a different mobile selector may evade these proofs. No complete construction in that class was found. These results are not a no-go theorem for all fixed-d chains, for all exact encodings, or for the proposed infrared phase.

[proved, theorem-application boundary] No Mastropietro-Porta or Porta-Scola-Singh theorem is invoked. There is no exact physical fermion model here with its hypotheses checked. The missing construction precedes the RG question. The full original S remains open, including microscopic matching and the infrared identification.

[numerical, checks] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes finite algebra controls using matrices of dimension at most 64. The adjacent-bond commutator identity and quartic extraction errors are zero in the recorded run; the determinant and head identities agree to roundoff. At b=4 the head leakage squared is 0.1701573216462289. These checks do not measure a central charge or simulate the microscopic mobile chain.

[proved, files and scope] Full proofs and limitations are in [PROOF.md](PROOF.md). [PROGRESS.md](PROGRESS.md) and [STATE.md](STATE.md) give the handoff, with [checks.py](checks.py) and [checks.json](checks.json) for reproduction. No network, subagents, cluster computation or files outside work_C52 were used. The project STATE was not modified because it is outside the permitted write scope. No independent human review is claimed.
