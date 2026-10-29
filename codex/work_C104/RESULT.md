# C104 result

[Proved, conditional on the WZW stiffness dictionary] There is a crossover-free polynomial bound for full simple WZW factors coupled to a conserved on-site semisimple symmetry,

    A_WZW <= (4 pi/3) d(d^2-1) M.

For only simply laced symmetry factors the coefficient improves to 2 pi/3. The proof permits any finite number of copies with different positive integer levels and velocities, including several copies coupled diagonally to the same microscopic simple factor. No upper bound on WZW levels is assumed or proved.

[Proved, limitation] The requested theorem for every charged U(1) infrared has not been proved. Charge quantization and locality of compact-boson vertices do not imply the proposed lower bound on Luttinger weights. This run does not provide a fixed-d nearest-neighbor lattice counterexample to such a lower bound either. The finite-d lattice question remains open.

## The WZW result

[Assumed] The counted sectors are full g_j,k WZW factors whose currents match the conserved microscopic g_j currents. Stiffness limits match their continuum flux response. A general CFT containing these currents can also have a coset stress tensor, which requires separate accounting.

[Proved from these inputs] Define the on-site index by tr_R(T^A T^B)=I_j delta_AB, with canonical generators normalized so that a long-root SU(2) has its usual spin matrices. For a spin-S representation, I_j=S(S+1)(2S+1)/3. For fundamental SU(N), I_j=1/2. With r_j=rank g_j and h_j its dual Coxeter number,

    S_j = sum over normalized Cartans D
        = r_j/(2 pi I_j) sum_{a in j} k_a v_a,

    A_j <= [2 pi I_j dim(g_j)/(r_j(1+h_j))] S_j.

All Cartans together obey sum_j S_j <= 4M. Restricting R to a long-root SU(2) gives I_j <= d(d^2-1)/12. The root-system identity dim(g_j)/[r_j(1+h_j)] <= 2 then proves the bound. The representation-sensitive version is

    A_WZW <= 8 pi M max_j
        { I_j dim(g_j)/[rank(g_j)(1+h_j^vee)] }.

[Proved] Large k helps this comparison, since c/k=dim(g)/(k+h^vee) decreases with k. Bounding the level by the representation content is unnecessary and is not supplied by the on-site Lie algebra alone.

## The missing U(1) step

[Proved within the compact-boson theory] For integer vertex labels m,n and charge Q=m,

    x_m,n = m^2/(4K) + n^2 K,
    spin_m,n = m n,
    D_Q = v K/pi.

The mutual locality condition is independent of K. Every K>0 is consistent with this local integer vertex lattice. The lightest charge-one vertex has dimension 1/(4K), with no upper bound from integrality. A finite on-site charge range does not establish such an upper bound on a scaling dimension.

[Proved] At most d integer values is weaker than a charge window of width d. Even imposing a window of width W<=d only bounds the normalization S_Q=tr((Q-tr(Q)/d)^2) by dW^2/4. The stiffness weight of the normalized charge is K/S_Q, so the missing lower bound on K remains.

[Proved, sufficient replacement] The exact surviving Gaussian inequality is

    sum_a v_a |Z_a|^2 <= 4 pi M.

A lower bound |Z_a|^2 >= 1/g(d) would suffice, but so would the weaker aggregate statement

    sum_a v_a <= G(d) sum_a v_a |Z_a|^2.

Neither has been established for all chains in the requested class. Nonzero coupling of each mode to a conserved charge is not a quantitative lower bound, and does not imply an independent microscopic charge for each mode.

## Neutral sectors and normalization corrections

[Proved] Independent neutral Ising, parafermion, and coset sectors contribute to A but not to the available charge stiffnesses. A d=4 product of a charged XX chain with coupling epsilon and a neutral critical Ising chain has D proportional to epsilon while the neutral contribution to A stays positive. Thus stiffness alone cannot bound the full A in the presence of an independent neutral factor. This does not refute an A/M bound.

[Proved, conditional on the full WZW identification] The parafermion part of SU(2)_k, including its Ising part at k=2, is covered when the full WZW current algebra and its common velocity are present. It is not covered merely by naming a Cartan U(1).

[Proved] The supplied pairing README needs three local corrections, documented without editing that folder. The torus double-commutator expectation has a minus sign. Its inverse-Z compressibility formula is inconsistent with its scalar convention; that convention gives kappa=Z^T V^(-1) Z/pi. For on-site trace-orthonormal WZW generators, the stiffness is k v/(2 pi I_R), not k v/pi unless I_R=1/2. The last issue changes higher-representation estimates.

[Numerical] Small deterministic checks passed for spin matrices through d=8, all 914 SU(2) representation decompositions through d=16, 401 Lie-algebra ratio checks, and the Heisenberg and spin-1 normalization comparisons. No heavy numerics or many-body simulations were run. These checks do not prove infrared matching.

[Proved, deliverables] Full derivations and source links are in PROOF.md. Reproduce the finite-dimensional checks with `python3 check.py`; output is in CHECKS.md. PROGRESS.md records the outcome, and STATE.md supplies local project memory. No files outside work_C104 were modified.
