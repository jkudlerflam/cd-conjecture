# C131 (draft): A compactness cover for C(2)

Draft spec from `contrib/rigorous_toolbox_map/`, for coordinator review. Paths are repository-relative. Your folder would be `codex/work_C131/`. Read `specs/_Cd_target.md`, `specs/_bound_background.md` with all updates, REPORT.md sections 9.73, 9.77 (C107) and 9.88, `contrib/rigorous_toolbox_map/TOOLBOX.md` Section 8, and the drafts `DRAFT_C129_d2_u1_bethe.md` and `DRAFT_C130_constructive_rg_casimir.md`. Rules of `_bound_background.md` apply. Large computations go in a SLURM-ready script, not submitted from the sandbox.

## Setting

Jonah has said that C(2) < infinity is a result he wants. The general spin-1/2 bond is XYZ (J_x, J_y, J_z) plus a DM vector D plus a uniform field b, up to on-site SO(3) rotation and gauge (TOOLBOX Section 8.1), nine parameters, eight after fixing M = 1. C107 showed the gauge quotient of the unit sphere is compact.

## Reduction lemma (to be proved first, it should be short)

For each fixed L, Q_L(h) = L(E_{2L} - 2E_L)/M(h) is continuous on the compact set S = {h : M(h) = 1}/gauge. Hence C(2) < infinity holds if and only if every h_0 in S has a neighbourhood U and a constant C_0 with Q_L(h) <= C_0 for all h in U and all L >= 2 in the chosen vacuum-compatible set. State the lemma with the vacuum-compatibility convention made explicit. In the U(1) class all L should be admissible (C129 Step 2). For general bonds, decide whether all even L, or all L, is the right convention, using C107's ferromagnetic XX example as a warning that the raw limit depends on the length subsequence.

## Strata of S and their status

Write a table classifying S, with for each stratum the local lemma required and the tool that might supply it.

(a) The U(1) class XXZ + DM_z + b_z: C129. Its neighbourhood needs a U(1)-breaking perturbation argument, which is C130's job only near free points.

(b) Quasi-free bonds (J_z = 0 in a JW frame): exact finite size (C130 free leg), and their neighbourhoods by constructive RG where nondegenerate.

(c) Gapped bonds. The required local lemma: if h_0 has ring ground states with a gap at least gamma > 0 above a ground space of fixed dimension q for all L >= L_0, then there are a neighbourhood U and C, xi with |E_L(h) - e_0(h) L| <= C e^{-L/xi} for h in U, so Q_L is uniformly small. Two separate inputs are needed: (i) exponential finite-size convergence of the ring energy at a gapped point (expected from Lieb–Robinson plus gap, via Hastings–Koma clustering, math-ph/0507008, and quasi-adiabatic continuation; I did not find this stated for frustrated chains, so state precisely what is proved and what is assumed); (ii) stability of the gap in a neighbourhood, uniformly in L, which is known for perturbations of classical or frustration-free systems with local topological order (Borgs–Kotecký–Ueltschi CMP 181 (1996) 409, Datta–Fernández–Fröhlich J. Stat. Phys. 84 (1996) 455, and the frustration-free stability literature, recalled and not verified in the survey) and appears open for general frustrated gapped chains. Near classical bonds with a Peierls condition, Borgs–Kotecký Theorem 4.1 gives both at once.

(d) Non-integrable critical surfaces far from free points (Ising-type surfaces of XYZ in a field, off-axis DM). No rigorous tool. Record how much of these surfaces lies within the C130 ball of some free point, as a numerical estimate.

(e) Multicritical points: the Heisenberg point (SU(2)_1 with a marginally irrelevant operator), the ferromagnetic point Delta = -1 (z = 2), the free Lifshitz points where Ising lines meet the U(1) line (XY at the saturation field), and any points where gapped and critical strata meet with diverging crossover length. Each needs its own local lemma. Note that only an upper bound on Q_L is needed, and that Q_L at Heisenberg decreases (4, 2.79, 2.62, 2.56 at L = 2, 4, 6, 8).

## First concrete step

Prove the reduction lemma. Then produce the stratification table with, for each stratum, the precise local statement, the literature theorem that would supply it (from TOOLBOX), and a status tag. Do not attempt (d) or (e) analytically in this run.

## Falsification test (numerical, small and decisive about where to look)

Maximize Q_L over S by multistart local optimization (Nelder–Mead or BFGS on the eight angles, gauge-optimized M computed by convex minimization as in `contrib/rigorous_toolbox_map/d2_sanity.py`) separately at L = 2, 4, 6, 8 (rings to 16 sites, sparse exact diagonalization), even and odd L recorded separately. The survey's 48-point sample never exceeded the Heisenberg value Q_2 = 4. Report sup_S Q_L for each L and the location of the maximizer. If sup_S Q_L grows with L, or the maximizers drift toward a stratum as L grows, that stratum is where C(2) could fail and must be attacked first. If the maximizer stays at Heisenberg with Q_L decreasing, record that as evidence that C(2) is attained at L = 2 by Heisenberg, which would be a striking sharp conjecture, C(2) = 4.

## What success means

A precise reduction of C(2) < infinity to a finite list of local lemmas, each tied to a stratum and a tool, with (a) and (b) handed to C129 and C130 and (c) stated sharply enough to be attacked next. Together with the numerical maximization this fixes the order of attack on C(2).
