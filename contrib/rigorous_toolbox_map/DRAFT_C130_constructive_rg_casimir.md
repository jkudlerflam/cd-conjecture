# C130 (draft): Casimir bounds in a constructive-RG neighbourhood of quasi-free chains

Draft spec from `contrib/rigorous_toolbox_map/`, for coordinator review. Paths are repository-relative. Your folder would be `codex/work_C130/`. Read `specs/_Cd_target.md`, `specs/_bound_background.md` with all updates, REPORT.md sections 9.71, 9.80, 9.81 and 9.88, `contrib/rigorous_toolbox_map/TOOLBOX.md` (Sections 5, 7, 8), and `codex/work_C128/RESULT.md` (the XX finite-size bound). Rules of `_bound_background.md` apply. Large computations go in a SLURM-ready script, not submitted from the sandbox.

## Why this route

The survey found one rigorous technology that produces the Casimir coefficient for non-integrable interacting systems: constructive fermionic renormalization group. Giuliani–Mastropietro, CMP 324 (2013) 179, arXiv:1211.3030, prove for the 2D Ising model with nearest-neighbour plus weak finite-range interactions, 0 < lambda < lambda_0, that at criticality the finite-size correction to the torus free energy is universal with central charge 1/2. Benfatto–Mastropietro, J. Stat. Phys. 138 (2010) 1084, arXiv:0907.2837, control quantum spin chains H = -sum[J_1 S^1 S^1 + J_2 S^2 S^2] - h sum S^3 + lambda sum v(x - y) S^3 S^3 at small lambda, including the renormalized sound velocity v_s = v_F + O(lambda). Inside such a ball every scale from the lattice to L is controlled, so the long-crossover obstructions (lesson 3) cannot occur, and the result is a rigorous finite C on an open set of non-integrable chains. It does not reach the general target. Its role is a stratum in a cover (C131) and a test of whether constructive methods give the L-uniform finite-size statement, which no one seems to have written for quantum rings.

## Target

Let h_0 be a quasi-free nearest-neighbour bond with m modes per site (d = 2^m), critical with finitely many simple Fermi points of nonzero velocity (a nondegenerate quasi-free critical chain), and v any parity-even bond with ||v|| <= 1. Prove: there are lambda_0(h_0) > 0 and C(h_0) < infinity such that for |lambda| <= lambda_0 and all even L,

    L (E_{2L} - 2 E_L) <= C(h_0) M(h_0 + lambda v),

and, as the sharper statement the method should give,

    E_L = e_0(lambda) L - (pi/(6L)) A(lambda) + (charge and twist terms of order v/L) + O(L^{-1-theta}),

with A(lambda) = sum_a c_a v_a(lambda), c_a fixed by h_0, uniformly in |lambda| <= lambda_0, in both Jordan–Wigner parity sectors.

Two concrete instances.
- d = 2, Ising surface. h_0 = XY with anisotropy gamma != 0 at its critical transverse field (c = 1/2), v = Sz Sz (in the JW frame this is the density–density term). h_0 + lambda v is XYZ in a field, non-integrable for lambda != 0. This is the quantum-chain analogue of Giuliani–Mastropietro and covers part of the non-integrable critical locus of the d = 2 cover.
- d = 4, Luttinger case. Two fermion modes per site with a non-integrable nearest-neighbour density–density coupling between the modes, h_0 two XX chains with different hoppings (two Fermi velocities). This is the smallest non-integrable nearest-neighbour Luttinger example (at d = 2 the U(1) nearest-neighbour class is integrable, TOOLBOX Section 8.2).

## The free-fermion leg (do first)

[Assessment] The repo's uniform finite-L ratio bounds cover XX, independent XX stacks and the C106 family (C118, C119, C128). A uniform statement for every nondegenerate quasi-free nearest-neighbour chain, with pairing and arbitrary bands, does not appear to be recorded. Prove

    L (E_{2L} - 2 E_L) <= C_free sum over Fermi points |v_F|,

for every quasi-free chain and every even L, with the parity sectors handled as in C128, using the exact single-particle sums and the convex quadrature argument of C128 band by band. The free check of `DRAFT_C129_d2_u1_bethe.md` (L Lambda = (pi J/4) sin(pi rho)(1 + O(L^{-2})) for XX at density rho) is the one-band case. Degenerate free points (a band with a quadratic touching, z = 2) must be identified and excluded or treated separately, since they are where lambda_0 will vanish.

## The interacting leg

1. Write Z(L, beta) for the ring, in each parity sector, as a Grassmann integral on the L x beta torus with the appropriate spatial boundary condition (JW gives periodic or antiperiodic according to parity) and antiperiodic imaginary time. Identify the free covariance and its decomposition into scales.
2. Check what in Giuliani–Mastropietro (and Greenblatt arXiv:1409.1852 for cylinders, Antinucci–Giuliani–Greenblatt for cylinder energy correlations) is uniform in the aspect ratio beta/L. The zero-temperature ring energy is the beta -> infinity limit at fixed L, the infinitely long cylinder, so the multiscale bounds must hold uniformly as the time extent goes to infinity. State exactly which estimate is needed: the effective potential at scale h with infrared cutoff at 2 pi/L in space and none in time.
3. Show that the flow of the running couplings is controlled down to scale 1/L with bounds independent of L, so that E_L = e_0 L + (scale-1/L effective energy), and that the scale-1/L effective theory is a Gaussian (Luttinger or Majorana) theory with renormalized velocities, whose ring energy is computed exactly. For the Ising instance the quartic term is irrelevant and this is the GM mechanism. For the Luttinger instance the quartic coupling is marginal and flows to a line of fixed points, which is the Benfatto–Mastropietro mechanism.
4. Track the M-dependence: the final constant must be stated as C(h_0) M, with M the gauge-optimized half-spread of the full bond.

## First concrete step

Read Giuliani–Mastropietro arXiv:1211.3030 and Benfatto–Mastropietro arXiv:0907.2837 (and the XYZ paper they cite, Rev. Math. Phys. 13 (2001) 1323, if obtainable) and record in PROOF.md: (i) the precise finite-size theorem in GM, its geometry (torus M x N, which limits), and its error term; (ii) whether either paper treats a quantum chain on a finite ring at zero temperature; (iii) the dependence of lambda_0 on the free model's data (Fermi velocity, mass, distance to degeneracy). That reading decides whether the Ising instance is a transcription or a new constructive proof.

## Falsification tests

1. Ising instance. For XYZ in a field near the Ising surface, h = (1 + gamma) Sx Sx + (1 - gamma) Sy Sy + lambda Sz Sz - g Sz, with gamma = 0.5 and lambda in {0, 0.05, 0.1, 0.2}, locate the critical field g_c(lambda) at each L by finite-size gap crossing, then compute Q_L = L(E_{2L} - 2E_L)/M by exact diagonalization for 2L <= 20 (no U(1), so sparse full-space, parity-resolved). The prediction is Q_L -> (pi/4)(1/2) v(lambda)/M with v(lambda) smooth in lambda. A non-monotone or growing Q_L as lambda increases at fixed gamma would indicate lambda_0 is tiny.
2. Free leg. On 200 random quasi-free nearest-neighbour bonds (m = 1, 2, with pairing), compute L(E_{2L} - 2E_L) / sum|v_F| exactly to L = 4096 and check the proposed C_free, with particular attention to bonds near a quadratic band touching.
3. Luttinger instance. Two XX chains with hoppings 1 and t, coupled by lambda n_1 n_2 on each site and each bond, d = 4, exact diagonalization to 2L = 12 and, if needed, a DMRG script for Typhon. Q_L should stay below the free value plus O(lambda).

## What success means

A rigorous theorem: finite C on an explicit open neighbourhood of every nondegenerate quasi-free critical chain, with the free-fermion bound uniform over all quasi-free chains. Partial success, the free leg alone, closes the quasi-free stratum of C131 and upgrades C119 and C128 from examples to a class.
