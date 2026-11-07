# C138: The ferromagnetic crossover of the XXZ ring

Your folder is `../work_C138/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md`, then all of `codex/work_C136/` (RESULT.md, PROOF.md, REFEREE.md, an independent referee who found the proof correct), `codex/work_C134/` (RESULT.md, PROOF.md, REFEREE.md) and `codex/work_C129/` (RESULT.md and PROOF.md, in particular the ferromagnetic band and the dilute bound). Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Your sandbox cannot reach Typhon; write any large computation as a SLURM-ready script and say so.

## Setting

For the XXZ ring H = J sum (Sx Sx + Sy Sy + Delta Sz Sz) at minimizing flux, write F_L(n) for the sector ground energy, e_bulk(Delta, rho) for the bulk energy density at density rho, and epsilon = 1 + Delta, Delta = cos gamma, gamma = pi - eta. C136 proved, refereed, that for Delta in [-1 + delta, 0] and all n, L,

    |F_L(n) - J L e_bulk(Delta, n/L)| <= K_all(delta) J/L,

which gives the whole U(1) spin-1/2 class with J_z/J in [-1 + delta, 0]. The constant diverges as delta -> 0, because the kernel's L1 norm 2 gamma/pi - 1 tends to one and the Neumann inverse degenerates. C129 proved the band L Lambda/M_0 <= epsilon L^2, which covers L up to order epsilon^(-1/2). C136 showed the existing estimates start working only at L of order epsilon^(-1), with a coefficient diverging as epsilon^(-3/2), and identified the exact missing input:

    (X) there are A, K_* independent of Delta and density with
        |F_L(n) - J L e_bulk(Delta, n/L)| <= K_* J/L   whenever (1 + Delta) L^2 >= A.

A one-sided version (only the lower bound on F_L, or a direct bound on the doubling defect) also suffices. Given (X), the attractive U(1) class closes with C = max{A, 5 K_*} + pi^2/4.

## Physics of the crossover

Near Delta = -1 the velocity is v = (pi J/2) sin(gamma)/gamma, of order J eta, and epsilon is of order eta^2, so the crossover length is of order 1/eta = epsilon^(-1/2): below it the ring looks like the isotropic ferromagnet (quadratic magnons), above it like a Luttinger liquid with small velocity. The Casimir energy pi v/(6L) is of order eta/L, which is small, so the target is not delicate at the level of the answer. The difficulty is only that the C134/C136 machinery normalizes errors by an O(1) scale that is wrong near the ferromagnet.

## Tasks

1. Rescale. Find the scaling of rapidities near gamma = pi in which the Bethe equations and the kernel have a nondegenerate limit (the isotropic ferromagnet's rational parametrization is the natural candidate), and redo the C134/C136 estimates in those variables. The aim is to replace the Neumann bound 1/(1 - q) with an invertibility estimate for 1 + K on the occupied interval that is uniform as eta -> 0 once L eta is large. Track every constant's dependence on eta and L separately.
2. Prove (X), or the one-sided version, with explicit A and K_*, uniformly for Delta in (-1, -1 + delta_0] for some explicit delta_0 > 0, and at all densities (combining with C129's dilute bound where convenient). Note that only errors of order J/L are needed while the Casimir term itself is of order eta J/L, so there is room.
3. If the Bethe route stalls, try a direct variational and comparison argument in the window epsilon^(-1/2) << L << epsilon^(-1): for example perturbation theory in epsilon around the isotropic ferromagnet in each sector (whose sector ground energies are -J L/4 exactly), with the first-order term computed exactly and the second-order term bounded with the correct L dependence.
4. Assemble with C136 and C129 into a single theorem on the whole attractive U(1) class, J_z/J in [-1, 0], every field and DM angle and every L >= 2, with one explicit C.
5. Numerics: with the Bethe solvers of C134 and C136, compute L |F_L(n) - J L e_bulk|/J and L Lambda/M_0 for Delta = -1 + 10^(-k), k = 2 to 6, densities 1/L to 1/2, and L across the window from epsilon^(-1/2) to well beyond epsilon^(-1) (several thousand where feasible), and check every intermediate inequality.
