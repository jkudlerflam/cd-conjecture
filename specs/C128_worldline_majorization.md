# C128: Crossed majorization for stoquastic chains in the world-line picture

Your folder is `../work_C128/`. Read `../../specs/_Cd_target.md` first, then `../../specs/C127_crossed_majorization.md` (a parallel run with the same target in operator language; do not wait for it), `../../specs/_bound_background.md` with all updates, `../../specs/_reconnection_lemma.md`, REPORT.md sections 9.80 to 9.88, `../work_C124/RESULT.md` and all of `../work_C125/`. Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; no network commands except a browsing tool for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Your sandbox cannot reach Typhon; write any large computation as a SLURM-ready script and say so.

## Setting

Take a stoquastic chain (nonpositive off-diagonal matrix elements of h in a product basis, possibly after an on-site unitary or a sublattice rotation) whose checkerboard transfer operator is positive (C125: XX, XXZ with Delta >= 0, the antiferromagnetic Heisenberg chain). In continuous imaginary time the path integral is a sum over world-line configurations with nonnegative weights. Cut the L x beta torus along one column. The column state is a decorated history (the site's state as a function of imaginary time, with each jump labelled by the bond that caused it), and the spatial transfer is a nonnegative kernel T(sigma, sigma') on such histories. Then Z(L, beta) = Tr T^L and R_L = log[(Tr T^L)^2/Tr T^(2L)].

## The route

The classical proof (C118, R_L <= log d) uses that T has rank at most d. Here the history space is infinite, but the events on a column form a point process in time whose intensity is of order M. The idea is a Poisson domination. Dominate the kernel T^L by an explicit positive kernel in which the history decomposes into independent events, each event costing a crossed decay that grows with its frequency content, so that Tr T^L is bounded by a free (Fock) crossed partition function with at most S(d) species per Matsubara frequency, of the form

    Z~(L) <= exp(K(d)) prod_(k>=1) [1 - exp(-2 pi k L/(V' beta))]^(-S(d)),   V' of order M,

which gives log Z~ <= K + S pi V' beta/(12 L) and A <= S V'/2. The two inequalities you need are an upper bound on Tr T^L and a lower bound on Tr T^(2L) (or a direct bound on the effective rank (Tr T^L)^2/Tr T^(2L)) that lose only exp(K + O(beta M/L)). Use the world-line structure and nonnegativity throughout. C124's commuting-chain proof (recorded histories as a nonnegative matrix) and §9.86's finding that, for stoquastic chains, the naive zero-cost step is a comparison of cross terms that fails numerically (ratio 1.055 for Heisenberg at beta = 6, L = 4) tell you which steps cannot be lossless.

## Tasks

1. Write the decorated-history transfer kernel explicitly for XX and for the Heisenberg chain, and verify numerically at small Trotter number that Tr T^L reproduces Z(L, beta) on even rings.
2. Prove a Poisson or Fock domination of the kind above for XX first, as a check against the exact crossed spectrum of C125 and against C119's constant. Then do Heisenberg and general stoquastic chains with a positive transfer. State precisely which property of the world-line measure (nonnegativity, the bounded intensity of events, the d-valued state between events, locality in time of the weights) each step uses.
3. Decide whether the frequency-dependent crossed decay can be obtained from the world-line representation directly, for example from a coupling of histories at distance L, or from a reflection or chessboard estimate in the spatial direction. C125's real-frequency bound is the operator version.
4. Test every inequality numerically on XX exactly and on Heisenberg and XXZ at Delta = 0.5 and 2 at small Trotter number. If the domination fails, identify the obstruction as a property of an explicit lattice chain, and state the minimal additional input.
