# C116 result

2026-10-29.

[Open] I did not prove A <= C M f(d) for general nearest-neighbor chains, for any finite explicit C and f. I found no counterexample to that target. The full-spectrum investigation identifies a suitable ratio and several obstructions to the proposed ways of bounding it. Derivations are in [PROOF.md](PROOF.md), with numerical details in [CHECKS.md](CHECKS.md).

[Proved conditional on the infrared trace limit] The equal-temperature ratio retains the desired information:

    R_L(beta) = log[Z(L,beta)^2/Z(2L,beta)],
    A = (4/pi) lim_(tau -> infinity) lim_(L -> infinity)
                       R_L(tau L)/tau.

Lengths must be vacuum compatible. A crossover-free estimate

    R_L(beta) <= K(d) + k M f(d) beta/L

would imply A <= (4k/pi) M f(d). This implication is proved; the required lattice estimate is not. The exact all-scale bound obtained is only

    |R_L(beta)| <= 4 beta M.

It loses a factor L in the scaling limit, as the ground-energy gluing route does.

[Proved] The unweighted ratio Z(L,beta)/Z(2L,2 beta) cancels the torus function and retains the bulk term. The corrected combination Q = 4 log Z(L,beta) - log Z(2L,2 beta) retains three times the torus function, but also retains its reciprocal-velocity contribution. Two decoupled rescaled XX sectors have fixed d = 4, M = 1/2 and A = 1 while B = sum c_a/v_a diverges. Their Q diverges at any fixed aspect ratio in the scaling limit. R cancels this thermal bulk term.

[Proved] A temperature-independent spatial transfer dimension D <= d^2 is impossible. For a positive transfer matrix, Cauchy-Schwarz would give R <= log D. More generally, any exact finite transfer matrix representing every ring trace with the correct pressure must satisfy, already for XX,

    log D_beta >= beta (sqrt(2)-4/pi) - 4 log 2.

At beta = 32 this exceeds log 4. This stronger obstruction does not require transfer positivity or CFT input. The elementary gate has at most d^2 operator Schmidt rank, but the spatial transfer space retains the entire Trotter history, whose dimension grows with Trotter depth.

[Proved and assumed published input] Rotating the indices of a positive two-site gate does not preserve positivity. Reflection positivity requires additional coupling conditions, as stated in [Jaffe and Janssens, Theorem 5.2](https://arxiv.org/html/1506.04197v2#S5.SS4). A twisted XX chain also gives an explicit even-ring example with R < 0. Thus positivity of R is not universal at arbitrary lengths. Golden-Thompson gives a valid two-layer bound but leaves extensive bulk slack.

[Proved under stated massive-sector hypotheses] A spectator with finite vacuum multiplicity, a fixed gap and exponentially small ring vacuum corrections contributes only log g_vac to the scaling limit of R at sufficiently large tau. Its contribution to R/tau vanishes. The exact massive Ising spectator obeys 0 <= R_massive <= log 2 at every L,beta. Nearly critical massive states can affect finite-size traces; the ordered infrared slope removes them without an assumption of a uniform crossover scale.

[Numerical] Full-spectrum checks on XX and Heisenberg through L = 8 and SU(3) through L = 6 satisfy the gluing bound and violate R <= log(d^2). The complete XX trace agrees with direct diagonalization to 5.7e-14. At L = 256, tau = 8, the R and Q estimators give 0.9994309 and 1.0000236 against A = 1. The small Heisenberg and SU(3) rings have substantial finite-size corrections, and no extrapolated A is claimed. A two-velocity XX test and a nonsymmetric-bond comparison also pass.

[Blocked computation] Typhon SSH timed out during banner exchange on the initial attempt and on a later retry. No cluster job was submitted. Larger full-spectrum tests are prepared in `run_typhon.slurm`, with a SLURM guard in the Python script, but remain unperformed. No heavy numerical computation ran locally.

[Open next step] The missing input is a bound on the effective spatial transfer rank at separation L with logarithm of order M beta f(d)/L, allowing a sublinear-in-tau remainder. A bound on elementary gate rank, a qualitative reflection inequality, or the dual coefficient B does not supply it. No general exclusion of other partition-function approaches is claimed.
