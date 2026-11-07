# C127: A free-boson majorization of the crossed spectrum

Your folder is `../work_C127/`. Read `../../specs/_Cd_target.md` first, then `../../specs/_bound_background.md` with all updates, `../../specs/_reconnection_lemma.md`, REPORT.md sections 9.80 to 9.88, and all of `../work_C125/` (RESULT.md, PROOF.md, CHECKS.md and its scripts). Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; no network commands except a browsing tool for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Your sandbox cannot reach Typhon; write any large computation as a SLURM-ready script and say so.

## Setting

For chains with a positive spatial transfer operator (the checkerboard construction of C125, valid for XX, XXZ with Delta >= 0, the antiferromagnetic Heisenberg chain at any spin and SU(2) Sym^k), write T_beta = exp(-H~) with crossed spectrum eps_j and crossed momenta p_j (eigenvalues of the Trotter-circle rotation), and Z~(L) = sum_j exp(-L (eps_j - eps_0)). Then log Z~(L) <= R_L <= 2 log Z~(L). The target is

    log Z~(L) <= K(d) + kappa(d) M beta/L   for all even L >= L_0 and all beta,

which gives A <= (6/pi) kappa M and hence a finite C(d).

## What C125 established, and the gap

Positivity, unitary Trotter rotations, an entire tilt, the cone condition eps - eps_0 >= |p|/V and the trivial bounds do not suffice. The n-boson countermodel satisfies all of them with A = nv. Its defect is a multiplicity, n independent oscillators per crossed momentum. On the lattice, the crossed states reached by one single-site insertion O in M_d at Matsubara frequency omega span at most d^2 - 1 directions per frequency (traceless O). C125 also proved the low-frequency decay |G_AB(i omega, x)| <= 2 ||A|| ||B|| V exp(-omega x/V)/[omega (V - omega)] for 0 < omega < V = 4eM, while the global cone fails at high frequency (XX has crossed energies growing like log m).

## Heuristic to make rigorous or refute

Suppose the low crossed spectrum were dominated by a Fock space of single-insertion excitations, with at most S(d) species per Matsubara frequency omega_k = 2 pi k/beta (counting both signs of the crossed momentum) and crossed energy at least omega_k/V each. Then

    Z~(L) <= exp(K) prod_(k>=1) [1 - exp(-2 pi k L/(V beta))]^(-S(d)),
    log Z~(L) <= K + S(d) pi V beta/(12 L) + (lower order),

so A <= S(d) V/2, of order 4e(d^2 - 1) M, a polynomial C(d). Normalization check: one c = 1 boson of velocity v has two oscillators per k and log Z~ = pi v beta/(6L). Also check against XX's exact crossed spectrum, eps_m - eps_0 = asinh(pi/beta) + asinh((2m - 1) pi/beta) in each fermion channel (C125, PROOF section 3). Note what the constant K must absorb. Light crossed states at zero crossed momentum (low-dimension primaries, light charge sectors near a ferromagnetic point) decay more slowly than the majorant at large L. For U(1) sectors their total is controlled by the stiffness, sum over charges of exp(-Q^2/(2 rho_s beta/L)), with sum_a rho_(s,a) <= 4M (Kohn budget). The majorization should therefore be stated with an explicit K and possibly a separate zero-frequency factor.

## Tasks

1. Make the single-insertion sector precise: the span of crossed states created by one insertion at each Matsubara frequency, the lower bound on their crossed energies from the C125 decay estimate (with all constants and its domain omega < V), and the multiplicity per frequency.
2. Multi-insertion states. Decide whether a majorization of the above form (any V' of order M, any S(d), any explicit K(d)) holds for every chain with a positive transfer operator. Candidate tools: a polymer or cluster expansion of the Trotter column in insertions with a crossed decay factor per polymer; a Gibbs–Peierls–Bogoliubov comparison of the crossed Gibbs state with a quasi-free crossed reference built from the insertion two-point functions; band-limiting (crossed states of crossed energy at most eps are concentrated on frequencies below V eps, then a Landau–Pollak–Slepian count of time-band-limited histories). A count that is only linear in eps (log N(eps) <= K + c beta eps) controls Z~ only for L > c beta, so the regime L << V beta needs the product (square-root, Cardy-type) structure. Say explicitly where interactions enter.
3. Test numerically every inequality you propose. Use XX exactly, and the Heisenberg chain and XXZ at Delta = 0.5 and 2 with a numerical transfer matrix at small Trotter number, resolving the crossed spectrum by crossed momentum and comparing it with the bosonic majorant.
4. If the majorization fails in general, construct the failure as a lattice chain or a lattice-realizable crossed structure (not an abstract spectral model), and identify the minimal extra lattice input that restores it.
