# C120 result

2026-10-29.

[Open] The universal estimate R_L(beta) <= K(d)+k M f(d) beta/L remains unproved. No counterexample to it or to A <= C M f(d) was found. The direct dimension-bounded-divergence replacements fail for identifiable reasons.

[Proved] Arbitrary-reference max-relative entropy, Petz and sandwiched Renyi divergences at every fixed order 0 < alpha < 1, and negative log fidelity are not bounded by Hilbert-space dimension. A faithful qubit pair suffices to make each arbitrarily large. Support min-relative entropy is instead zero between the finite-temperature ring Gibbs states, so its direct use loses the ratio information.

[Proved] Conditional min-entropy does give a dimension bound and an exact replacement identity. With X = exp(-beta H_(2L)), Y = exp(-beta(H_L tensor I+I tensor H_L)), rho = X/Tr X, sigma = Y/Tr Y, and a seam region C of dimension q,

    mu(W) = inf {Tr Omega_D : W <= I_C tensor Omega_D},
    R = Hmin(C|D)_sigma - Hmin(C|D)_rho + log[mu(Y)/mu(X)],
    R <= 2 log q + log[mu(Y)/mu(X)].

At fixed L the last term is beta[E_(2L)-2E_L]+O(1). All of the unknown slope remains there. The entropy capacity does not supply a bound on that normalization ratio.

[Proved] Petz and sandwiched Renyi divergences also give exact mixed-trace identities. For example,

    R = D_alpha(rho||sigma)
          + [log Tr(X^alpha Y^(1-alpha))-log Tr X]/(1-alpha).

Replacing the Umegaki divergence by D_alpha in the original forward identity gives a lower bound on R. The reverse identity gives a valid upper bound, but at any fixed order below one it still loses a factor proportional to L on XX. Fidelity similarly leaves an unnormalized overlap. Fuchs-van de Graaf does not bound its logarithm or remove the normalization term.

[Proved] Cone positivity, even strict positivity and complete positivity, does not imply the proposed quantum trace Cauchy-Schwarz inequality. A depolarized qubit rotation has superoperator dimension four and

    log[(Tr Phi^m)^2/Tr Phi^(2m)] = log[(1+a^m)/(1-a^m)],

which diverges as a approaches one. Its Birkhoff contraction coefficient is a. This is a counterexample to the abstract map inequality, not to the physical chain target.

[Proved] For a real Trotter transfer matrix, the classical proof extends with an explicit negative-pair penalty,

    R^(n) <= log D_n + log[1+Nminus/Tr B^2],
    B = T_n^(L/p),   Nminus = sum_(i != j) max(0,-B_ij B_ji).

Even a zero sign penalty leaves the full time-history dimension D_n = d^(pn). A beta M/L estimate still requires control of its effective dimension at separation L. The actual XX/Heisenberg transfer samples have zero penalty, so non-positivity is not the only obstruction.

[Numerical] All exact identities and valid comparisons passed the completed checks. Direct dimension ceilings and forward lower-order replacements failed on XX, Heisenberg and SU(3). Full ED covers lengths 4 and 8 for the qubit chains, and 3 and 6 for SU(3). Petz trace-identity residuals are below 2e-13. Twenty-four Trotter networks were checked against independent temporal contractions on both lengths. Full XX thermal traces reach paired length 16384. Ground-state free-fermion determinants test fidelity, Renyi divergences and conditional min-entropy through paired length 1024.

[Scope] No large finite-temperature conditional-min-entropy optimization or finite-temperature large-ring fidelity was performed. Those parts are covered by analytic inequalities, small-ring witnesses, and the explicitly identified zero-temperature determinant checks. No large local diagonalization, cluster submission, or independent referee review was performed.

[Work record] Derivations and primary references are in [PROOF.md](PROOF.md), failed candidates and numerical scope in [CHECKS.md](CHECKS.md), and reproduction in [checks.py](checks.py). The next useful input would be an independent beta M/L bound on one of the retained unnormalized quantities. Bounding only normalized-state distinguishability does not supply it.
