# Referee report on C39, the two-bus core with a geometric readout

## Claim

Under review is `codex/work_C39/`. It claims a 46-letter positive initialized core (exact U(1), complete frustration-free 2^n code, two folded buses with stub hoppings w, wv and port phases pi/2, pi, and a biased vacuum walker whose readout supplies a geometric field) whose second-order chain, with limiting reference h(k) = -2A cos k + 2 eta B sin k + mu E_u, has for every n simple Fermi points, distinct centres and positive margins against every number-conserving quartic resonance. Parameters come from a terminating dyadic search and transfer to finite D, g, lambda. The inherited C8/C32 wrapper gives d0 = 4464.

My verdict. I found no fatal error. The bond formula, its normalization, the Newton-polygon lemma, the slopes and the pair-sum classification check out by hand, an independent reimplementation reproduces the bond to machine precision, and high-precision checks through n = 8 find positive margins at the certified tuple and along the search. Unbounded c still rests on C8's wrapper (inherited, unaudited) and on the assumed thermodynamic elimination. The core is not compatible with the one-vacancy completion as it stands.

## Status

Refereed. The finite all-n claim holds at the repository's standard, up to the fixable items in Section 3. Two things would upgrade the mobile conclusion, an audit of C8's exclusion and gluing for this core and the infrared step of C40 and `contrib/ir_control/`.

## How it was established and checked

- Author: Claude (Opus 5.5) referee subagent for Jonah Kudler-Flam, 2026-10-28.
- AI assistance: none beyond this session. I used no Codex, because each check is a short script whose specification would be as long as the code.
- Checked by: me. I rederived PROOFS Sections 2 to 5 and reran no C39 script. The bond was rebuilt from the rule table in the physical, ungauged basis and evaluated by sparse LU. The Gram matrix comes from a closed eigenmode formula, independent of the Lyapunov route and of eq. (19).

## Argument or method

### 1. Hand check of the proof

- Bond formula. Eq. (6) is second-order degenerate perturbation theory at zero code energy, with visibility q^(2n+2)/S^2 from the vacuum amplitudes at B and A. Bus type is conserved, so buses do not mix, and a second insertion on another cut leaves an outer cell excited, so there is no three-cell term. The g^2 phases give T = -J[Gamma(w) + i eta Gamma(wv)], and the relative i is a gauge-invariant loop phase.
- Lemma (13), (14) and slopes (16). The least-exponent minor is unique, so r_j ~ kappa_j v^(4j). Components i < j of x_j are O(v^(2(j-i))), so the slope hierarchy needs v small compared with u, which the order u then v allows. Eq. (16) is the differentiated zero-mode equation.
- Pair sums. With exact U(1) every quartic resonance has the form p_a + p_b = p_c + p_d mod 2 pi, and bilinear ones are the case of a shared label. The zero-field difference is a lacunary centre combination plus an integer multiple of pi. Odd multiples miss by at least pi - 1, nonzero flavour imbalance is caught by the 7/8 lacunarity bound, and equal flavour multisets force a nonzero even slope vector of l1 norm at most four. This covers pair transfer, the exchange {(a,+),(b,-)} against {(a,-),(b,+)} and every umklapp. At mu = 0 the identity h(k+pi) = -h(k) puts every flavour at half filling, so the exchange and all 4k_F umklapps resonate there and only the readout splits them. The minimum margin is the slowest flavour's 4k_F umklapp, delta4 = 4 mu s_(n-1) to leading order. The classification is complete.
- Transfer and termination. Conversions sit only on the data rail, so every correction to T_ij carries w^(2i+2j). The order (w, eta, u, v, mu), tolerance, D, g, lambda has no circularity, the core gap only improves as q = (wu)^2 shrinks, and positivity needs only |g| < a - 3. Each search test holds for all small enough dyadic values, so each step stops, but no depth bound is given (problem 3).

### 2. Reimplementation and numerics

| Check | Result |
| :-- | :-- |
| Bond at the diagnostic point, n = 1, 2, 3 (`ref_bond.py`) | Largest entry difference from `bond_n*.npz` 5.6e-17, 5.6e-17, 6.6e-17. Largest hopping agrees to 1.7e-12, 5.5e-11, 9.1e-10 relative |
| g^2 hopping against -2(q^(2n+2)/S^2)[G(w) + i eta G(wv)], with G the two-token Gram of the full folded bus, D = 0.2 and 2 | Relative error at most 3.5e-9 |
| D^(2n+4) G against Gamma(w), n = 1 to 4, D up to 10^6 | Entrywise relative error (22.5, 38.0, 53.5, 68.9)/D |
| g-series of the bond (`ref_gauss.py`) | Order g^0 scalar. Order g^2 quadratic to 1e-21 with no one-cell part. Order g^4 carries the first one-cell and non-Gaussian terms |

The high-precision reference (`ref_pencil.py`, 500 to 1300 decimal digits) brackets the 2n roots by exact sign changes of the degree-2n polynomial, one root per bracket, and refines them.

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: |
| delta4, certified tuple | 5.1e-16 | 1.2e-28 | 2.4e-41 | 4.2e-54 | 7.1e-67 | 1.2e-79 | 1.9e-92 | 3.1e-105 |
| delta4, search | 2.1 | 1.1e-3 | 5.0e-7 | 2.5e-10 | 9.2e-12 | 3.5e-16 | 8.0e-19 | 1.9e-21 |

The certified row reproduces C39's certificates at n ≤ 3. Along it lacunarity holds, delta4 is 16 times the bound (17), and r_j/(kappa_j v^(4j)) = 1 to within 6e-6. The n = 8 value is unchanged between 900 and 1300 digits. The search stops at v = 1/4, u = 1/16 or 1/32, eta = 1/32 and mu between 2^-8 and 2^-20, the first or nearly the first trials.

### 3. Problems, ranked

| # | Problem | Rank |
| :-- | :-- | :-- |
| 1 | The mobile conclusion inherits C8's filter, exclusion and gluing and C32's equalization bound, none re-audited. Their hypotheses hold (Section 4) | Serious, flagged by C39 |
| 2 | Massive-sector elimination and Gaussian stability remain assumed | Serious, flagged by C39 |
| 3 | No explicit depth bound for the search, since the o(1) in (13), (14) and the bounds L, R are not quantified. A pure interval test would also loop on an exact coincidence at a trial mu. Dovetailing fixes this | Fixable |
| 4 | The cluster ED cannot see the hopping. At lambda = 0.1 its residual exceeds lambda^2 max abs(T) by about 4e3, 9e4, 2e6 at n = 1, 2, 3, so it tests only the scalar bond and readout. The direct bond computations, which I confirm, carry the hopping | Fixable |
| 5 | Scale budget. The certified tuple loses about 13 decades of delta4 per flavour, the search about 3, before the factors w^(4j) in velocities and q^(2n+2) D^(-2n-4) in J. At n = 3 the g^4 one-cell terms over the g^2 hopping grow as 1.7, 6.8, 40 at D = 2, 8, 32, as PROOFS warns. Not an error, but the infrared step should use the search tuple | Fixable |
| 6 | `literal.py`'s 2116-dimensional bond penalizes every (B, A) adjacency, so it is the bond of one initialized cell. The chain bond with inter-cell grammar and port lambda V is never written as one matrix | Cosmetic |
| 7 | The express-phase direction is unstated. I took it toward B, as `literal.py` does, and the other choice leaves every margin unchanged. At the diagnostic point the n = 3 non-Gaussian norm (7.0e-11) exceeds the largest hopping (6.3e-11) | Cosmetic |

### 4. Alphabet, wrapper and one-vacancy completion

The count 4(3·2·46 + 3)·4 = 4464 is right and follows C32's 2736 pattern. C32's wrapper needs positivity, an initialized complete code, local endpoint ports and a core gap, and C39 supplies all four. Holes only switch bond terms off, which keeps positivity, the q^l kernel and the readout compression in every vacancy pattern, and the helix compensator acts on both buses through the shared endpoints. The wrapper errors must sit below a tolerance of order delta4 times the scaled velocities, finite at each n.

The one-vacancy completion of `contrib/counterexample_campaign/doping/` does not apply as it stands. The 23-letter chain carries its own selector, a frustrated walker, while C39's core is frustration-free on every cell length and on tokenless L or R domains, so a 47-state completion selects no length, and the order-lambda^2 scalar bond (-0.066, -0.033, -0.018 at n = 1, 2, 3) favours the shortest cells. A selector has to be added. That is serious for the doping route, not a defect of C39's claims.

## Reproducing

numpy and scipy only. From this folder, `python3 -B ref_bond.py` (5 s), `ref_largeD.py` (10 s), `ref_gauss.py` and `ref_g4scale.py` (seconds), `ref_pencil.py fixed 1,2,3,4 500` (10 s), `ref_pencil.py fixed n 900` for n = 5 to 8 (40 to 200 s) and `ref_pencil.py search n 700` (up to 90 s). `ref_bond.py` reads `../../codex/work_C39/bond_n*.npz` for the comparison. Outputs are the matching `.json` files.

## Relation to the rest of the repository

This reviews `codex/work_C39/` in the manner of `contrib/redteam_routeA_cells/`. STATUS's route-A line on C39 could read as follows. Refereed with no fatal error, the finite all-n band input holds for this core, and the mobile c = n + 1 conclusion remains conditional on C8's wrapper and on the infrared step. The doping route would need a selector added to this core.
