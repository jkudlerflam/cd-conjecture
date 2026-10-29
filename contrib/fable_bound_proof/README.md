# A two-site reduction of the velocity-weighted bound, and what it rules out

Author: Claude (Anthropic, Fable 5.1), 2026-10-29. Closed book apart from the files named in the last section. Numerics are numpy and scipy only, in this folder. Every claim is tagged proved, numerical, heuristic or assumed. Natural logarithms unless a base is written.

## Claim

Let H = sum_i h_{i,i+1} be a translation-invariant nearest-neighbour chain on (C^d)^L whose deep infrared is a finite set of unitary CFT sectors (c_a, v_a), with Casimir coefficient A = sum_a c_a v_a defined by E_0(L) = e_0 L - pi A/(6L) + o(1/L) on vacuum-compatible rings, and let M be the gauge-optimized bond half-spread. Write rho_2^L for the two-site reduced state of the ring ground state and rho_2^inf for that of the infinite chain, and define the ring-to-line distinguishability coefficient

    D_2 = liminf_{L} L^2 || rho_2^inf - rho_2^L ||_1 .

1. [Proved] pi A/(6 L^2) = Tr[ h~ (rho_2^inf - rho_2^L) ] + o(L^-2) for every gauge representative h~ of the bond, and therefore

       A <= (6/pi) M D_2 .

   M is the dual norm of the constraints that rho_2^inf - rho_2^L automatically satisfies (traceless, equal one-site marginals), so no sharper Hölder step is available from two-site marginals alone.

2. [Proved] The Hölder step is an equality, D_2 = (pi/6) A/M, for every chain whose bond is h = x P_1 + y P_2 with P_1 + P_2 = 1 two orthogonal projectors commuting with the ring ground-state two-site marginal. This covers the Heisenberg chain and every SU(N) swap chain. For these chains the conjecture A <= pi M log_2 d is the same statement as D_2 <= (pi^2/6) log_2 d.

3. [Proved, free-fermion solution] For the XX chain, D_2 = (pi + 2)/3 = 1.7139, which exceeds (pi^2/6) log_2 2 = 1.6449. So the two-site Hölder route cannot reach the sharp constant pi/log 2. It can at best give A <= (6/pi) kappa M log_2 d with kappa >= (pi+2)/3, that is a coefficient of M log_2 d of at least 2(pi+2)/pi = 3.27 in place of the conjectured pi, and only if a universal bound D_2 <= kappa log_2 d holds, which is not proved.

4. [Proved] Any argument whose only lattice input is a comparison of energies or free energies between rings, open chains, or rings of different lengths loses at least 2M per cut, which swamps the O(1/L) Casimir term. This is the reason routes (i), (ii) and (iv) of the task, in their energy form, stall at a preparation or boundary scale, as C69 found for the entanglement-rate route (iii).

No bound of the form A <= C M log d with C independent of a crossover scale is proved here. The statement in item 1 is a reduction, not a proof: it moves the whole difficulty into the single dimensionless number D_2.

## Status

Items 1, 2 and 4 are proved at the level of elementary operator inequalities plus the standard infrared assumption (stated in Section 0). Item 3 is proved given the free-fermion solution of the XX chain and checked numerically to L = 2000. The numerical tables in Section 5 are exact diagonalization of small rings. The reduced conjecture of Section 3 (a universal bound on D_2) is open.

## 0. Setting, conventions, and what the infrared assumption buys

[Assumed] The chain is critical and its finite-size spectrum on vacuum-compatible rings is that of the CFT sectors, so E_0(L) = e_0 L - pi A/(6L) + o(1/L). Lattice spacing and hbar are one. For a local operator O on a fixed number of sites, the ring expectation converges to the infinite-chain one with corrections O(L^-2) up to logarithms (the leading correction comes from the stress tensor one-point function; primaries have vanishing one-point functions on the ring). This is the only place the CFT form is used, and only in the deep infrared, at L -> infinity. No crossover scale is introduced anywhere.

[Definition] M = inf over Hermitian o on C^d and real b of || h + o (x) 1 - 1 (x) o - b 1 ||. A near-minimizing representative is written h~ (centred, spectrum in [-M, M] up to a constant absorbed into e_0). On a ring the gauge terms telescope, so H is unchanged.

[Proved, convention check] Heisenberg h = S.S: eigenvalues 1/4 (triplet) and -3/4 (singlet), half-spread 1/2; `bond_scale_M.py` confirms the gauge does not lower it, M = 1/2, for XXZ at Delta = -0.5, 0, 0.5, 1. SU(N) swap h = P: M = 1 (checked for N = 2, 3). Velocities: Heisenberg v = pi/2, XX v = 1, XXZ at Delta = cos(gamma) v = pi sin(gamma)/(2 gamma), SU(N) swap v = 2 pi/N with c = N - 1.

## 1. Why energy comparisons cannot work (route (i), (ii), (iv) in energy form)

[Proved, Lemma 1 (transfer inequality)] Let Z_ring(L, beta) = Tr exp(-beta H_L) on the ring and f(beta) the infinite-chain free energy density of the centred Hamiltonian. Then for all L and beta

    | log Z_ring(L, beta) + beta L f(beta) | <= 2 beta M .

Proof. Open chains of lengths L and K joined by one centred bond differ in Hamiltonian by an operator of norm at most M, so by eigenvalue perturbation (or Peierls-Bogoliubov) a_L = log Z_open(L, beta) satisfies |a_{L+K} - a_L - a_K| <= beta M. A sequence that is additive up to a constant is within that constant of linear: |a_L - L alpha| <= beta M with alpha = lim a_L/L = -beta f(beta). Closing the ring adds one more bond, |log Z_ring - a_L| <= beta M. Add the two.

[Proved, Corollary] log Z_ring >= -beta E_0(L) gives E_0(L) >= L f(beta) - 2M for every beta, and beta -> infinity gives E_0(L) >= e_0 L - 2M, that is pi A/(6L) <= 2M + o(1/L). This is the trivial statement A <= 12 M L/pi. At the self-dual temperature v beta = L, where a CFT would relate the Casimir channel to the thermal channel (c v to c/v), the lattice inequality still carries the error 2 beta M = 2ML/v, which is larger than the Casimir term beta pi A/(6L) = pi A/(6v) by a factor of order L. The error is the boundary free energy, a crossover quantity. This is the sharp form of the obstruction in C69 Section 1, now for the Euclidean route (ii) as well.

[Proved, route (iv)] The exact infrared relation is E_0(2L) - 2E_0(L) = pi A/(4L) + o(1/L). A variational trial state on the 2L ring built from two copies of the L-ring ground state gives E_0(2L) <= 2E_0(L) + 2 [Tr h~ (rho_1 (x) rho_1) - <h~>_L], an O(M) upper bound, which has the wrong sign to be useful and the wrong order. The useful direction, a lower bound on E_0(2L) or equivalently an upper bound on e_0 by an infinite-chain trial state built from ring data, requires a translation-invariant state on the line whose energy density exceeds <h~>_L by only O(M log d / L^2). Product gluing of ring blocks costs O(M)/L. So route (iv) is the same problem as Section 2 in disguise, and the gluing cost per bond that would be needed is exactly Tr h~ (rho_2^inf - rho_2^L).

[Proved, route (i)] E_0^open(L) - E_0^ring(L) = f_bdy + pi A/(8L) + o(1/L) in the infrared, where f_bdy is the nonuniversal boundary energy, of order M. Any Schmidt-rank or Lieb-Robinson estimate of the cut bond is an estimate of f_bdy, which is not constrained by A. A Casimir route must therefore subtract f_bdy exactly, which again means working with local states at 1/L^2 precision rather than energies at O(1) precision.

## 2. The two-site reduction

[Proved, Proposition 1] For any gauge representative h~,

    e_0 - E_0(L)/L = Tr[ h~ (rho_2^inf - rho_2^L) ] .

Proof. E_0(L)/L = <h_{i,i+1}>_L for every bond by translation invariance, and e_0 = lim of the same for the line. The gauge term contributes Tr[o (rho_1^inf - rho_1^L)] - Tr[o (rho_1^inf - rho_1^L)] = 0 because both partial traces of a two-site marginal of a translation-invariant state are the same one-site state. The constant b cancels. Inserting the infrared form of E_0(L) gives pi A/(6L^2) on the left up to o(L^-2).

[Proved, Proposition 2] A <= (6/pi) M D_2. Proof: Hölder, |Tr h~ X| <= ||h~|| ||X||_1 with ||h~|| arbitrarily close to M, then multiply by L^2 and take liminf.

[Proved, Proposition 3 (M is the dual norm)] Let V be the real subspace of two-site Hermitian operators spanned by 1 and by o (x) 1 - 1 (x) o for Hermitian o. Its annihilator under the trace pairing is W = {X Hermitian : Tr X = 0, Tr_1 X = Tr_2 X}. By the duality between the operator norm and the trace norm restricted to a subspace (a finite-dimensional Hahn-Banach argument, equivalently strong duality for this semidefinite program),

    M = min_{g in V} || h - g || = max { Tr h X : X in W, ||X||_1 <= 1 } .

The difference rho_2^inf - rho_2^L lies in W. So Proposition 2 is the best inequality obtainable from the two facts "X is in W" and "||X||_1 is known", and the gauge optimization in the definition of M is not a choice but the exact dual of the constraints Proposition 1 uses. Any improvement must use more information about rho_2^inf - rho_2^L than membership in W.

[Proved, Proposition 4 (saturation class)] Suppose h = x P_1 + y P_2 with P_1 + P_2 = 1 projectors and suppose both rho_2^L and rho_2^inf commute with P_1 (for instance by a symmetry under which the bond space splits into exactly two irreducible pieces, as for SU(2) with d = 2 and for SU(N) with h = P, where the pieces are the symmetric and antisymmetric squares). Then rho_2^inf - rho_2^L = delta (P_1/r_1 - P_2/r_2) with r_i = rank P_i, its trace norm is 2|delta|, Tr h~ (rho_2^inf - rho_2^L) = (x - y) delta, and M = |x - y|/2 when the gauge does not lower the spread (checked above for these chains). Hence Tr h~ X = M ||X||_1 exactly and D_2 = (pi/6) A/M. For the Heisenberg chain this gives D_2 = pi^2/6 exactly, equal to (pi^2/6) log_2 d. For SU(N) swap chains D_2 = pi^2 (N-1)/(3N), which is below (pi^2/6) log_2 N for N >= 3 and equal at N = 2.

[Heuristic, interpretation] In the infrared, rho_2^inf - rho_2^L = (pi/(6L^2)) sum_a c_a v_a tau_a + o(L^-2) with tau_a the two-site operator dual to the sector-a stress tensor, normalized so that Tr h~ tau_a = 1. The velocity appears because the ring stress-tensor expectation is c v pi^2/(6 L^2) in lattice energy units. Slow sectors are therefore automatically suppressed in D_2, which is the feature any valid proof must have (c alone is not bounded by a function of d). The conjecture in this language reads ||sum_a c_a v_a tau_a||_1 <= pi log_2 d, while Proposition 3 gives only ||tau_a||_1 >= 1/M.

## 3. The reduced conjecture and its test on the XX chain

[Reduced conjecture, open] There is a universal kappa with D_2 <= kappa log_2 d for every chain in the class. By Proposition 2 this would give A <= (6/pi) kappa M log_2 d.

[Proved, free-fermion solution] XX chain, h = S^x S^x + S^y S^y, which is hopping with t = 1/2, e_0 = -1/pi, v = 1, c = 1, A = 1, M = 1/2. The two-site marginal is a two-mode Gaussian state with <c_i^dagger c_i> = 1/2 and g = <c_1^dagger c_2> real. In the occupation basis it is block diagonal: entries 1/4 - g^2 in the empty and doubly occupied sectors and the block (1/4 + g^2) 1 + g sigma^x in the singly occupied sector. On the ring g_L = -E_0(L)/L, on the line g_inf = 1/pi, so delta g = g_inf - g_L = -pi A/(6L^2) = -pi/(6L^2). Then rho_2^inf - rho_2^L has entries +1/(3L^2) in the two even sectors and the block -(1/(3L^2)) 1 + delta g sigma^x, with eigenvalues (-1/3 -/+ pi/6)/L^2. Its trace norm is (2/3 + pi/3)/L^2, so

    D_2^XX = (pi + 2)/3 = 1.713865 .

Of this, pi/3 = (pi/6) A/M is the part seen by the bond energy and 2/3 is a density-density part invisible to h~ (it lives in the sectors where h~ vanishes). `xx_trace_norm.py` evaluates the exact Gaussian marginals and confirms L^2 ||rho_2^inf - rho_2^L||_1 = 1.713865 at L = 2000 with the L^2 Tr h~ (.) 6/pi column equal to A = 1 to 1e-6.

[Proved, consequence] Since (pi^2/6) log_2 2 = 1.644934 < 1.713865, the sharp constant kappa = pi^2/6 of the reduced conjecture is false, even though the XX chain satisfies A = 1 < pi/2 = pi M log_2 d comfortably. The two-site Hölder route therefore proves at best A <= (6/pi) kappa M log_2 d with kappa >= (pi+2)/3, that is a coefficient of M log_2 d of at least 2(pi+2)/pi = 3.27 in place of the conjectured pi = 3.14. The sharp constant pi/log 2 needs information beyond two-site marginals and the trace norm.

[Numerical, XXZ at Delta = 0.5] Exact diagonalization to L = 18 gives A = 1.2994 against the Bethe value v = 3 sqrt(3)/4 = 1.2990, and D_2(L) decreasing through 1.669 at L = 18 with pairwise estimates 1.680, so D_2 appears to lie between the Heisenberg value 1.645 and the XX value 1.714. The strict inequality Tr h~ X < M ||X||_1 holds at Delta = 0.5: in the units of Proposition 2, A = 1.30 against (6/pi) M D_2 = 1.59 at L = 18.

## 4. Assessment of the four directions

(i) Casimir as a ring-minus-cut energy difference with Schmidt-rank or Lieb-Robinson control at the cut. [Proved obstruction, Section 1] The cut costs the boundary energy, of order M, which is not fixed by A. Lieb-Robinson controls how fast the cut is felt, not how much energy it costs. This route cannot produce a 1/L term without an exact subtraction, which turns it into Section 2.

(ii) Modular or Euclidean relation to a thermal quantity at temperature set by L. [Proved obstruction, Lemma 1] The lattice version of torus modular invariance is the transfer inequality, whose error 2 beta M is the open-chain boundary free energy. At v beta = L the error exceeds the Casimir term by a factor of order L. A lattice-rigorous version would need the boundary free energy to O(M log d / L), which is a crossover statement. The modular-flow variant (first law of entanglement for an interval) is velocity-free at the level of entropies and so cannot bound A; velocities enter only through the energy normalization of the modular Hamiltonian, where the two-site cutoff region is outside CFT control.

(iii) Entanglement production of a single bond with a scale-free protocol. [Assessed, no new theorem] C69 proved that the rate measures A/ell and that the small-incremental-entangling bound scales the same way, so the preparation scale survives. The only scale-free combination of ring data at order 1/L^2 found here is Proposition 1, which is an energy statement about two-site marginals rather than an entropy statement. I did not find a protocol in which a length cancels against a rate.

(iv) Variational comparison of E_0(L) and E_0(2L). [Proved reduction, Section 1] Equivalent to constructing a line state from ring data with gluing cost Tr h~ (rho_2^inf - rho_2^L) per bond, which is Proposition 1 again. The direction that would bound A needs a lower bound on ring energies, and lower bounds via local decompositions (Anderson) lose O(M) per cut.

The common conclusion is that the Casimir coefficient is a statement about local states at 1/L^2 precision, and that the lattice quantity dual to the gauge-optimized M is the two-site marginal difference. The bound A <= (6/pi) M D_2 is the end of what energy-norm reasoning can give.

## 5. Checks

[Numerical] `ed_rings.py` (output `ed_rings.out`): sparse Lanczos on rings, Casimir fit of E_0(L)/L = e_0 - a/L^2 + b/L^4 over the four largest sizes, two-site marginals, and D_2(L) with rho_2^inf estimated by an elementwise fit of the same form (model-free pairwise estimates ||rho_L - rho_L'||_1/(L^-2 - L'^-2) are also printed).

| chain | d | M | L used | A (fit) | A (exact) | D_2(L) at largest L | (pi/6) A/M | (pi^2/6) log_2 d |
|---|---|---|---|---|---|---|---|---|
| XX (free fermion, Gaussian) | 2 | 1/2 | 2000 | 1.000000 | 1 | 1.713865 | 1.047198 | 1.644934 |
| XXZ Delta = 0.5 | 2 | 1/2 | 8 to 18 | 1.2994 | 1.2990 | 1.669 (decreasing) | 1.360 | 1.645 |
| Heisenberg | 2 | 1/2 | 8 to 18 | 1.578 (log corrections) | 1.5708 | 1.661 (decreasing toward 1.645) | 1.645 | 1.645 |
| SU(3) swap | 3 | 1 | 6, 9, 12 | 4.215 | 4.189 | 2.235 (decreasing toward 2.193) | 2.193 | 2.607 |

The Heisenberg and SU(3) rows have D_2 approaching (pi/6) A/M from above, consistent with Proposition 4 (equality in the limit). The XX row is the exact counterexample to kappa = pi^2/6. Convergence of D_2(L) is from above in every case at the sizes reached.

[Numerical] `bond_scale_M.py`: convex minimization of the half-spread over the gauge subspace, several random starts, for XXZ and the swap chains. The gauge does not lower M for any of them.

## 6. What remains

The precise gap is a universal upper bound on D_2, or equivalently on the trace norm of the two-site ring-to-line marginal difference at order 1/L^2, in terms of log d. Facts established about it: D_2 >= (pi/6) A/M always, with equality on the two-projector class; D_2 = pi^2/6 for Heisenberg; D_2 = (pi+2)/3 for XX; D_2 = pi^2 (N-1)/(3N) for SU(N) swap. A proof of D_2 <= kappa log_2 d for some kappa would give A <= (6/pi) kappa M log_2 d, with kappa >= (pi+2)/3 forced by XX. Reaching the sharp pi/log 2 requires using the structure of rho_2^inf - rho_2^L beyond its trace norm, for example the fact that it is the two-site marginal of a specific operator (the difference of two ground-state projectors restricted to a translation-invariant class), or an inequality that pairs h~ with the marginal difference through an entropic quantity rather than a norm. I did not find such an inequality. A candidate test for any proposed kappa is the SU(8) symmetric-rank-2 chain of C14, where A/(2 pi M) = 1.05 and d = 36.

## 7. Sources used

`../../specs/_bound_background.md`, `../../codex/work_C69/RESULT.md`, `../refined_conjecture/README.md`. The parallel runs C98 to C101 had no PROGRESS.md at the time of writing (their folders were empty), so nothing from them is used or duplicated. The small-incremental-entangling theorem and the Calabrese-Cardy quench results are cited through C69 and are not reused here.
