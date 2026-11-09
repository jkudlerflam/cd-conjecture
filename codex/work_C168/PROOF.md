This work was done by Claude Opus 5.5 (agent C168, 2026-11-08), not by Codex. Astra was not used (out of credits).

# C168 derivations: the symmetry-free stratum N

Tags: [proved] (complete argument, or an identity checked to machine precision), [derived] (physics-standard derivation with named inputs), [numerical] (script and output named), [heuristic], [assumed]. A tag applies to the paragraph or displayed claim it opens. Conventions are those of C154: H = sum_j [sigma_j^T J sigma_{j+1} + B . sigma_j], N is the 8-dimensional set of bonds with no on-site symmetry, translation is exact. c is defined per clustering ground state (C174 R8), so a cat or a first-order coexistence never counts as c = 2.

## 0. Referee summary of the key step

C174 R1 shows that no connectivity argument can give (H), because upper bounds on c travel only from higher c to lower c. So this run attacks the dynamics directly. It classifies the frontier of a would-be c > 1 component of N into the softening events that create it, and closes those events that have a small parameter.

The new closed event is the generic double onset (Lemma D). Take a point of N where two independent branches of a gapped phase soften together. This has codimension 2 in N, because nothing relates the two branches. Near such a point both species are dilute, and in one dimension dilute particles with any short-range interaction become impenetrable. Their cyclic order along the chain is then fixed up to exchange processes. In N the two effective masses differ generically, with relative imbalance eta, and that makes the energy depend on the order. The cost is 0.31 eta^2 hbar^2 n^2/M per unlike interface (numerical, imp_mix.out, with a second-order argument for the sign). Passing requires low-energy transmission, which vanishes like k a_t in one dimension, so the exchange is of order hbar^2 n^3 a_t/M (imp_swap2.out). The pseudospin chain of the two species therefore has |Delta_ps| of order eta^2/(n a_t), which diverges at the onset. Either the species segregate, so a single species wins and the transition is first order, or they alternate, so the pseudospin is gapped. In both cases c = 1, and no c = 2 region touches the double onset. Equal masses, a zero-energy transmission resonance, or a bound state between the species each cost one more tuning.

Together with a codimension classification of onsets out of arbitrary gapped references (Lemma O, which replaces C154 Lemma 3 in N, as C174 R9 asks) and the clock-model bound for tongue melting (Lemma P, the pattern P4 of C174 R5 at weak locking), this leaves a residual (R_N) of five named events. None of them has a small parameter.

## 1. Scope and what has to be excluded

[derived, C154 4.2 and 5.2, C174 R5] A c > 1 conformal point of N is one of the following.
- P1, Ising × floating boson with sigma at momentum pi (codimension 1). It needs an adjacent "floating + SSB" phase (FS).
- P3, two floating bosons (codimension 0), an open region.
- P2 and the other reducible theories (codimension at least 2), for example Ising × boson with uniform sigma.
- P4, irreducible c > 1 theories (codimension at least 2): Z_m parafermions with m ≥ 5 (c = 8/7 at m = 5), N = 1 superconformal minimal models with 1 < c < 3/2, W_3 minimal models (6/5, 10/7, ...), and SU(2)_3 off the cyclic stratum.

[proved, C174 Lemma R1, given (A2) and gap stability] The closure of a c > 1 component contains no gapped bond and no conformal bond of smaller c. So (H) cannot be obtained from connectivity or deformation to the U(1) class. Lead 3 is dropped in that form. What survives of it is Proposition F (section 6), which only sorts the frontier into dynamical events.

K_II is outside the scope of this run. It carries an antiunitary that maps k to -k, so its two-valley onsets have equal masses by symmetry and Lemma D does not apply there (that is C162's case). As C174 R6 asks, K_II is left under (A6)/(H*).

## 2. Lead 1: there is no kinematic route at one qubit per site

Lemma 1. (a) [proved, standard] On N the only symmetry is translation Z. Its anomaly group is H^3(Z, U(1)) = 0, and a filling (Oshikawa or LSM-type) constraint needs a microscopic U(1). So there is no filling constraint on the emergent charges of a floating liquid in N.

(b) [proved for H0; derived for H0 + W; numerical support] The range-3 qubit chain

    H0 = sum_j [ -t1 e^{i a1} s+_j s-_{j+1} - t2 e^{i a2} s+_j (-Z_{j+1}) s-_{j+2} + h.c. ] - (mu/2) sum_j Z_j

is Jordan-Wigner free fermions with band eps(k) = -2 t1 cos(k - a1) - 2 t2 cos(2k - a2) - mu. On an open set of (t2, a1, a2, mu) the band has two Fermi seas with four Fermi momenta that are rationally independent of 2 pi at generic points, so c = 2 exactly.

Now add any small nearest-neighbour two-site term W, which removes every on-site symmetry. Its fermion-parity-even part is local in the fermions. Its number-changing pieces carry nonzero integer combinations of the Fermi momenta, which do not vanish off a union of resonance hypersurfaces, so they are forbidden by translation as in C154 4.1. Its parity-odd part (single X, Y and XZ-type terms) is a vertex operator with a Jordan-Wigner string in the two-component bosonization, and it also carries an incommensurate momentum. So c = 2 survives on an open set of a symmetry-free qubit chain with range 3.

[numerical, calib_c2.py, calib_c2_local.out] At rings 12/16 the median c_S is 1.99 over 59 two-sea samples (W of norm 0.1), against 1.00 over 208 one-sea samples. For Ising × boson, the range-3 and range-4 qubit models with c = 3/2 listed in C156 (Rahmani-Zhu-Franz-Affleck, cluster-Ising) play the same role.

(c) [derived, consequence] Both suggested routes of lead 1 use only the qubit dimension and translation, and (b) shows these are compatible with two floating bosons. Those routes are a filling constraint on each emergent charge, and compatibility of two vertex lattices with one spin per site. So any exclusion of P3 or P1 on B_2 must use the range-2 structure of the bond, and it must be dynamical.

## 3. Lemma O: onsets out of gapped phases in N, by codimension

This replaces C154 Lemma 3 in N, which C174 R9 shows is controlled only at exact saturation manifolds and near h_F. It also covers periodic references (C174 R8).

Lemma O. [derived; inputs: transversality, and the absence of any symmetry relating momenta or branches in N] Let G be a gapped phase of N, with a unique vacuum or p vacua related by broken translation, whose gap closes continuously at h_0, at finitely many momenta. Write the low excitation branches of a vacuum as E_a(k; h), with k in the zone of that vacuum. The gap is min_{a,k} E_a. Since N has no antiunitary and no reflection, nothing forces E_a(k) = E_a(-k). No unitary relates distinct branches, so degenerate minima are non-generic.
- At a generic point of the onset set (codimension 1) the minimum is attained by one branch at one momentum k_1. If k_1 is not 0 or pi, the soft mode is complex: an emergent U(1), a Pokrovsky-Talapov or KT onset into one floating boson, c = 1. If k_1 is 0 or pi, k_1 = -k_1 and the mode is self-conjugate: one real field, Ising, c = 1/2.
- At codimension 2 the possibilities are:
  - (O1) two incommensurate minima, from two branches or from one branch, which is Lemma D;
  - (O2) a real mode at 0 or pi together with an incommensurate one, the "P1 corner";
  - (O3) real modes at 0 and pi, two Majoranas, c ≤ 1 (Ashkin-Teller or Gaussian, C68);
  - (O4) a Lifshitz point of one branch (quartic minimum, or a minimum leaving 0 or pi), one species, non-conformal.
- Codimension 3 and higher adds equal masses or transmission resonances in (O1), and triple onsets.

Proof sketch. min_{a,k} E_a = 0 is one real equation. Requiring it at two distinct (a, k) is two equations, independent in the absence of a symmetry that maps one solution to the other. The real or complex character follows from k = -k mod 2 pi. QED.

Two refinements. [derived]
- Momenta 2 pi/p with p = 3 or 4. A complex mode exactly at k_1 = 2 pi/p carries an allowed clock term sigma^p at zero momentum, so it gives Z_p clock physics. Reaching k_1 = 2 pi/p is one extra condition, unless the momentum is locked inside a commensurate tongue, which Lemma P covers.
- Charge-conjugation actions. A boson on which translation acts by charge conjugation times a shift (the "C row" of C171) arises at an onset only from the (O3) pair. Write psi = chi_0 + i chi_pi. Translation fixes chi_0 and flips chi_pi, so it maps psi to psi^dagger. That gives one boson (c = 1). A boson pair that translation swaps is the same structure for the relative boson. So both enter only through a second, independent sector, which is (O2) or a softening inside a liquid.

[proved, special case] For a product reference this reproduces C154 Lemma 3. Its quadratic g(x) cannot vanish at an endpoint and at an interior point, so (O1) and (O2) are absent there.

## 4. Lemma D: generic double onsets are immiscible or alternating

Setting. h_0 is a type-(O1) point: two species A and B with band minima k_A ≠ k_B, masses M_A ≠ M_B, and relative imbalance eta = (M_B - M_A)/(M_A + M_B) ≠ 0, over a gapped vacuum with correlation length xi_0. The deficits mu_A and mu_B are two independent local coordinates transverse to the codimension-2 set. Hypotheses: (h1) at h_0 there is no A-A, B-B or A-B bound state, so the lowest composites are A and B themselves; (h2) the A-B transmission has no zero-energy resonance.

Lemma D. [derived; inputs: dilute-gas universality in one dimension, Lemma O, and the numerical checks below] There is a ball U around h_0 in which no point has c = 2. In the quadrant mu_A, mu_B > 0 of U the ground state is a single-species liquid (c = 1), with a first-order line between the A and B liquids, or an alternating ABAB crystal with gapped pseudospin and c = 1.

Steps.
1. Emergent number conservation. A conversion A <-> B carries k_A - k_B ≠ 0, and pair creation AA or BB carries 2 k_{A,B} into a branch whose minimum is not at -k. So these are off shell at low energy, and they only renormalize the couplings, as in C154 4.1 and C162. N_A and N_B are conserved in the infrared.
2. Impenetrability. For densities n_A, n_B ≪ 1/xi_0, the gas is a two-species nonrelativistic gas with interactions of range xi_0. In one dimension the relative wave function at low momentum k has transmission amplitude t(k) = O(k a_t), apart from a zero-energy resonance, which (h2) excludes. So the species are impenetrable at leading order, and the cyclic order of A's and B's along the chain is conserved. Bound states are excluded by (h1).
3. Order energy from mass imbalance. [numerical, imp_mix.py, imp_mix.out] For impenetrable hard-core particles with hoppings t_A = 1 + eta and t_B = 1 - eta on a ring, E(ABAB) - E(AABB) = 19.6 eta^2/L^2 for L = 16 to 32. The coefficient divided by eta^2 is 19.75, 19.65, 19.62 at eta = 0.1 and L = 16, 24, 32, and its eta-dependence is weak (20.15 at eta = 0.3). The energy is exactly zero at eta = 0, which is the Girardeau degeneracy. This is 0.31 eta^2 hbar^2 n^2/M per unlike interface with n = 4/L and hbar^2/(2M) = t, in agreement with C164's 0.3 from six particles.

   [derived] The sign and the eta^2 scaling follow from second-order perturbation in the kinetic imbalance. At first order every particle carries the same kinetic share on a ring, independent of the order. At second order segregated domains can trade density between heavy and light regions, which alternating orders cannot. Scattering-length differences contribute at O(hbar^2 n^3 Delta a/M) per interface.
4. Exchange. [derived, numerical check imp_swap2.out] Passing costs the transmission probability at relative momentum k of order n, so J_perp = O(hbar^2 n^3 a_t/M), as in the Gaudin-Yang strong-coupling expansion. On the lattice with a swap amplitude s = 0.2 and 2 + 2 particles, the exchange energy satisfies dE_x L^3 = -104, -112, -116, -118, -119, -120 for L = 12, 16, 20, 24, 28, 32. That approaches a constant, which is the n^3 law.
5. Pseudospin chain. On the squeezed lattice the species form an XXZ chain with J_z = kappa eta^2 hbar^2 n^2/M + O(n^3) and J_perp = O(hbar^2 n^3 a_t/M), together with a pseudospin field from mu_A - mu_B and from the energy difference. Its pseudospin U(1) is the emergent one of step 1. So Delta_ps = J_z/|J_perp| ≥ C eta^2/(n a_t), which diverges as n → 0.
   - For J_z > 0 (the mass-imbalance sign): an Ising ferromagnet in a field, fully polarized, so one species wins. Coexistence happens only on the line P_A(mu_A) = P_B(mu_B), where the transition is first order.
   - For J_z < 0 (possible only if scattering terms dominate, which they do not as n → 0): pseudospin Néel, an alternating crystal with gapped pseudospin and c = 1.
   - c = 2 requires -1 < Delta_ps ≤ 1. That is excluded for n < n_* of order eta^2/a_t, and so on the ball U where both densities stay below n_*.

   [numerical check, imp_swap2.out] The ground-state weight on alternating orders at eta = 0.5, s = 0.2 is 0.327, 0.283, 0.242, 0.206, 0.175, 0.149 for L = 12 to 32. It falls steadily, which is segregation. At eta = 0 it stays at 0.514 to 0.517, an L-independent mixed state. With 4 particles this is a trend, not a phase diagram.

Consistency checks.
- [derived] g_AB → 0 (free fermions with two valleys): a_t → ∞ and n_* → 0, so the c = 2 two-sea phase of free fermions is recovered continuously. Lemma D says only that any interaction opens an immiscible sliver at the double onset, of density width of order eta^2 M g_AB.
- [derived] eta = 0 by symmetry: inversion in the J1-J2 chain near saturation, or the antiunitary in K (C162). Then J_z comes from scattering lengths alone, and c = 2 near the onset is possible. It occurs in J1-J2 (C156) and was argued against in K by C162's rearrangement lemma. This is why Lemma D needs N.
- [derived] With an A-B bound state, the molecule's onset at E_A + E_B = E_bind comes first (codimension 1). Lemma O gives one species there, and the remaining single-particle onsets are guest onsets inside the molecule liquid, which is the residual R(i).

Consequence. [derived] In N a c = 2 region can touch an onset of a gapped phase only at (a) codimension-3 points with eta = 0, or with a transmission resonance, or with bound-state degeneracies, inside wedges |eta| ≲ (n a_t)^{1/2} and only if the scattering-length pseudospin is easy-plane there; or (b) through a guest onset inside an already dense liquid.

## 5. Lemma P: P4 at weak locking

C154 2.3 and C163 Proposition P predict, near the U(1) class, a staircase of commensurate tongues where a vertex cos(m Phi) locks the floating boson. C174 R5 asks whether a tongue can melt through an irreducible c > 1 point.

Lemma P. [derived; inputs: Coulomb-gas analysis of clock models (José, Kadanoff, Kirkpatrick and Nelson 1977; Elitzur, Pearson and Shigemitsu 1979), self-dual sine-Gordon criticality (Lecheminant, Gogolin and Nersesyan, Nucl. Phys. B 639 (2002) 502), and chiral clock melting (Huse and Fisher 1984)] Suppose that near a tongue the infrared is one compact boson perturbed by vertex operators with near-Gaussian dimensions. In C159's convention, x(V_{p,n}) = p^2/(4K) + n^2 K. Then every transition out of the tongue has c ≤ 1 or is non-conformal, for every m.

Steps.
- The locking vertex V_{m,0} has dimension m^2/(4K). The minimal dual vertex V_{0,1} has dimension K. The product is m^2/4.
- Both are relevant (dimension < 2) only if m^2 < 16, that is m ≤ 3. At m = 4 both are marginal together.
- For m ≥ 5 there is a finite window K in (2, m^2/8) where neither is relevant. So a c = 1 phase separates the locked tongue from any phase disordered by V_{0,1}, and the transitions are KT (c = 1).
- With the chiral drift that N allows (spin-±1 operators, an incommensurate mismatch delta), the transition into the tongue is Pokrovsky-Talapov (z = 2). For m ≥ 5 the chiral clock melts through a floating phase.
- For m ≤ 4 the self-dual points are Ising (m = 2, c = 1/2), the Z_3 parafermion, which is three-state Potts (c = 4/5), and the Z_4 parafermion (c = 1), or chiral non-conformal transitions.
- In N, V_{0,1} itself carries momentum 2 k_F and is usually forbidden, which only removes the gapped disordered side.

So the Z_m parafermion with m ≥ 5 (c = 2(m-1)/(m+2) > 1), and more generally any irreducible c > 1 fixed point, is not reached while the vertex dimensions stay near their Gaussian values. Reaching it needs clock amplitudes large enough to move the dimensions far from those values, as at the Fateev-Zamolodchikov points.

[derived, count] At such a fixed point N allows the translation-invariant thermal operators: for Z_5, epsilon_1 of dimension 4/7 and epsilon_2 of dimension 12/7. It also allows chiral spin-±1 operators. So P4 has codimension at least 2 (C154 table: at least 2 for Z_5), and the same holds for the N = 1 and W_3 minimal models.

Check. [proved, arithmetic] m = 2, 3, 4 give c = 1/2, 4/5, 1, all ≤ 1, so P4 can only start at m = 5, where the window K in (2, 25/8) is nonempty.

## 5b. The candidates added by C171

[derived; input: the relevant chiral perturbation of three-state Potts, the spin-1 field of dimensions (7/5, 2/5) and total dimension 9/5 (Huse and Fisher 1984; Cardy 1993)] Potts × floating boson, with Z_3 order at momentum 2 pi/3 (c = 9/5), is a new row for the N table.
- The thermal epsilon is allowed. So is the chiral operator, which is translation invariant, spin 1 and relevant, because N has no reflection or antiunitary to forbid it. sigma V and epsilon V carry incommensurate momenta.
- So the conformal Potts × boson fixed point has k ≥ 2 in N, not 1. At generic points of the would-be one-tuning surface the transition is chiral: Huse-Fisher, believed non-conformal, or split by an intermediate floating phase. Both sides then have c ≤ 1 + 1 with a non-conformal or c = 1 transition in between.
- The conformal c = 9/5 point sits on a codimension-2 locus and belongs to the reducible codimension-2 class with P2.

[derived] The other C171 items are covered by mechanism rather than by table.
- The C row and swapped boson pairs come only from (O3) pairs (Lemma O refinement).
- Two-boson phases whose second boson is a composite with wavevector in Z theta + Z pi (where C163 Lemma L fails) are, in N, either a bound state at a double onset, excluded by hypothesis (h1) of Lemma D and sent back to Lemma O, or a guest inside a dense liquid, which is (R_N)(i). This run uses no Lemma-L-type locking.
- Codimension-2 candidates are generic and need not lie near symmetric strata (C171 R2). Proposition F does not use proximity to strata. P2 has only c = 1 neighbours and enters through (F2), as a uniform Ising mode softening inside a liquid.

## 6. Lemma H and Proposition F

Lemma H (the AF Heisenberg neighbourhood in N). [derived; inputs: SU(2)_1 with currents and marginal current-current couplings, Garate and Affleck 2010] A small N-perturbation of 4 S.S gives SU(2)_1 plus the following terms.
- tr g and N^a are odd under the centre (momentum pi), so exact translation forbids them to all orders.
- The chiral fields h_L.J_L + h_R.J_R with h_{L,R} = c1 B ± c2 D.
- The marginal couplings g_ab J_L^a J_R^b.

For |h_L| ≠ |h_R|, both nonzero (generic in N), rotate the left movers so that h_L is along z_L and the right movers so that h_R is along z_R. Then phi_L and phi_R shift by different incommensurate wavevectors. Every marginal term except J^z_L J^z_R oscillates, and it does so on the length v/|h| of order 1/delta, long before the marginal flow acts (on a scale of order e^{1/delta}). The infrared is one floating boson, c = 1. On the special loci (|h_L| = |h_R|, for instance D perpendicular to B to leading order, or h_L = 0) more marginal terms survive and give XXZ-type physics with c in {0, 1}. Resonances that would re-admit N^a or tr g need h of order the bandwidth.

So c ≤ 1 on a full punctured neighbourhood at the marginal level. This is the uniform-in-direction form of C159 3.2 (zero field) and C165 4 (hidden Heisenberg), and it partly answers C174 R3 there.

Proposition F (frontier classification). [derived, given (A1) with c per sector, (A2), C174 R1, Lemmas O, D, P, H and the controlled set S_ctrl of C174 R2] Let Sigma be a c > 1 component of N, and let Sigma^ be the union of all c > 1 components joined to it through conformal frontier points. Every point of the outer frontier of Sigma^ is non-conformal and is one of the following.
- (F1) A first-order point.
- (F2) A point where, coming from a c = 1 liquid (with any translation action: floating, half-floating, commensurate, or charge conjugation), a sharp excitation independent of the liquid's boson reaches zero energy. This is a guest onset (P3, composites included), an Ising mode at pi (P1 from the floating side), a uniform Ising mode (P2), or a Z_3 mode at 2 pi/3 (Potts × boson).
- (F3) A type-(O2) corner, or more generally a c ≤ 1 transition of a gapped phase (Ising, two-Majorana Gaussian, Potts or chiral clock) meeting the onset of an independent branch.
- (F4) A codimension-3 exception of Lemma D.
- (F5) A point of S_ctrl.
- (F6) The end of a P4 locus at strong locking.

Proof. Points with a gapped side and a single soft species have c ≤ 1 on the other side (Lemma O), which rules them out. Generic double onsets are excluded by Lemma D. The remaining ways a sector can become gapless next to a c ≤ 1 region are listed. Tongue boundaries at weak locking have c ≤ 1 by Lemma P. QED.

By C174 R1 a component touching only (F5) satisfies (H*). So with (A5) restricted to S_ctrl, c ≤ 1 on N follows from the residual (R_N) below.

## 7. The residual (R_N)

(R_N) [assumed], on N outside S_ctrl and the neighbourhoods of Lemmas D, H and P:
- (i) Guest, Ising or clock mode inside a liquid. In a c = 1 liquid of N (any translation action), no sharp excitation independent of the liquid's boson reaches zero energy while the liquid survives. This applies at momentum 0, pi, 2 pi/3 or an incommensurate k, and includes composites of the host's quanta. It is A_S of C165 and the two halves of A_W (C163) restated for N.
- (ii) P1 corner. No codimension-2 point where a c ≤ 1 transition of a gapped phase (Ising at 0 or pi, two-Majorana, Potts or chiral clock) coincides with the onset of an independent branch of the same phase.
- (iii) Lemma D exceptions. No equal-mass, resonant or bound-state-degenerate double onset with an easy-plane pseudospin at the scattering-length level (codimension at least 3).
- (iv) P4. No irreducible c > 1 fixed point at clock amplitudes of order the exchange.
- (v) First-order islands. No c > 1 component whose outer frontier is entirely first order (F1).

Where it bites, and why the method stops.
- Items (i) and (ii) need a lattice-scale branch to come down at harmonics of order the exchange. In the U(1) class every excitation below the lattice scale belongs to the one boson, by Bethe-ansatz completeness (C163 Fact 2 and Fact 3). Near h_F there is one phason (C160 M1 to M3). Near the AF and hidden Heisenberg points the only sector is SU(2)_1 (Lemma H, C165). So no small parameter is left.
- The dilute analysis of Lemma D does not extend to a guest inside a dense host. There the guest is transparent to the host, because host backscattering costs the guest a momentum 2 k_F that it cannot absorb on shell. So the dilute guest is a stable Tonks gas, and a two-component Luttinger liquid is consistent.
- I checked and rejected a Wentzel-Bardeen-type instability argument for this case. The exact compressibility matrix stays positive semidefinite, as the Bethe-ansatz Hubbard chain near saturation shows.
- At the P1 corner the dilute particle couples to the critical Ising sector only through epsilon, an exactly marginal line defect, so the P1 surface t + g n = 0 emanates without obstruction.
- Items (iii) and (iv) are codimension at least 3 and at least 2, and both need coupling-sized amplitudes.
- Item (v) is invisible to any local softening analysis and is addressed only by sampling (section 8).

(R_N) is implied by (A6) and is strictly narrower. The following parts of (A6) are now derived:
- generic double onsets out of any gapped phase, periodic ones included (Lemma D);
- onsets at codimension 1 (Lemma O);
- the full neighbourhood of the AF Heisenberg point in N (Lemma H);
- tongue melting at weak locking for every m (Lemma P).

## 8. Lead 1, second part: what the uniform sample bounds

[numerical, stageA_analyze.py] In C155 stage A (12,000 uniform points on S^8, rings 12/16), 636 points have crit < 0.15, and 28 of those have c_obj > 1.1. C155 stage D never refined these 28, because refine.py collects only search and gap-closing outputs. So the statement that every rising reading traces to the XZ Néel cat did not cover stage A.

Proxy calibration. [numerical, job 17602943 stage P, typhon_out/summary.txt] The proxies are 1600 range-3 qubit chains with known c (section 2): free fermions with a symmetry-breaking W of norm 0.1 ("free"), and the same with V Z Z, V uniform in [-0.6, 0.6], and W of norm 0.25 ("int"). Fraction flagged:

| family, c | N | F1: c_S(12,16) > 1.5 | F1 and F2: also c_S(16,20) > 1.5 | C155 rule (crit < 0.15 and c_obj > 1.1) |
| :-- | --: | --: | --: | --: |
| free, c = 1 | 513 | 0.170 | 0.016 | 0.008 |
| free, c = 2 | 167 | 0.641 | 0.341 | 0.018 |
| int, c = 1 | 521 | 0.107 | 0.012 | 0.006 |
| int, c = 2 | 163 | 0.331 | 0.178 | 0.006 |

The C155 rule detects about 1 percent of c = 2 chains. The reason is that crit < 0.15 fails on small incommensurate rings for c = 1 and c = 2 alike. So the stage-A statement "no c > 1 among 12,000 points" carried almost no weight against an open c = 2 region (a 95 percent bound near 3/(0.006 × 12000) = 4e-2).

Refinement. [numerical, job 17602943 stages S1 and S2, typhon_out/out/]
- All 394 stage-A points with c_S(12,16) > 1.5, together with the 28 flags (416 in total), were run at rings 16/20. Nine have c_S(16,20) > 1.5.
- The 45 points with c_S(16,20) > 1.5, or c_obj(16,20) > 1.1, or a C155 flag, were run by open-chain DMRG at L = 32, 64, 128 (chi up to 128, truncation below 3e-10). 43 read c_prof(128) ≤ 1.16, most between 0.8 and 1.05.
- Two meet the first stage of the pre-fixed rule: c_prof(128) and the two-size slope 6 dS_mid/d ln L from 64 to 128 both above 1.2, with truncation below 1e-7.
  - idx 10193: c_prof 0.51, 0.61, 1.94 at L = 32, 64, 128. Its L = 128 energy lies 0.06 above the linear continuation of L = 32 and 64. Its 64 -> 128 energy slope is 9e-4 above the 32 -> 64 slope, the largest such anomaly in the set. That points to an unconverged L = 128 state.
  - idx 450: c_prof 0.09, 1.62, 1.21, falling from L = 64. Its energy slopes agree to 3e-5.

The second stage of the rule (FES at chi up to 256 in two- and four-site cells, and DMRG reseeds at L = 128 and L = 192) is job 17603367. Until it reports, these two points are open candidates.

The local refinement (refineA_table.out, 18 of the 28 flags at L = 32 and 64 before I stopped it to reduce load) agrees with the Typhon L = 64 energies to 1e-6.

Bound. [numerical calibration, statistics] Suppose the c > 1 set has normalized volume p on S^8 and a point in it is flagged by F1 and F2 with probability eta (the c = 2 proxy rate), and that a true c > 1 point passes the DMRG stage with probability about 1 (C164's two-XX control reads 1.73 at L = 64). Then zero confirmed detections among N = 12,000 points give

    p < -ln(alpha)/(eta N):  p < 7.3e-4 (95%), 1.1e-3 (99%) with eta = 0.341 (free);  p < 1.4e-3 (95%), 2.2e-3 (99%) with eta = 0.178 (int).

This holds if 10193 and 450 fail the second stage. If one of them survives it, the sample contains a candidate and the bound is void.

Caveats.
- The proxies have range 3 and are close to free, so their finite-size behaviour need not match an interacting nearest-neighbour c = 2 liquid.
- The bound refers to the S^8 measure of C155's coordinates (J diagonal, D, f).
- It says nothing about sets of very small volume, which is where (R_N)(ii) to (iv) and the thin wedges of section 4 would live.

Side remark. [numerical] Only 6 to 17 percent of c = 1 proxies pass crit < 0.15, and the median c_S over all stage-A points is 0.93. So stage A's 5.3 percent "gapless" fraction likely undercounts the floating volume substantially.

## 9. Consistency with earlier runs

- C154: Lemma O reduces to Lemma 3 for product references. The P1 and P3 counts are unchanged, since Lemma 1 confirms they are kinematically allowed.
- C159: Lemma A does not apply in N, and both shifts are free, as C165 4.5 notes.
- C162 and C164: Lemma D contains C164's mass-imbalance result for soliton valleys as a special case and fails exactly where C162 works (K, equal masses).
- C165: Lemma H is consistent with the SU(2)_1 multiplet argument, and (R_N)(i) contains A_S.
- C174: the frontier classification is the dynamical content that R1 allows. P4 is included, and K_II is assigned to (A6)/(H*).
