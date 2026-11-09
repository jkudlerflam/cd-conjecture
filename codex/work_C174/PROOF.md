# C174 derivations supporting REDTEAM.md

This work was done by Claude Opus 5.5 (agent C174, 2026-11-08), not by Codex. Astra was not used (out of credits).

Conventions are those of C154 PROOF. Tags: [proved], [derived], [numerical], [heuristic], [assumed]. Three short derivations are collected here. REDTEAM.md uses them.

## 1. What the closure lemma can and cannot transport

Setting. Σ ranges over connected components of the set of bonds flowing to one conformal manifold (C154 PROOF 7). c(h) is the central charge of a superselection sector of h (see REDTEAM R8 for why this has to be per sector). S_ctrl is the set of controlled singular points, and C_ctrl the set of controlled conformal points (the conformal points of the strata of (A3), (A4), (A5') and of the λ-neighbourhood).

Lemma R1. [proved, given (A2) and stability of a finite gap] Let Σ* be a component with c(Σ*) > 1. Then the closure of Σ* contains no gapped bond and no conformal bond h' with c(h') < c(Σ*). In particular it contains no point of C_ctrl.

Proof. A gapped bond with a finite gap (unique or finitely degenerate ground space) keeps a gap on a neighbourhood, so it is not a limit of gapless bonds. If h' is conformal and lies in the closure of Σ*, then (A2) gives a neighbourhood U of h' on which every conformal bond has c ≤ c(h'). U meets Σ*, so c(Σ*) ≤ c(h'). If h' is in C_ctrl, c(h') ≤ 1 < c(Σ*), a contradiction. QED.

Corollary R1. [proved, given Lemma R1] Restricted to components with c > 1, assumption (H) of C154 is equivalent to

(H*) the closure of every component with c > 1 contains a point of S_ctrl.

If the λ-neighbourhood is counted in the controlled set it changes nothing, because it is open, so a component whose closure meets it intersects it, and c = 1 there. Given (A5) on S_ctrl (emergent c ≤ 1 on a full punctured neighbourhood of each point), (H*) is equivalent to "there is no component with c > 1 outside the strata". So (H) together with (A5) is equivalent to the conclusion of the theorem on N, and the closure lemma moves no upper bound from a c ≤ 1 region into a c > 1 component. Bounds travel only from higher c to lower c, so they reach the floating c = 1 phases from P1-type surfaces, never the reverse.

Remark. The closure lemma has one real consequence for c > 1: a c > 1 component in N cannot touch a symmetric stratum, not even at a stratum point that is conformal with c ≤ 1. It must avoid the strata altogether, or meet them only at points where c > 1 on the stratum (which would themselves be counterexamples), or at singular points.

## 2. Periodic cells: the rank-two cut lemma

Lemma R2. [proved] Let a nearest-neighbour chain with one qubit per site be put in a cell of p sites. In any quadratic description (fermionic BdG after Jordan-Wigner in the parity class, or linear spin waves about a period-p product reference), the Bloch matrix is h(K) = A + B e^{iK} + B† e^{-iK}, where B is supported on the two rows of site p of a cell and the two columns of site 1 of the next cell. Then det h(K) is a trigonometric polynomial of degree at most 2 in K, for every p.

Proof. In the Leibniz expansion of det h(K) each term takes one entry from each row and column. A factor e^{iK} comes only from entries in the 2 × 2 block (rows of site p, columns of site 1), and a term can use at most two such entries, one per row. Likewise e^{-iK} comes only from the transposed block, at most two entries. Hence det h(K) = Σ_{n = -2}^{2} d_n e^{inK}. QED.

Consequences.
- [proved] Fermions. det h(K) has at most four zeros on the circle counted with multiplicity, and each gapless Majorana-Weyl pair costs two (one Ising Majorana gives a double zero, as at the critical transverse-field Ising chain, where det = -8(1 - cos K) in suitable units). So every Hartree-Fock-Bogoliubov state of the parity class with any period p has c_MF ≤ 1. This extends C152's Lemma P2 from p ≤ 2 to all p.
- [proved, kinematics only] Bosons. For a stable period-p product reference h(K) ≥ 0, and τ_x h(K)* τ_x = h(-K) gives det h(-K) = det h(K). A nonnegative degree-2 trigonometric polynomial vanishes at no more than two points, so the soft set at onset is {k0, -k0}, {0}, {π} or {0, π}. This is the soft-set count of C154 Lemma 3 for any period. It does not settle whether a ±k0 pair is one boson or two (Lemma 3 needs pairing for that, and C160's Lemma 3' that excludes interior minima was proved only for uniform references).

[numerical] rank_check.py draws 200 random models for each p in {1, 2, 3, 4, 5, 7}, fermionic and bosonic, and Fourier-transforms det h(K) on 64 points. The largest harmonic with |n| > 2 relative to the largest with |n| ≤ 2 is at most 1.8e-14 (rank_check.out).

Limitation. The lemma is quadratic. Beyond mean field a period-p broken phase has p dressed bands per cell, and nothing here bounds the number of soft dressed modes.

## 3. Site-dependent symmetries with a degenerate pair of singular values

C154 Theorem 1 is proved for distinct singular values and the degenerate case is marked "sketched". The sketch closes as follows.

[derived, short algebra] Let s_1 = s_2 ≠ s_3 (or s_1 ≠ s_2 = s_3) with distinguished axis e, so S = s_e e e^T + s (1 - e e^T). From (1.1), S W_j = Q_j S with W_j = O Q_{j+1} O^{-1}, and as in Lemma 1 both Q_j and W_j commute with S^2, so they lie in the commutant G' = O(2) × Z_2 of S (rotations about e, π rotations about axes perpendicular to e, and their negatives). Since they commute with S, S W_j = S Q_j, so W_j = Q_j on the range of S. If S is invertible this is everything. If one singular value vanishes, the two agree on the range and the common determinant fixes the remaining sign. Either way Q_{j+1} = O^{-1} Q_j O, and every O^{-n} Q_0 O^n must lie in G', that is, Q_0 must preserve each of the lines O^n e.
- O e = e (O a rotation about e). Then Q_0 = R_e(α) works for every α and is uniform: the U(1) class (B along e, or B = 0).
- O e = -e (O a π rotation about an axis u perpendicular to e). Then O^{-1} R_e(α) O = R_e(-α), so Q_j = R_e((-1)^j α): the staggered U(1) class (B along e is needed, which with O e = -e and Q_j^T B = B for all j means B along e and invariant, as in a(XX - YY) + uZZ + bZ).
- O e not equal to ±e. Q_0 has to preserve the finite or infinite set of lines {O^n e}. If the orbit is infinite or the lines are not mutually orthogonal, only Q_0 = ±1 survive (Θ at B = 0, stratum T). If O^2 e = ±e with O e perpendicular to e, the three axes e, O e, e × O e form a frame and one obtains the parity class with a staggered extra symmetry, as in Theorem 1 (iii). If O cycles three mutually orthogonal axes, it is the cyclic stratum at β = γ (or the analogous degenerate edge).
- s_1 = s_2 = s_3. J = s O is a twisted Heisenberg bond, a site-dependent rotation of ±Heisenberg, which lies in the U(1) class (XXZ plus DM at the twisted isotropic point), with site-dependent SO(3) symmetry when B = 0.

So no stratum beyond C154's table arises from degenerate singular values. The list is complete for on-site product symmetries, unitary or antiunitary, uniform or site-dependent. It omits spatial reflections composed with on-site maps (C154 PROOF 1.4 notes them, dimension at most 5) and non-on-site symmetries, see REDTEAM R10.
