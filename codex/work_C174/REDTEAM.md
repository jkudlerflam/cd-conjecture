# C174 red team: the conditional theorem "c ≤ 1 on B_2"

This review was done by Claude Opus 5.5 (agent C174, 2026-11-08), not by Codex. Astra was not used (out of credits).

Target: C154 PROOF 7 (closure lemma, (A1) to (A5'), (H)), the stratification of C154 PROOF 1, and the claim in STATUS and specs/C169_C174_swarm.md that "c ≤ 1 on every symmetric class" plus (H) (with (A_W) on T and K_I, and A_S near the hidden Heisenberg point) gives c ≤ 1 on all of B_2. Read in full: C154, C141. Read RESULT.md of C152, C153, C159, C160, C162, C163, REPORT §9.110 to §9.126, the C143 singular list. Derivations are in PROOF.md, the one numerical check in rank_check.py and rank_check.out.

Severity scale. Critical: the implication fails as stated. High: the implication holds only after an assumption is restated, and the current wording hides where the content is. Medium: a needed input is missing from the stated critical path, or is supported only at generic points. Low: definitional or presentational.

## Summary

No critical hole. The implication of C154 PROOF 7 is formally valid. But it is not sound modulo the three assumptions the critical path names. Two findings change how the result should be read.

1. (R1) Given (A2), the closure lemma moves no upper bound into a component with c > 1. Such a component can only accumulate on controlled singular points. So (H), restricted to the components that matter, is equivalent to the conclusion on N away from the controlled singular set. (H) is not a mild topological input, and no connectivity argument can prove it.
2. (R2) Read literally, (A5) covers every singular point, and C143's amended singular list contains the band-onset hypersurfaces of every stratum, N included. That would put the exclusion of the P3 mechanism (a second onset inside a floating liquid) inside (A5), which is tagged derived. The derived content of (A5) is a much shorter list.

The rest are scope corrections and missing inputs.

## Holes

### R1. (H) carries the whole N problem [high]

Claim. [proved, PROOF 1, given (A2) and gap stability] If c(Σ*) > 1, the closure of Σ* contains no gapped bond and no conformal bond with smaller c. In particular it contains no controlled conformal point (U(1) class, free locus, XYZ, parity, staggered, cyclic, λ-neighbourhood). So for c > 1 components (H) is equivalent to (H*): the closure meets the controlled singular set S_ctrl. With (A5), (H*) is equivalent to "no c > 1 component in N". The theorem therefore reads "c ≤ 1 on the strata, near S_ctrl, and on N implies c ≤ 1 on B_2".

Why it matters. C154 PROOF 7 says upper bounds "propagate from frontiers inward", and RESULT says (H) "fails exactly for a c > 1 island". Both are true, but every c > 1 component that does not touch S_ctrl is an island by Lemma R1. Bounds flow from higher c to lower c (the P1 surface bounds the adjacent floating phase), never the reverse. The real gain is weaker: a c > 1 component in N cannot touch a symmetric stratum anywhere, not even at a conformal point with c ≤ 1. C154's own remark that "(H) is weaker than (A6)" is true only in the uninteresting sense that (A6) also forbids harmless c < 1 theories such as Potts.

Consequence for C168. An attack on (H) through topology or connectivity of critical loci cannot succeed. It has to be a dynamical exclusion of c > 1 fixed points in N away from S_ctrl, which is (A6) in all but name.

Narrowest fix. Restate PROOF 7 as a reduction: c ≤ 1 on B_2 holds if and only if it holds on the strata, on punctured neighbourhoods of S_ctrl, and at the conformal points of N. Replace (H) by (H*) or by (A6), and drop the "topological" framing in STATUS.

### R2. The scope of (A5) [high]

Statement in PROOF 7: "(A5) near every singular point the emergent central charge is at most 1 in every direction". By (A1), singular means z > 1 or extensive degeneracy. C143's amended list (COMP) includes "band-onset hypersurfaces with their junctions" in every stratum, N included. These are 7-dimensional.

What literal (A5) then contains. (i) The boundary of a P3 region is a codimension-1 surface where a second species has a Pokrovsky-Talapov onset inside a c = 1 floating liquid. It is a band-onset point, so literal (A5) forbids the c = 2 region next to it. (ii) The Lifshitz junctions of C154 PROOF 5.4, where Ising surfaces end on PT surfaces in N, are singular and lie in the list, but PROOF 5.4 is conditional on a one-boson infrared there. (iii) PT onsets out of gapped entangled phases (dimerized, Néel, locked tongues) in N. None of these is derived. So literal (A5) is the P3 exclusion relabelled, and (H) becomes nearly vacuous because most components touch some onset surface.

What is derived. Emergent c ≤ 1 on a neighbourhood holds at:
- the Fibonacci points [derived, C141, C143];
- the Bravyi-Gosset set, near h_F given (SC) (C160), and along γ = {-R(φ)} away from soliton-onset layers (C163 Proposition G);
- onsets on the free locus [proved];
- U(1)-class singular points, including saturation and h_c1, with weak harmonics [derived];
- polarized onsets where the one-magnon problem is controlled (R9);
- XZ, in the directions and hyperplane treated (R3).

Narrowest fix. Replace (A5) by an explicit list S_ctrl with a reference for each entry. Move every other singular point (generic onset surfaces in N and their junctions) out of the controlled set, so that (H*) has to reach S_ctrl itself.

### R3. Direction-wise control is not neighbourhood control [high]

The closure lemma at a singular frontier point h* needs emergent c ≤ 1 on a full punctured neighbourhood. The inputs give "for every direction there is a radius", which does not give one radius for all directions. A c > 1 component can approach h* along a cusp tangent to a direction where the radius shrinks to zero. This is a quantifier interchange, and it bites at three places.
- XZ. C143 and C144 control "every direction checked". C157 treats the tuned hyperplane with the weighted variable λ = t/β², and the important regime is that curved, non-radial approach. No weighted blow-up covers the full punctured neighbourhood. The one numerically unresolved region (open-chain c ≈ 1.3 at L = 256, the Néel-cat reading of C166, job 17602270 pending) sits in the XZ neighbourhood, at distance of order 0.2.
- h_F. C160 is a weighted blow-up, but C161 (iii) shows the semiclassical correspondence is not uniform at soliton onsets. There K → 1 while ħ_eff → 0, so the limits ε → 0 and L → ∞ do not commute. C162 covers protected two-valley onsets. Accidental ones near h_F rest on (A_sol′).
- γ. C163 Proposition G holds "away from soliton-onset layers". Those layers are not covered, apart from (A_sol′).

Narrowest fix. For each point of S_ctrl, name a finite set of weighted blow-up charts that covers a punctured neighbourhood, and mark the uncovered layers. Either control those layers (soliton onsets near h_F and γ, curved approaches at XZ) or take them out of S_ctrl, so that they fall under (H*).

### R4. Stratum results are generic-point results, and the theorem needs every point [medium-high]

By Lemma R1, a c > 1 component of N can only meet a stratum at a point where c > 1 on the stratum. Such a point is generically a multicritical point of high codimension inside the stratum, possibly isolated. These are exactly the points that sampling numerics miss and where finite-size crossovers diverge. Specific gaps:
- C152 (parity class). Tricritical Ising points are not located (Set 3 of job 17597010). The Peotta et al. meeting region is unresolved (Set 4). The XZ window is pending. Intermediate couplings far from controlled regions are numerical only, and strong attraction at Δ ≠ 0 is heuristic.
- C153 (staggered box). The large-b corner is heuristic, and the b = 6 and 8 neutral level is undecided (job 17596442).
- Cyclic stratum, (A5'). The Schur argument holds only inside the SU(2)_1 basin. The boundary between the SU(2)_1 region and the ferromagnet-like ordered region was not resolved (PROOF 3.5). The grid is 20 points with spacing 0.3, and an SU(2)_3 curve is "not excluded between grid points". The 53-point grid (job 17594765) was never collected.
- T and K_I, given (A_W). C163 still has three first-stage path flags (KI22, KI26, N34) awaiting FES at χ ≥ 256.

STATUS's sentence "c ≤ 1 is derived at the physics standard on every symmetric class" is stronger than C152 and C153 claim for themselves.

Narrowest fix. For each stratum, an inventory of the junctions (points where two or more critical loci of the stratum meet), with c and the standard of each, plus an in-stratum codimension count showing which c > 1 theories could sit at isolated points. For the cyclic stratum, decide the SU(2)_1/ordered boundary, first order or not.

### R5. The candidate list misses irreducible c > 1 theories [medium-high]

C154 PROOF 5.2 says every candidate with k ≥ 2 "contains P1, P2 or P3 type sub-loci" in its neighbourhood, and PROOF 7 says "every such island needs two independent soft sectors". That is false for irreducible c > 1 CFTs, which form one species and have only gapped or c ≤ 1 neighbours. Examples:
- Z_N parafermions with N ≥ 5, c = 2(N − 1)/(N + 2), which gives 8/7 at N = 5;
- N = 1 superconformal minimal models with 1 < c < 3/2 (81/70 at m = 5);
- W_3 minimal models (6/5, 10/7, ...);
- SU(2)_3 off the cyclic stratum.

The Z_5 parafermion is in C154's own table with k ≥ 2, so its locus would be 6-dimensional in N. Its deformations lead to Z_5-ordered and disordered gapped phases and, in the chiral-clock picture, to a c = 1 floating phase. That is no P1, P2 or P3. Its natural home is the melting of the Z_m clock-locked tongues at m ≥ 5 that C154 PROOF 2.3 and C163 Proposition P predict, at harmonic amplitudes of order one [heuristic, by analogy with chiral clock models]. C159 PROOF (line 194) already notes that its list is not exhaustive for T and K. As worded, (A6) excludes these theories ("a single boson or at most two Majoranas"), but its one-species rationale does not reach them.

Narrowest fix. Add a pattern P4 (an irreducible c > 1 multicritical point bounding a locked tongue) with its codimension count in N, including the chiral spin-±1 operators that N allows. Then either show that tongues with m ≥ 5 melt through c = 1 phases at every amplitude, or name P4 explicitly as part of (A6).

### R6. The sheet K_II is not assigned to any assumption [medium]

(A_W) is stated for T and K_I. (H)/(A6) is stated "on N". K_II (half of the 5-dimensional stratum K, finite-field parents) is in neither. C159 says its liquids "have not been mapped". By C159's own counting (x = 1/(4K) for the commensurate vertex, Proposition C), K_II keeps a c > 1 window in strongly repulsive liquids, K < 2/15 for Ising × boson or 1/(4K_1) + 1/(4K_2) > 2 for two bosons. That window is far from the attractive window (A_W) addresses. C163 §1.2 groups K_II with N.

Narrowest fix. Say explicitly that K_II is under (A6)/(H*), or state an (A_W)-type assumption for the repulsive K_II window.

### R7. The critical path lists three assumptions, and the chain uses more [medium]

Besides (A_W), (H)/(A6) and A_S, the conclusion also rests on the following.
- (A1), the trichotomy gapped, conformal or singular, with finitely many sectors.
- (A2), C141's weak c-theorem. C141 marks the multi-velocity case as assumed, though R11 shows the use here is mild.
- (SC), C160's semiclassical correspondence, which is the basis of the h_F neighbourhood.
- (A_sol′), C162 and C164 on accidental two-valley soliton onsets, for the h_F and γ layers of R3.
- The asymptotic Schrieffer-Wolff step of C163 Proposition G.
- Lemma 3 at spin 1/2, see R9.
- The every-point stratum statements of R4.
- The full-neighbourhood control at XZ of R3.
- Pending numerics: jobs 17597010, 17596442, 17594765 (never collected), 17602270, and the C163 FES second stage.

Narrowest fix. List these in STATUS's "still assumed" block, or fold them into named assumptions.

### R8. Translation broken spontaneously with a larger cell [medium]

The single-species tools assume a one-site uniform reference: C154 Lemma 3, C160 Lemma 3', the single-harmonic band of the free locus, and C152's Lemma P2 (period ≤ 2). In a period-p broken phase there are p bands per cell.

Partial fix supplied. [proved, PROOF 2; numerical check to 2e-14 for p ≤ 7] Only one bond crosses each cell boundary, so the inter-cell Bloch block is supported on 2 × 2 entries, and det h(K) is a degree-2 trigonometric polynomial for every p. This gives c_MF ≤ 1 for every Hartree-Fock-Bogoliubov state of the parity class at any period, and at most two soft momenta at the onset of any periodic product reference. Open: beyond mean field, and whether a ±k0 pair over a periodic reference is one boson or two.

Definitional point [low, but needed]. c must be defined per superselection sector, that is per clustering ground state. Otherwise first-order surfaces between two c = 1 liquids, and cat states, read c = 2. The open-chain c ≈ 1.3 near XZ (C166) is this artefact in practice. Write the definition into (A1).

### R9. Lemma 3 is linear spin waves about a product state [low-medium]

The Lemma 3 band is exact only for the one-magnon problem over an exact eigenstate (U(1)-type saturation) and is controlled near h_F, where ħ_eff is small. A generic polarized phase of N is dressed, and its instability occurs at |B| of order |J|. There spin-1/2 linear spin waves are not controlled, and the dressed quasiparticle band has all harmonics. C154's check (error falling as 1/|B| at |B| = 8 to 32) is deep inside the gapped phase, not at onset. In N the conclusion survives anyway by codimension: two independent minima softening together need an extra tuning, and without reflection or Θ nothing ties k0 to -k0.

Narrowest fix. Re-derive "one soft species at a generic onset in N" by codimension. Tag Lemma 3 as controlled only at exact-saturation manifolds and near h_F.

### R10. Completeness of the strata list [low]

[derived, PROOF 3] The list is complete for on-site product symmetries, unitary or antiunitary, uniform or site-dependent. The degenerate-singular-value case that C154 sketched closes with no new stratum. Not covered:
- spatial reflections composed with on-site maps, dimension ≤ 5, which C154 PROOF 1.4 mentions but does not analyse, and which can protect a chirality order and lower the codimension of P2-type candidates;
- non-on-site symmetries. Kramers-Wannier keeps nearest-neighbour range only at the critical transverse-field Ising point, so this is nearly empty, and the integrable submanifolds of de Leeuw-Pribytok-Ryan are all c ≤ 1.

These matter for the codimension heuristics behind (H), not for the implication.

Narrowest fix. Add the reflection strata to the table with their codimension counts.

### R11. (A2) and the closure lemma's other steps [low, checked]

- "c is constant on Σ" holds by the definition of Σ.
- Gapped points cannot lie in the closure of a critical component, by gap stability.
- The λ-neighbourhood is redundant once (A2) holds at the conformal points of the U(1) class.
- The incommensurate-torus count survives the dense set of rational points, because only finitely many relevant vertices exist at fixed K.
- (A2) is used only at controlled points with c ≤ 1, at most two Majorana species, where it reduces to the finite gap of everything outside the CFT. Its radius is not uniform: it goes to zero near XZ, where v = 5.7ε and the A Majorana has a gap of order ε. This feeds R3.

No fix needed beyond noting these.

### R12. Numerical support for (H) in N [low]

The only followed symmetry-free case, C68's point 1282 joined to its U(1) parent, was never completed. There is an unexplained c_mid = 0.51 at λ = 0.25, L = 64, and job 17594765's path1282 output is not in typhon_out. So "every followed case is connected" has no completed symmetry-free instance.

## Verdict

[derived] The chain is valid as an implication. It is sound modulo a longer list than the critical path names:
- (A1) with c defined per sector, and (A2);
- every-point versions of the stratum results (R4), with K_II assigned (R6);
- (A5) restricted to an explicit S_ctrl with full punctured neighbourhoods (R2, R3), which brings in (SC), (A_sol′) and the asymptotic Schrieffer-Wolff step;
- (A_W) and A_S;
- (H*) or (A6), which must also exclude the irreducible P4 candidates (R5).

(H) is not weaker in substance than the statement being proved on N (R1). The present theorem reduces B_2 to N and proves nothing inside N. The open problem is a dynamical exclusion of c > 1 fixed points in N away from S_ctrl, together with the junction inventories of R4.
