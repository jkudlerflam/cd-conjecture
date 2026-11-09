This run was done by a Claude Opus 5.5 agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C171 result: referee of C159 and C163 (T and K_I)

The full report is REFEREE.md, the derivations are in PROOF.md, and the log is PROGRESS.md. The scripts are opcount.py and lemmaL_check.py, with outputs in opcount.out and lemmaL_check.out. All runs were local and took seconds.

## Verdict

The individual lemmas of C159 are sound: Lemma 1, Lemma A for one boson, Lemma B, the SU(2)_k exclusions, and every row of table 4.1, which I reproduced with an independent operator enumerator. So are C163's Fact 2, Proposition P and Proposition G. The reduction of T and K_I to (A_W) is not complete. Five routes to c > 1 avoid both (A_W) as stated and every derived exclusion.

- R1 [derived, numerical count]. Three-state Potts order at momentum 2π/3 × a half-floating boson, c = 9/5, needs only ε. That is codimension 1 at every K > 1/2, with no attraction. It is missing from C159 table 4.1, and also from C154 table 4.2, where it has k = 1 in N.
- R2 [derived]. Codimension ≥ 2 candidates are dismissed by a non-sequitur ("must lie near symmetric strata"). The replacement H_T → A6_TK quotes C154 5.2 as "P1/P3" and drops P2, whose T version has only c = 1 neighbours at any K.
- R3 [derived, numerical count]. Lemma L is wrong on its Q2 = nQ1 + π branch. For n = 2 no spin-0 locking operator is relevant at any 1 < K1 ≤ 2.5 tested, and composites with n ≥ 3 are unlocked at moderate K1. So second bosons whose wavevectors lie in Zθ + Zπ survive at codimension 0, and (A_W)(b) covers only wavevectors outside that lattice.
- R4 [derived, numerical count]. A second boson whose θ is Θ-reversed has x2 = 1/(4K2). With a repulsive host (K1 = 0.9, K2 = 0.2) it gives a codimension-0 c = 2 candidate. C159's "x = K, so attractive" translation and C163's restriction to hosts with K > 1 miss it.
- R5 [derived]. A sector on which translation acts by charge conjugation is absent from the table. It reaches codimension 1 only for K1 > 3/2, but at momentum 0 or π, outside (A_W)'s concrete list.

Smaller corrections.
- C163 Fact 1 should be tagged [derived].
- "Every soft wavevector of the U(1) parent lies in Zα + Zπ" is false at finite field, though harmless on T and K_I.
- "Dimerization needs frustration" is contradicted by C154's cyclic stratum and by C159 3.2.
- "A strongly repulsive composite would be gapped" holds only for a commensurate density.
- KI22 is a first-stage pass inside W (K_I, Δ = −0.75, c_two 1.76) sitting on an S_mid jump. It awaits C166's FES (job 17602270, output not yet collected).

## Corrected statement of the residual

(A_W') [assumed] In T and K_I, outside the controlled neighbourhoods, no gapped sector becomes critical or forms an independent liquid while a half-floating liquid is present, at any K and any codimension. This covers Ising modes at 0 or π, Z_N orders with N ≥ 3, and second bosons outside Λ, inside Λ with an irrelevant conversion vertex, θ-commensurate, or translated by charge conjugation.

This is C154's (A6) restricted to T ∪ K_I. What survives of C163's narrowing is "harmonics of order the exchange", which is backed by Proposition P, Proposition G and C160. "Attractive liquids only" does not survive. K_II is not covered by any of this.

## Open

Two lemmas would restore part of the narrowing. (i) Exclude Z_N order (N ≥ 3) independent of the boson at nearest-neighbour bonds, the period-3 analogue of C165. (ii) Show that a codimension-2 c > 1 point in T forces a codimension ≤ 1 c > 1 structure nearby, which is kinematically false for P2, so it would have to come from dynamics.
