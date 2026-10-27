# Gapped entanglement bootstrap: what it offers for c ≤ d − 1

Scope. The entanglement bootstrap for gapped states (Shi, Kato, Kim, McGreevy and collaborators) and the gapless 1+1D offshoots by the same group, read for whether any piece can bound the central charge of a critical nearest-neighbour chain on (C^d)^{⊗L} by something linear in d.

Tags. [read] = I opened the paper (arXiv number, section or equation). [recalled] = memory or search snippet only. [own] = my arithmetic or inference, to be checked.

Convention warning. 2404.03725 and 2403.18410 use c_tot = c_L + c_R, so a non-chiral CFT with Virasoro c has c_tot = 2c. Lin-McGreevy (2303.05444) and Li-Lin-McGreevy (2509.04596) use Virasoro c.

## Summary

The gapped bootstrap turns two exact entropy identities on small balls into sectors, fusion spaces, quantum dimensions and a chiral-central-charge formula. Almost every quantitative output is an entropy or a modular-Hamiltonian moment, bounded by the log of a dimension, so alone they reproduce the K log d obstruction. Exact Markov structure also yields a few counting statements bounded by a dimension rather than its logarithm. Those are the right size for d − 1, but none is tied to c. The 1D descendant, the vector fixed-point equation, is the closest tool (2509.04596 cites c-d explicitly), and it inherits the log obstruction in a clean, computable form (§7).

## 1. Axioms A0 and A1

Statement. [read 1906.09376 §2.1, eqs. (2)-(4)] For a reference state σ and every ball of radius ≤ r (a "µ-disk"), A0 says S(σ_BC) + S(σ_C) − S(σ_B) = 0 with C a disk and B the shell around it, and A1 says S(σ_BC) + S(σ_CD) − S(σ_B) − S(σ_D) = 0 with the shell cut into two arcs B, D.

Hypotheses and mechanism. Exact equalities, 2D, finite-dimensional sites, shapes fixed up to topology, r constant. A0 implies I(A:C) = 0 and A1 implies I(A:C|B) = 0 for A outside the ball, and both are strictly stronger (an infinite-temperature Gibbs state has I(A:C) = 0 but violates A0) [read §2.1, eqs. (5)-(6)]. SSA and recovery extend both to all disks [read Prop. 3.3].

Transfer. In 1D the A1 combination on contiguous B C D is the weak-monotonicity combination Δ. For a CFT ground state on contiguous a, b, c, Δ(a,b,c) = −(c_tot/6) ln η and I(a:c|b) = −(c_tot/6) ln(1 − η) [read 2404.03725 §1.2, Example 4.2]. So the axioms fail at criticality by an amount proportional to c, and everything downstream that needs exact A0/A1 is unavailable. The violation is itself an entropy, hence O(ℓ ln d) on ℓ-site blocks [own].

## 2. Merging lemma, Markov block structure, commuting supports

Statement. Merging lemma [read 1906.09376 Lemma 3.8, eqs. (16)-(19); restated 2112.08398 Lemma 2.6]. If ρ_BC = σ_BC and I(A:C|B)_ρ = I(B:D|C)_σ = 0, there is a unique τ_ABCD with marginals ρ_ABC and σ_BCD, with I(A:CD|B)_τ = I(AB:D|C)_τ = 0, and ρ ↦ τ preserves entropy differences. Structure theorem [read 2404.05867 Thm 3, eqs. (51)-(52), citing Hayden-Jozsa-Petz-Winter]: I(A:C|B) = 0 iff H_B = ⊕_j H_{b_j^L} ⊗ H_{b_j^R} and ρ_ABC = ⊕_j q_j ρ_{A b_j^L} ⊗ ρ_{b_j^R C}. Consequence: the support projectors of ρ_AB and ρ_BC commute [read 2404.05867 Prop. 3, eq. (54)].

Hypotheses and mechanism. Exact zero CMI, finite dimensions. The merged state is the Petz map σ_CD^{1/2} σ_C^{−1/2} (·) σ_C^{−1/2} σ_CD^{1/2} applied to ρ_ABC [read 1906.09376 Lemma 3.2 cites Petz; explicit form recalled]. Prop. 3 plus A0/A1 gives a local commuting-projector parent Hamiltonian with a stable gap [read 2404.05867 Thm 1].

Transfer. The block structure is a counting statement. With B a single site there are at most d blocks and Σ_j dim b_j^L · dim b_j^R ≤ d [own, immediate]. That is the scale c-d wants. The obstacle is that a critical chain is never Markov on contiguous intervals (CMI ∝ c above), and approximate recovery bounds (Fawzi-Renner type [recalled]) live in the entropic currency, which returns log d. A possible use would be exact or near-exact Markov structure on a subalgebra, for instance a neutral sector, which would bound its block count by d. I know of no paper doing this [own].

## 3. Information convex set and isomorphism theorem

Statement. Σ(Ω) is the set of states on Ω obtained by partial trace from states on a thickening Ω' that agree with σ on every µ-disk inside Ω' [read 1906.09376 Def. 3.1, eq. (13)]. It is compact and convex [Prop. 3.4], a single point for a disk [Prop. 3.5], and {|ψ⟩⟨ψ|} for S² [Prop. 3.7]. Isomorphism theorem: regions joined by a path of elementary extensions and restrictions have information convex sets related by a linear bijection preserving distances, fidelity and entropy differences [read Thm 3.10, Cors. 3.9.1-3.9.2]. Huang-McGreevy-Shi extend this to immersed regions and regular homotopies [read 2112.08398 Def. 2.1, Thm 2.5]. Shi-Huang-McGreevy immerse punctured closed manifolds into a ball and heal the puncture ("completion trick"), deriving remote detectability, i.e. unitarity of a pairing matrix, as a theorem [read 2301.07119 abstract, Lemma 4.4].

Hypotheses and mechanism. Exact A0 and A1 along the path. Each elementary step is a merge with σ on a ball or a partial trace, both injective under the Markov conditions.

Transfer. Poor. In 1D an interval admits no nontrivial immersion into a line, and every step needs exact Markov [own].

## 4. Simplex theorem, sectors, fusion spaces, quantum dimensions

Statement. For an annulus X, Σ(X) is a simplex whose finitely many extreme points σ^a_X are mutually orthogonal; these define the superselection sectors [read 1906.09376 Thm 4.1, eq. (31)]. For a 2-hole disk, each Σ^c_ab(Y) is the state space of a fusion space V^c_ab, and N^c_ab = dim V^c_ab obeys Kitaev's fusion axioms [read Thms 4.4-4.5, Def. 4.3, Props. 4.6-4.11]. The entropy difference f(a) = [S(σ^a_X) − S(σ^1_X)]/2 equals ln d_a with d_a d_b = Σ_c N^c_ab d_c [read Def. 4.2, eq. (36); Def. 5.1, eqs. (58), (60)]. The Levin-Wen CMI equals 2 ln D and the Kitaev-Preskill γ equals ln D [read Props. 5.1-5.2, eqs. (61)-(62)], and the maximum-entropy element of Σ(X) is Σ_a (d_a²/D²) σ^a_X [read proof of Prop. 5.1]. On gapped domain walls the same machinery yields "parton" sectors with S_topo,N = 2 ln D_N and S_topo,U = 2 ln D_U [read 2008.11793 eq. (1), Fig. 12, Thm III.1].

Hypotheses and mechanism. Exact A0 and A1; finite local dimension only makes the label set finite. Cut X = LMR into thin annuli. Extreme points have I(L:R) = 0, so the fidelity factorizes, F_LR = F_L F_R. The isomorphism theorem gives F = F_L = F_R = F_LMR ≤ F_LR = F², so F ∈ {0, 1} [read eqs. (32)-(33)].

Transfer, and the counting flag. This is the main source of dimension-type bounds. Because F_L = F, distinct sectors are orthogonal already on the thinnest sub-annulus, so N ≤ dim H_L ≤ d^k for a width-k annulus [own corollary]. Also N ≤ D² = exp I(A:C|B) ≤ d^{2|A|} [own, from Prop. 5.1, SSA and d_a ≥ 1]. Fusion multiplicities are bounded by the dimension of the 2-hole disk. These are polynomial in d at fixed region size, the right scale for d − 1. The trouble is that sector data does not control c. E8 level 1 has one sector and c = 8, while SU(d)_1 has d sectors and c = d − 1 [recalled, standard]. In 1D the sector analog is the DHR sectors of the chiral algebra, whose index μ = Σ_a d_a² is a two-interval quantity [recalled, Kawahigashi-Longo-Müger], so a 1D simplex theorem would bound μ by d^{O(ℓ)}, not c. It could help only inside a class where c is itself a count, such as c = rank for Gaussian or WZW theories (compare the susceptibility-rank result in REPORT.md).

## 5. Modular commutator

Statement. J(A,B,C)_ρ = i Tr(ρ_ABC [K_AB, K_BC]) with K = −ln ρ [read 2110.06932 eq. (2)]. For A, B, C partitioning a bulk disk and meeting at a point, J = (π/3) c_- [read eq. (3)]. If I(A:C|B) = 0 then J = 0 [read eqs. (4)-(5)]. Under the modular flow of K_BC, dS(ρ_AB(s))/ds at s = 0 equals J [read 2110.10400, around eq. (15)]. J agrees on the two sides of a gapped domain wall [read 2110.10400 Thm 1]. For contiguous intervals on a time slice of a 1+1D CFT vacuum, J = (π c_-/6)(2η − 1) with c_- = c_L − c_R [read 2206.00027 eq. (1)], and for a "chiral thermal state" e^{−β_L H_L − β_R H_R} with c_L = c_R = c, J = (π c/3)(η_eff^{(β_L)} − η_eff^{(β_R)}) [read eqs. (10)-(11)].

Hypotheses and mechanism. Deformation invariance is rigorous given A1, since Markov decompositions K_ABC = K_AB + K_BC − K_B let one slide boundaries away from the triple point [read 2110.06932, Fig. 2 argument]. The identification with c_- is a physical argument, not a theorem: K_D decomposes into local terms (eq. (7)), its conserved "modular current" (eq. (10)) is matched to the edge current (π/12) c_- T² by treating K_D as a physical Hamiltonian [read pp. 2-3].

Transfer. On a time slice J sees only c_-, which vanishes for the chains in scope (lit.md agrees). The chiral-thermal formula sees c itself, but building such a state needs H_L and H_R separately, and these are not local lattice operators [own]. Dimension bounds on J are polylogarithmic: |J| ≤ 2 min(log d_AB, log d_BC) Λ(I(A:C|B)) [read 2609.09019 Lemma 2.8, eq. (2.25)].

## 6. Finite dimension versus chirality

Statement. No state on a tensor product of finite-dimensional sites satisfies bulk A1 with J ≠ 0, and none satisfies A1 with on-site U(1) and a nonzero Hall estimator [read 2408.10306 Thms 3.1-3.2]. Finite dimension can be relaxed to finite local entropies [read App. F]. In the approximate setting, approximate local Markov makes J approximately deformation invariant, and if J ≠ 0 the strip CMI cannot decay faster than exponentially in the width of B [read 2609.09019 abstract, Thm 2.9, Cor. 4.8]. Kim-Lin-Ranard-Shi show that strict area law, or A0/A1, implies a commuting-projector parent Hamiltonian [read 2404.05867 abstract, Thm 1].

Mechanism. The instantaneous modular flow I_A(t)|ψ⟩ = ρ_A^{it}|ψ⟩ uses ρ_A from the current state [read 2408.10306 Def. 4.1], and Markov states admit a decomposition move [Lemma 4.2]. Along |Ψ(t)⟩ = I_AB(t)|Ψ⟩, A1 keeps J = π c_-/3 for all t [Lemma 5.1], so S_BC(t) = S_BC(0) + π c_- t/3 [eq. (5.4)], which is unbounded and contradicts S ≤ ln dim H_BC.

Transfer. This is the one place where finite d enters as a contradiction rather than a ln d estimate, so the size of ln d is irrelevant. It needs a quantity that drifts linearly forever under a flow preserving the hypotheses. For a non-chiral 1D CFT the time-slice J vanishes and the modular flow of an interval is geometric, so S_BC(s) saturates as the flowed endpoint approaches the tip of the causal diamond [own, from the vector field in 2206.00027 Fig. 2]. I see no 1D quantity with a drift proportional to c [own].

## 7. The 1+1D descendants: vector fixed-point equation and stationary c_tot

Statement. For contiguous intervals A, B, C with cross ratio η, a 1+1D CFT ground state satisfies K_Δ|ψ⟩ ∝ |ψ⟩ with K_Δ = (K_AB + K_BC) − η(K_A + K_C) − (1 − η)(K_B + K_ABC), and S_Δ = ⟨K_Δ⟩ = (c/3) h(η) with h the binary entropy. Equivalently, |ψ⟩ is a critical point of S_Δ, which is ≥ 0 as a convex combination of weak monotonicity and SSA [read 2303.05444 eqs. (1)-(4), (9)]. Interval EHs then decompose into 1- and 2-cell EHs, and the global conformal generators follow, e.g. H = (1/π) Σ_j K̃'_[j,j+2] [read eqs. (14)-(17)]. The reconstructed Hamiltonian H_rec ∝ Σ_i (K_{i,i+1} − K_i) is nearest-neighbour and nonnegative by the Lin-Kim-Hsieh operator weak monotonicity [read 2509.04596 eqs. (1.4)-(1.5), citing 2211.13372]. Kim-Li-Lin-McGreevy-Shi define c_tot(D) and η(D) by e^{−6Δ/c_tot} + e^{−6I/c_tot} = 1 on a "conformal ruler" D and prove that δc_tot = 0 under all state variations iff the VFPE holds, a locally checkable condition [read 2404.03725 Def. 4.1, Assumption 4.6, Def. 4.7, Thm 4.8, eq. (4.13)]. With bulk A1 and c_tot > 0, c_tot is the same on every ruler iff the η(D) form a consistent set of cross ratios [read Props. 5.5-5.6]. The non-chiral case needs a genericity assumption [read Assumption 6.1]. In a limit c_tot reduces to Casini-Huerta's 6r∂_r S [read Example 4.3].

Hypotheses and mechanism. The CFT statements are continuum ones; on lattices the VFPE holds approximately, with power-law decaying error [read 2303.05444 p. 2]. The interval modular Hamiltonian is a parabolic-weight integral of T_00, and the weights in K_Δ cancel identically [read eqs. (6)-(7)]. Stationarity ⇔ VFPE is first-order perturbation theory of entropies. RG monotonicity of c_tot, and c_tot ≥ |c_-|, are open [read 2404.03725 §7.2].

Small-space search. Li-Lin-McGreevy look for critical points of c_Δ = ⟨ĉ_Δ⟩, ĉ_Δ = (3/(2 ln 2)) Σ_i (K_{i,i+1} − K_{i+1}), among translation- and time-reversal-invariant states of 4 qubits and 4 qutrits [read 2509.04596 eqs. (2.1)-(2.2)]. Qubits give Ising 0.526, XX 1.211, Heisenberg 1.245, and a perfect-tensor "maximum point" at 6.000 [read Table 1]. Qutrits give tricritical Ising 0.744, Potts 0.844, and unidentified candidates with 1 < c < 2 [read §3]. They call the search a stringent test of c-d, report that c_Δ always overestimates c, and find inferred c below d − 1 [read §3].

Transfer, and the log obstruction made explicit. With single-site A, B, C on a line (η = 1/4), Δ ≤ 2S_B ≤ 2 ln d and I(A:C|B) ≤ 2 ln d, so S_Δ ≤ 2 ln d, and an exact lattice-scale VFPE would force c ≤ 6 ln d / h(1/4) ≈ 10.7 ln d. With ℓ-site blocks, c ≤ 6ℓ ln d / h(η) [own]. On the 4-site circle, S_{i,i+1} − S_{i+1} ≤ S_i ≤ ln d gives c_Δ ≤ 6 log₂ d, saturated at d = 2 by their maximum point [own, consistent with Table 1]. So their test can return c_Δ > d − 1 only for d ≤ 30 (6 log₂ 31 ≈ 29.7 < 30), and their empirical c_Δ ≥ c must fail for SU(d)_1 at L = 4 once d ≥ 31 [own]. In this language c-d becomes a statement about the onset scale ℓ*(d) at which the VFPE becomes accurate: SU(d)_1 needs ℓ* ≳ (d − 1) h(η)/(6 ln d) [own]. Controlling ℓ* from the nearest-neighbour bond h is not attempted in any of these papers.

Other papers. 2606.12540 extends H_rec to any dimension, numerics only [read abstract]; 2510.23720 treats chiral gapped corners [read abstract]; 2309.17155 and 2506.04221 seem irrelevant [recalled]. "On the strength of the vector fixed-point equation" (cited as to appear, 2509.04596 ref. [20]) reportedly shows a VFPE state is the ground state of its H_rec; I could not find it on arXiv.

## 8. Ranking for c-d

(i) VFPE and stationary c_tot with the nearest-neighbour H_rec (§7). It is the only native 1D tool, it defines c from a single state, the reconstructed Hamiltonian is nearest-neighbour like the c-d hypothesis, and it already serves as a c-d search engine. It is entropic, so c ≤ 6ℓ* ln d / h(η). Best uses are a counterexample search at d ≤ 30 and a reformulation of c-d as a bound on ℓ*(d).

(ii) Markov block structure and commuting supports (§2). Counts ≤ d per site. It needs exact or controlled Markov structure on some subalgebra, and no route exists yet.

(iii) Simplex-theorem counting (§4). Sectors are orthogonal on the thinnest region, N ≤ d^k. Right scale, wrong observable, unless paired with a class where c is a count.

(iv) The instantaneous-modular-flow no-go (§6). The only argument in which finiteness yields a contradiction independent of the size of ln d. It needs a 1D quantity with drift ∝ c under a hypothesis-preserving flow, and I know of none.
