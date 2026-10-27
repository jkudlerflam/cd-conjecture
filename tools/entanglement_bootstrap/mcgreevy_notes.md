# McGreevy's lecture notes and reviews as a toolkit for the c–d conjecture

Digest of John McGreevy's course notes, TASI lectures and review, read for anything that bounds a 1+1d central charge from lattice data. Tags: [read: source, section, page] means I found the statement in the PDF. [recalled] means from memory, not checked against these texts. [computed] means my own short derivation from the stated ingredients, to be rechecked before use.

## Sources

- QI: "Quantum information is physical," Phys 213, Winter 2023, mcgreevy.physics.ucsd.edu/w23/2023-213-lectures.pdf. The same chapters appear in the Phys 239 notes (f19, s18).
- TASI: "TASI Lectures on Symmetries in Quantum Matter" (2023, updated Jan 2024), mcgreevy.physics.ucsd.edu/talks/2023-TASI-lectures.pdf. Section 5 is the entanglement bootstrap for 1+1d CFT.
- GS: "Generalized Symmetries in Condensed Matter," arXiv:2204.03045, Ann. Rev. CMP 14 (2023) 57.
- QPM21, QPM24: Quantum Phases of Matter notes (s21/2021-211C, s24/2024-230).
- QFTsrc: "Where do QFTs come from?", Phys 239a, Spring 2014 (s14/239a-lectures.pdf).
- 215C: QFT part 3, Spring 2025, ch. 5 (generalized symmetries).
- s-src: "Hierarchical growth of entangled states, or s-sourcery," talk slides (2018).

Phys 217 (RG), Phys 220 (group theory), 215C (2026) and the 2020 UQM lecture add nothing relevant.

## Bottom line

Nothing in these notes bounds c linearly in the on-site dimension d. Every entropic or information-theoretic tool he presents gives c = O(log d) when evaluated on a fixed number of sites. Three things are useful. (i) A lattice-only argument from strong subadditivity (SSA) shows how the log d obstruction degrades to c ≤ 3ℓ* log d, where ℓ* is the length at which conformal behavior sets in. A thermodynamic argument gives the same scale independently. (ii) The entanglement-bootstrap c-functional comes with a stationarity theorem, and it is the sharpest local characterization of "CFT ground state" in his notes. It ignores charges, so it sees neutral sectors. (iii) LSMOH, anomaly and non-invertible-symmetry arguments give only lower bounds, but they are the only place where integer data tied to d enters.

## 1. SSA on the lattice: monotonicity, concavity, extensivity

Statement [read: QI §5, p.116; §6.4, p.136]. For a translation-invariant state of an infinite system, SSA makes the slab entropy S(ℓ) nondecreasing and concave in ℓ, and S(ℓ)/ℓ has a limit. No Hamiltonian, gap or Lorentz invariance is assumed. He notes that the raw inequalities hold instance by instance and therefore survive disorder averaging.

Mechanism. Apply S(X)+S(Y) ≥ S(X∖Y)+S(Y∖X) and S(X)+S(Y) ≥ S(X∩Y)+S(X∪Y) to overlapping translates. On a finite ring the argument fails once 2ℓ > L, since the translates touch on the far side.

For c–d [computed]. Concavity plus subadditivity give S(ℓ+1) − S(ℓ) ≤ S(1) ≤ log d for every ℓ. If S(ℓ) = (c/3) log ℓ + const holds for ℓ ≥ ℓ*, then (c/3) log(1 + 1/ℓ*) ≤ log d, so

c ≤ 3 log d / log(1 + 1/ℓ*) ≈ 3ℓ* log d.

With ℓ* = 1 this is c ≤ 3 log₂ d. For SU(d)₁ to be consistent, its conformal onset length has to satisfy ℓ* ≳ (d−1)/(3 log d). Attacked entropically, the conjecture becomes a statement about how late the CFT regime is allowed to set in. This is an upper bound, but only once ℓ* is controlled.

## 2. Casini–Huerta entropic c-theorem

Statement [read: QI §9, p.188-189; TASI §5.4, p.78]. In the ground state of a Lorentz-invariant 1+1d QFT, c(r) = 3r S′(r) is nonincreasing in r and equals c at a fixed point.

Mechanism. A Minkowski analogue of Ptolemy's relation gives bc = ad for the boosted intervals. SSA is applied to the two boosted intervals whose union and intersection are the outer and inner intervals, and Lorentz invariance identifies their entropies with those of intervals of equal proper length. Area-law terms cancel in the differences.

For c–d. This gives c_IR ≤ c(r) at any r where boost invariance holds. On the lattice, c(r) at r of a few sites is at most about 3r log d (by §1), so the bound is again O(log d). The hypothesis that fails at the lattice scale is exact boost invariance. A Lieb-Robinson cone (§5) cannot stand in for it, because the step relating states on two slices of one causal diamond needs a unitary map between them. He stresses that this is a relativistic statement [read: QI §6.4, p.134]. Upper bound, O(log d).

## 3. Entanglement-bootstrap c-functional and the vector fixed-point equation

Statement [read: TASI §5.1, p.66-69; theorem §5.2, p.73; §5.4, p.78]. Take contiguous A, B, C with cross ratio η and define

S_Δ(ψ) = S_AB + S_BC − η(S_A + S_C) − (1−η)(S_B + S_ABC) = η Δ(A,B,C) + (1−η) I(A:C|B) ≥ 0,

with Δ = S_AB + S_BC − S_A − S_C. In a 1+1d CFT vacuum, the matching combination K_Δ of modular Hamiltonians is a c-number, K_Δ|ψ⟩ ∝ h(η)|ψ⟩ with h the binary entropy. This follows from the single-interval modular Hamiltonian (Bisognano–Wichmann plus a conformal map). The theorem: K_Δ|ψ⟩ ∝ |ψ⟩ holds if and only if the extracted c(ψ) is stationary under every norm-preserving variation of the state. The conjecture: stationary points for all η are RG fixed points. The same identity merges ρ_AB and ρ_BC into ρ_ABC (it reduces to Markov merging at η = 0). It also yields a parent lattice Hamiltonian H_rec = Σ_i (K_[i,i+2] − K_[i,i+1]), which he says works even for small local Hilbert spaces. For one choice of regions, c(ψ) reduces to the Casini–Huerta function.

Mechanism. dS = ⟨dψ|K|ψ⟩ + c.c., because the δK term drops out when tr ρ = 1 is preserved. Stationarity of S_Δ is therefore the vector equation.

Numerics [read: TASI p.69]. Take A, B, C and the complement each to be one qubit on a 4-site ring. The stationary points then include approximants to Ising (c ≈ 0.53), XX, Heisenberg and a Fibonacci chain, and the maximum is c = 6.000. A note on normalization: the table uses c(ψ) = 3S_Δ/h(η), so Ising gives about 1/2, while eq. (5.2) and §5.4 label the coefficient c_L + c_R. I follow the table.

For c–d [computed]. On the 4-site ring with d-dimensional sites (η = 1/2), we have Δ ≤ 2S_B ≤ 2 log d and I(A:C|B) ≤ 2 log d. So S_Δ ≤ 2 log d and c(ψ) ≤ 6 log₂ d, with equality for Bell_d(1,3) ⊗ Bell_d(2,4). For d = 2 this is the table's "maximum point." Blocks of ℓ sites give O(ℓ log d). SU(d)₁ exceeds 6 log₂ d once d ≥ 31, so its ground state cannot be resolved at four sites. Even so, this looks like the best tool here for neutral sectors. It never refers to on-site charges, it sees Majorana and parafermion content only through c, and it gives a variational characterization of c. A linear bound would need a non-entropic constraint on which states of (C^d)^⊗ℓ can be stationary points with a nearest-neighbor H_rec. Upper bound, O(ℓ log d).

## 4. Small incremental entangling, entropy-rate bounds, first law

Statement [read: QI §6.3, eq. 6.20, p.131-132]. For a pure state on aABb evolved by e^{iH_AB t}, the rate obeys |∂_t S(Aa)| ≤ c‖H‖ log D with D = min(dim A, dim B) (Bravyi, with ancillas by Van Acoleyen et al.). Related results [read: QI §6.4, p.133-136] are the Afkhami-Jeddi–Hartman rate bound |∂_t Ŝ| ≲ (area) × s_T from relative-entropy monotonicity plus Lorentz invariance, the first law δS = δ⟨K⟩, and Casini's Bekenstein bound. In the Bekenstein bound, joint convexity of relative entropy stops species proliferation from violating it.

Mechanism. Maximize the trace norm of [log ρ_A, |Ψ⟩⟨Ψ|] over Schmidt spectra.

For c–d. The entangling rate through one nearest-neighbor bond is O(‖h‖ log d), so any construction of the critical ground state by local or adiabatic dynamics inherits log d per bond per unit time. The Bekenstein remark bounds entropy, not c. O(log d) only.

## 5. Lieb-Robinson bounds

Statement [read: QI §6.5, p.139; QPM24 §1.4]. For local H with terms decaying exponentially in diameter, ‖[A_X(t), B_Y]‖ ≤ c e^{−a d_XY}(e^{2st} − 1), and the Lieb-Robinson velocity is 2s/a. He uses it for Hastings' 1d area law and for the finite-depth-circuit definition of gapped phases.

For c–d. It gives an emergent cone with v ≤ v_LR = O(‖h‖), and for nearest-neighbor chains this is independent of d [recalled]. That caps the CFT velocity, but it does not count modes. Neutral.

## 6. MERA and s-sourcery

Statement [read: QI §6.1, footnote 45, p.123-124]. S(A) is at most the number of MERA bonds cut times log χ. In 1d the minimal cut crosses O(1) bonds per layer over log R layers, and the coefficient is set by the hyperbolic radius. [read: QI §6.1, p.125; s-src] A size-doubling "s-source" circuit implies S(2R) ≤ sS(R) + kR^{d−1}, and CFTs are expected to be s = 1. Swingle–McGreevy entanglement thermodynamics (1505.07106, in s-src) bounds S(ρ_A) by the entropy of the maximum-entropy local Gibbs state with the same ⟨H_x⟩.

For c–d [computed]. A binary MERA gives S(R) ≤ 2 log₂ R · log χ, hence c ≤ 6 log₂ χ, so SU(d)₁ needs χ ≥ 2^{(d−1)/6}. The 1d s = 1 recursion gives c ≤ 3k/log 2, with k set by the entanglement a doubling circuit adds per boundary, which is again a log of local dimensions times depth. These are consistency checks on the obstruction, not ways around it.

## 7. CFT thermodynamics versus lattice entropy

Statement [read: QFTsrc §4.3.2, p.131]. Scale invariance forces F = −c V T^{D} in D spacetime dimensions, with c a rough count of degrees of freedom that is proportional to Virasoro c in 1+1d. In 1+1d, s(T) = πcT/(3v) [recalled].

For c–d [computed]. The lattice entropy density satisfies s(T) ≤ log d and is increasing in T. So if CFT thermodynamics holds up to T*, then c ≤ 3v log d/(πT*). Equivalently, the thermal length at crossover obeys v/T* ≥ (π/3)(d−1)/log d for SU(d)₁, whatever v is. This matches the entanglement onset scale of §1. Two independent lattice routes say the same thing: SU(d)₁ only shows up on scales of order d/log d sites. A proof of c ≤ d−1 must either be non-entropic or control this onset length.

## 8. LSMOH and anomaly matching (lower bounds)

Statement [read: QPM21 §1.8, p.35-38]. With U(1) and translations, gauge invariance of the term (ν/2π)∫A∧dθ forces ν ∈ ℤ particles per unit cell, unless (i) extra gapless modes make the effective action nonlocal or (ii) k-fold degeneracy (topological order or translation breaking) allows ν ∈ ℤ/k. Flux threading gives the spectral-flow version. [read: GS §3, p.16-17] An anomaly is a property of the Hilbert space and the symmetry action, independent of H and RG-invariant, and it forbids a trivially gapped symmetric ground state. SPT edges must be gapless, symmetry-broken or topologically ordered, "the same conditions arising from" LSMOH. [read: QPM24 §3, p.84] Translation can be anomalous on a tensor-product lattice because it is not on-site. [read: GS §3.2] Oshikawa's argument generalizes to 1-form symmetries.

For c–d. These force something nontrivial in the IR, which is a lower bound. For an SU(d) chain with a fundamental on each site, the PSU(d) × translation anomaly is where d enters as an integer. Matching it in an SU(d)-symmetric critical phase is argued to require at least SU(d)₁, so c ≥ d−1 [recalled: Yao–Hsieh–Oshikawa, PRL 2019, not in McGreevy, check]. Together with the conjectured upper bound, this would pin the ULS chain. A generic chain with no on-site symmetry has no LSM anomaly, so the tool says nothing there.

## 9. Generalized and non-invertible symmetries

Statement [read: 215C §5.1, p.146-148]. Of his five properties of symmetry operators, the only essential ones are [H, U] = 0 and locality preservation. Dropping invertibility gives fusion-category symmetries [read: 215C §5.9, p.167-168; GS §5, p.21-24], for example the Ising Kramers–Wannier defect with N² = 1 + η. RG flows driven by local operators can only generate operators that pass through the defect. Modular-invariance argument [read: GS eq. 38, p.24]: if a line's expectation value cannot be a non-negative integer, a unique gapped vacuum is excluded (for Fibonacci W² = 1 + W, the number of vacua must be even). Every irrep of a finite G appears in the spectrum of a 1+1d CFT (GS ref. 110). Emergent higher-form symmetries have exact consequences, unlike emergent 0-form ones [read: GS §2.6; TASI §1]. Under Hilbert-space fragmentation the commutant algebra grows exponentially with size, while for ordinary symmetries it grows polynomially [read: GS §7, p.27].

For c–d. This is the natural language for the neutral sectors, which are labeled by Verlinde-type lines (KW for Ising, ℤ_N duality for parafermions) rather than by U(1) charges. McGreevy gives only statements of the "no trivial gap" type. A linear-count route, which is my speculation and not in the notes: the IR lines acting on neutral sectors may have to be realized by locality-preserving lattice operators (MPOs) commuting with H. On (C^d)^⊗L with nearest-neighbor H, the commutant of the bond algebra might then bound how many independent neutral sectors can appear.

## 10. Counting gapless modes: Luttinger as an anomaly, additivity of c

Statement [read: GS §6.2, p.26]. Fermi and non-Fermi liquids carry an emergent loop-group symmetry, an independent U(1) at each Fermi point. Luttinger's theorem is the anomaly of this symmetry, matched to the microscopic filling (Else–Thorngren–Senthil). In 1d this is U(1)^n_L × U(1)^n_R for n pairs of Fermi points. [read: QFTsrc §4.2.1, p.128] c is additive, a free boson has c = 1 and a Majorana c = 1/2. [read: QPM24 §2.2, p.48-49] Two chiral Majoranas bosonize to one chiral boson, and the SU(N)_k WZW model has c = k dim G/(k+N), which gives c = d−1 for SU(d)₁. [read: QFTsrc §0.1, p.5-6] He lists six mechanisms for gaplessness, with criticality and "CFT with no relevant operators" as separate entries.

For c–d. The loop-group anomaly counts independent emergent U(1)s. It is the closest thing to a rank bound in these notes and lines up with the known susceptibility-rank result. It cannot see neutral sectors.

## On gapless degrees of freedom per d-state site

I found no quantitative statement. He calls c "a measure of the number of critical degrees of freedom" (QI §6, p.123). The closest thing to a per-site count is the 4-qubit bootstrap table, whose maximum c = 6 is the entropic ceiling 6 log₂ d of §3 at d = 2.

## Ranking for c–d

1. Entanglement-bootstrap c-functional with the stationarity theorem (§3). It ignores charges and gives a local, variational handle on c, which makes it the right object for neutral sectors. Its log d ceiling is explicit, so any improvement has to come from a constraint tied to H_rec.
2. SSA concavity plus the onset-length reformulation (§1, cross-checked by §7). It recasts the obstruction as c ≤ 3ℓ* log d and suggests proving ℓ* = O(d/log d) or finding a non-entropic route.
3. Non-invertible and generalized symmetries (§9). They are the only framework here that labels neutral sectors, but so far they give lower bounds only.
4. LSMOH and anomaly matching (§8). Linear-in-d lower bounds where on-site symmetry exists, which could pin the saturating case.
5. Luttinger as an anomaly of the emergent U(1)s (§10). A rank-type count for charged sectors, complementing the known susceptibility result.
