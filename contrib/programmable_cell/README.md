# Programmable cells at fixed d: port-idle history states

## Claim

Class: translation-invariant, nearest-neighbour, Hermitian chains H = Σ h_{i,i+1} on (ℂ^{d0})^⊗L, h depending on a target size N but not on L. Question (REPORT §9.1 to 9.3, C9 §5): can spontaneously formed cells reproduce a prescribed gadget (code, port compressions and, for the second-order route, the boundary resolvent) with the same h on every bond?

1. History-state cells need not dilute their ports. If no gate acts on the port and the port enters only through energy terms, the code and its port compression are exact, and the boundary resolvent R(z) and two-cell exchange K2 converge to the target's as O(L⁴h_max²/(J_c Δ²)), L = T + 1. One gate on the port leaves an error that survives J_c → ∞. [proved by Astra P1 and checked, numerical]
2. The published 1D universal Hamiltonians (KPBC, 2003.13753) already preserve input-track observables to O(η), removing the C9 dilution term T_PE/T by idling. Dilution was a cost, not an obstruction. [read]
3. Two gadgets fit this interface. (A) The G1/C9 second-order gadget, code |Ω⟩_{p,s} ⊗ |i⟩, enforced by one frustration-free projector. (B) A new "dressed Bell-pair" cell for the first-order route, with an exactly degenerate 2^n-dimensional ground space for every coupling and ports compressing exactly to αγ_a and αiγ_aP on one Clifford algebra (C8's joint form) when the channels are odd under a port grading. [proved, numerical]
4. A chain is specified at the level of tracks and terms, at d0 = 4k_c q, with k_c the alphabet of a Gottesman–Irani-type computation track (probably 10² to 10³), independent of q. Route A gives c = N − 1 at d0 = 8k_c for every N. Route B gives c = q²/4 > d0 − 1 once q > 16k_c. [heuristic]
5. The cell is no longer where the construction stops. What remains is writing out the computation track, infrared control, and the commensurability ruling of `contrib/unbounded_c/` §4(v). Boundary-visible complexity is not bounded by any f(d), and the true tradeoffs bound c only through the conformal onset scale (the C7 dichotomy).

## Status

Heuristic overall. Items 1 and 3 are proved and checked numerically, item 2 is read, items 4 and 5 are an outline with the open steps named. Writing out the computation track would make the cell step a theorem. A counterexample would still need infrared control.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), sub-session of the coordinating session, 2026-10-27.
- AI assistance. gpt-6-astra (high), closed book via `gpt_solve.sh` on `spec_P1_port_idle.md`, output in `astra_P1/`. Literature read by me. [read, abstract] means only the abstract was seen.
- Checked by me. My `port_idle_check.py`, written before seeing Astra's numerics on a different model, agrees. I reran Astra's `experiment.py` (tables reproduced to all digits) and checked its T = 1 port-gate example by hand. The dressed Bell-pair cell and the interface-span argument are mine (`majorana_bank_check.py`, `interface_span_check.py`).

## Argument

### 1. What the literature supplies

- Cubitt–Montanaro–Piddock (1701.05182), restated in [read, 2003.13753 §2.1]. Local encodings obey E(A_j ⊗ 1) = (A'_j ⊗ 1)E(1), so V = ⊗V_i. A (Δ, η, ε)-simulation has ‖Ṽ − V‖ ≤ η, ‖H'_{≤Δ} − Ẽ(H)‖ ≤ ε, and preserves spectra, partition functions and dynamics, but nothing is said about matrix elements reaching above the cutoff. If the port qudit is encoded as itself, F_phys Ṽ = ṼF + O(η), and the boundary resolvent is preserved to O(η/z + ε/z²) (my addition).
- Kohler–Piddock–Bausch–Cubitt (2003.13753) [read, §3.1, §3.3, §3.4, Lemma 4.3, Thms 4.5, 4.6]. A translation-invariant nearest-neighbour chain of constant local dimension simulates any finite target. The description sits in the binary expansion of the chain length (Gottesman–Irani counter) or of a coupling phase read by phase estimation. Input qudits occupy their own track, and the phase-estimation machine acts on them only after idling, so V = Σ_i Φ_idling(i)⟨i| is exactly local with ‖V − V'‖ ≤ 2T_PE/T. Then η = O(T³/Δ + T_PE/T), ε' = ε + O(T⁴/Δ), at the price Δ ≥ Δ'T² + T³/η + T⁴/ε. An explicit Δ dependence of h is unavoidable (Cor. 4.8).
- Gottesman–Irani (0905.2419) [read, abstract]: chain length as the only input of a QMA_EXP-complete 1D problem. Bausch–Cubitt–Ozols (1605.01718) [read, abstract]: the same at local dimension about 40. Watson–Bausch–Gharibian (2012.12717) [read, abstract]: local observables of translation-invariant ground states are P^{QMA_EXP}-hard to estimate. BCLP (1810.01858) and GI's segments: `contrib/unbounded_c/` §1 [read, agent]. Caha–Landau–Nagaj (1712.07395) [read, abstract]: idling clocks.
- Edge codes. A nearest-neighbour parent of an MPS injective on one site has edges of dimension at most √d [recalled, quant-ph/0608197], and bosonic chains without symmetry have no stable edge degeneracy [recalled]. Neither applies, since our codes are fine-tuned, not protected, and frustration-free history Hamiltonians have ground spaces as large as their witness register [read, 2003.13753 Lemma 4.1]. I found no bound on edge degeneracy by d alone. The area law [recalled, 1301.1162] gives g ≤ exp(Õ(log³d/ε)) at relative gap ε.
- Bravyi–DiVincenzo–Loss (1105.0675) [read, abstract]: Schrieffer–Wolff and linked clusters for lattices.

### 2. Port-idle history states

Inputs ℂ^q_p ⊗ ℂ^q_s ⊗ H_w (port, partner, workspace), clock t = 0..T, gates on H_w only, W_t = U_t⋯U_1, and H = H_prop(J_c) + H_init(J_c) + Σ_t |t⟩⟨t| ⊗ L h_t with h_t on (p, s, one workspace qudit c). Since F commutes with every W_t and with the clock, Fη(ψ) = η(Fψ) for every port operator, so the source F_aG never leaves the history space. The first-order target is A = Σ_t W_t†h_tW_t. Astra proves (`astra_P1/RESULT.md` §1 to §5):

- If every h_t annihilates W_tG, the code histories are exact zero modes, and with positive h_t the whole ground space for every J_c. Signed h_t need a clock-scale condition, and a counterexample shows it is necessary.
- ‖R_sim(z) − R_A(z)‖ ≤ δ/((Δ + z)(Δ + z − δ)), δ = ξ²/(γ − M + z), γ = J_c/2L², M = L max‖h_t‖, ξ ≤ M, and similarly for K2. Clock-excited states enter through a Schur complement, not discarded.
- Syndrome-sector targets are reachable by Pauli expansion and Clifford routing to c. Literal universality on a fixed workspace fails, and one initialized ancilla of dimension q² + 1 repairs it. Terms K_{ps} ⊗ Y_c with K|Ω⟩ = 0 keep the code exact under any approximate compilation.
- Every cell has an exact clock tower J_c[1 − cos(kπ/L)], so at fixed ‖h‖ the normalized gap shrinks like 1/L², the price of the construction.

In `port_idle_check.py` (T = 6, two workspace qubits, random gates) the code is exactly two-fold, the port compression is 2 × 10⁻¹⁷, and the errors in R(z) and K2 fall from 0.16 and 0.98 at J_c = 10² to 1.7 × 10⁻⁵ and 9.4 × 10⁻⁵ at J_c = 10⁶. A gate on the port, undone one step later, leaves errors 0.16 and 0.83 at every J_c while the spectrum still converges. Astra's model gives slopes −0.999 and limiting errors 0.75, 0.69 [numerical].

### 3. Two gadgets

(A) Second order (G1, C9 §2). Code |Ω⟩_{p,s} ⊗ |i⟩, enforced frustration-free by Δ(1 − |Ω⟩⟨Ω|) at t = 0, so G†F_aG = 0 and the first-order logical term vanish exactly. Spectral weights come from terms K ⊗ L with K|Ω⟩ = 0, and only they carry error ε_H, entering K2 at relative size ε_H/(ρ²Δ).

(B) First order, dressed Bell pair (new). With |Ω_a⟩ = (F_a ⊗ 1)|Ω⟩, Clifford generators γ_a on the logical register and B = Σ_{a≤r}|Ω_a⟩⟨Ω| ⊗ γ_a, put H_L = Δ(1 − |Ω⟩⟨Ω|) + g(B + B†). Since B†B = r|Ω⟩⟨Ω| ⊗ 1, each logical ψ spans with r^{-1/2}Bψ a two-level block with the ψ-independent matrix [[0, g√r], [g√r, Δ]]. The ground space is exactly 2^n-fold for every g, and

G†F_bG = −(sin 2θ/√r)γ_b + (sin²θ/r) Σ_{a,a'} q⁻¹tr(F_aF_bF_{a'}) γ_aγ_{a'},  tan 2θ = 2g√r/Δ.

The bilinear vanishes when the F_a are odd under a port grading (q₊ + q₋ = q), since closed walks on a bipartite graph have even length, leaving r = 2q₊q₋ channels (C9's graded count). The mirror term with η_a = iγ_aP on a second pair gives C8's joint form. The two ends do not commute, but on the span of {ψ, B_Lψ, B_Rψ, B_LB_Rψ, B_RB_Lψ} every Gram entry and matrix element reduces to sums of γ_a² and η_a², so the joint ground space stays exactly 2^n-fold. `majorana_bank_check.py` finds splittings below 10⁻¹⁴ and compressions exactly αγ_b and αη_b (residual below 10⁻¹⁴, equal α at both ends) for graded cells up to q = 4, r = 8, n = 4. Each term is port-idle (a Clifford W_t maps γ_a to Z_c). The first-order logical bond α²λ Σ_a η_{a,j}γ_{a,j+1} is r uniform Majorana chains, c = r/2 = q₊q₋ (C8 §2), untuned at this order.

### 4. The chain

Site space ℂ⁴ (selector #, L, X, R) ⊗ ℂ^{k_c} (computation track with logical and ancilla qubits) ⊗ ℂ^q (register), with terms

1. The C7 selector (its Eq. 13), cell length b(N) chosen by μ.
2. A standard-form, frustration-free computation track with # as boundary symbol, blind to L, X, R (so it commutes with the walker). It decodes N from b (counter) or from a coupling phase and produces the schedule (W_t, h_t).
3. Registers. On # and on the site before # (the ports) the register is never gated. On the next site inward it is the partner. Elsewhere it is an initialized ancilla, so no free spectators remain.
4. The interface. One menu item, Π_Ω on (port, partner) conditioned on a head letter, plus its mirror. For route A the partner stays static, as in Astra's proofs. Route B rotates it with a finite gate set, which keeps the port idle but loses compilation-proof exactness, so its code splitting must stay below α²λ. Pure-port terms F ⊗ 1 are unreachable (partial traces off the port are invariant) and unneeded (γ_a and T_μ are traceless). `interface_span_check.py` shows that the conjugates (1 ⊗ W†)(Π_Ω ⊗ Y)(1 ⊗ W), W on partner plus one ancilla qubit, span every Hermitian operator with scalar pure-port part (61 = 64 − 3 and 253 = 256 − 3 at q = 2). Without partner rotations the menu would need about q² head letters, and route B would fail its dimension count.
5. The bond λ Σ_a F_a ⊗ F_a on (port_i, port_{i+1}) ⊗ [ℓ_{i+1} = #], commuting with everything above.

Legal cells of any length carry zero computation energy, so selection is unchanged. For route A the hierarchy is J = ρ²λ²/Δ ≪ λ ≪ ρ²Δ, ε_H ≪ ρ²Δ, J_c ≥ 2L²(2M + D_req) (Astra Eq. 13), walker and defect gaps (3π²t/b², ~t/b⁴) above λ, and K > 4t. All are finite for each N and independent of L, though J/‖h‖ is tiny (ρ possibly 3^{-O(N²)}, G1 §8).

d0 = 4k_c q. Route A (q = 2) has d0 = 8k_c for all N, so c = N − 1 > d0 − 1 once N > 8k_c, and c is unbounded, both conditional on §5. Route B with q₊ = 8k_c, q₋ = 8k_c + 1 has c = 64k_c² + 8k_c > d0 − 1 = 64k_c² + 4k_c − 1.

### 5. Where it fails

1. The computation track is not written out, and k_c is a guess. Astra's audit (RESULT §12) lists the obligations. Every illegal sector (no head, several heads, cycles, histories cut by a marker) needs a gap bound. The Bell pair must be an energy term, never a legality constraint, or F leaves the legal sector. Coefficients must come from a finite menu at stated precision. The selection margin (~t/b⁴) must beat the inter-cell shift λ²‖K2‖. Each looks routine [heuristic].
2. Infrared control is the decisive gap, and it has nothing to do with cells. For route A, SU(N)_1 receives marginal SU(N)-breaking current-current terms at O(λ³/Δ²) and O(ε_H), which must land in the marginally irrelevant basin, an open condition. A target without relevant operators (incommensurate Luttinger, n ≥ 9, C8 §6) can be synthesized by C9 §2, but its stability is RG lore. Route B looks better. Each flavour has one Majorana per cell, right movers at k = 0 and left movers at k = π, so every Majorana mass carries momentum π and one-cell translation forbids it. Small translation-invariant quadratic corrections only shift Fermi points, and the interacting problem is again C8 §6's Luttinger liquid [heuristic]. No exact symmetry is needed, which matters because C8 Theorem 3 excludes the large ones.
3. The thermodynamic limit needs a local Schrieffer–Wolff expansion with per-bond errors for decoupled gapped cells, standard-looking but unwritten.
4. Commensurability. Rings with b ∤ L carry frozen defects, as in `contrib/unbounded_c/` §4(v).

### 6. Is there an obstruction theorem?

Not of the form asked. KPBC reproduce any finite gadget on open fixed-d chains with end-site observables preserved to O(η), and the port-idle variant is exact on the code. Astra's elementary example (RESULT §13, checked by me for b ≤ 5) is a translation-invariant d = 4 chain with a unique vacuum, gap at least μ − 2|t| at every length b, zero port compression, and b distinct boundary poles. Boundary-visible complexity is unbounded at fixed d, even at fixed gap.

True statements are tradeoffs. Astra proves two [proved, RESULT §14 and §15, the first checked by me]. Without any locality, ‖K2‖ in operator-Schmidt nuclear norm is at most mN/2Δ, so a swap coefficient obeys |β|Δ ≤ mN/(2(N² − 1)), and the exchange through a fixed port falls at least like 1/N. With a Lieb–Robinson bound, a swap of visibility η_vis = |β|g/‖v‖² needs (1 − ε²)N² ≤ d^{2ℓ} + 1 with ℓ ~ (u/εη_vis g) log(1/εη_vis). Others are the area law, KPBC's growth of Δ with precision, the clock tower of §2, and exponential clustering [recalled, Hastings–Koma], which puts C8's code (non-commuting left and right compressions) only in cells with gap below about (v/b) log(1/α²). Our cells satisfy all of these because their gaps and logical scale shrink with N. A tradeoff bounds c only together with a lower bound on the relative energy scale of the conformal window, the onset-length problem of C7 Theorem 3. The dichotomy of REPORT §9 survives, with its second branch now hinging on infrared control rather than on cells.

## Reproducing

- `port_idle_check.py` (numpy, seconds): idle versus touched port, output `port_idle_check_stdout.txt`.
- `majorana_bank_check.py` (numpy, scipy, six minutes, mostly the q = 4 joint cell), output `majorana_bank_check_stdout.txt`. Use OPENBLAS_NUM_THREADS=1 (multithreaded LAPACK failed to converge here).
- `interface_span_check.py` (numpy, under a minute).
- `spec_P1_port_idle.md` is the Astra spec. `astra_P1/` holds RESULT.md, `experiment.py`, `checks.py` (seconds) and outputs, and `astra_P1.out`/`.err` are the raw run.

## Relation to the rest of the repository

This answers the missing programmable cell of C9 §5 and REPORT §9.1 to 9.3. Cell realization is no longer an obstruction (modulo writing out the computation track), and infrared control is the open step for a c–d counterexample, with commensurability awaiting Jonah's ruling. "Where things stand" could say that both block routes now reduce to infrared control of an explicitly constructible logical chain.
