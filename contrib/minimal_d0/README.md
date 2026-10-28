# How small can d0 be? Alphabet accounting for the cell construction

## Claim

Class: one-site translation-invariant, nearest-neighbour, Hermitian chains on (ℂ^{d0})^⊗L with L-independent couplings, built as in `contrib/programmable_cell/` (spontaneous cells, port-idle history states, first-order route B or second-order route A). Question: the smallest d0 at which the construction gives c > d0 − 1, or c unbounded.

1. The estimate d0 = 4k_c q overcounts. Registers need not be tensored onto every letter, only two letters need the q-dimensional port register, and the cell program can run on the smallest known translation-invariant history machine (42 dimensions).
2. One static Bell-check term on the (port, partner) bond is the only term that has to touch the port. Every port-touching operator a cell needs is a positive time-average of it, conjugated by gates on the partner and its neighbour, plus terms away from the port [proved, Astra; checked]. No head visits a port-type site.
3. Route A (unbounded c, qubit ports) needs d0 ≈ 130 with the C7 walker tensored onto a BCO-class machine, and d0 ≈ 45 if the machine's clock selects the cell length. Route B (c = q₊q₋) needs d0 = 2q + R, R the bulk alphabet, about 76 to 88 with a merged selector and 174 to 208 with a tensored walker.
4. The swing factor is the selector. Merging it with the history clock needs a machine whose legal history is the longest penalty-free path in every cell, with no free internal lengths. This sharpens the obstruction of C10 §9.

## Status

Heuristic. Item 2 is proved and checked. The counts rest on the assumptions in §5. One machine written out with the longest-path property of §4 would make the merged counts a construction, infrared control aside.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), 2026-10-27.
- Literature: a Claude Sonnet 5.5 agent extracted the papers of §1 from the arXiv texts.
- gpt-6-astra (high), closed book, on `spec_frames.md` (the span and positivity question of §2). Output in `astra_frames/`.
- Checked by me: `frame_checks.py`, written before Astra's output, reproduces the span dimensions at q = 2 and 4, synthesizes the q = 2 route-B cell from one static term exactly, and checks that partner-Pauli junk keeps the q = 4 code exactly degenerate. `replay_astra.py` replays Astra's recipes. I read Astra's proofs of span and positivity but did not rederive the general-m gate lists. C10 (`codex/work_C10/`) was read, not duplicated.

## Argument

### 1. Literature: the smallest history-state alphabets

| Construction | d | Alphabet | TI | Input |
|---|---|---|---|---|
| AGIK, 0705.4077 | 12 (13 for a rigorous gap) | 4 qubit letters + 4 one-dimensional | no | gates in position-dependent terms |
| Hallgren–Nagaj–Narayanaswami, 1312.1469 | 8 | 2 qubit + 4 one-dimensional | no | same |
| Nagaj–Wocjan, 0802.0886 | 10 | 5 ⊗ 2 | yes | program in the initial state (dynamics only) |
| Gottesman–Irani, 0905.2419 | not stated (BCO estimate ~10⁶) | six tensored tracks | yes | N in binary, counter |
| Bausch–Cubitt–Ozols, 1605.01718 | 42 | 3 qubit + 36 one-dimensional | yes | length → base-6 program by counter |
| Watson–Bausch–Gharibian, 2012.12717 | 44 | BCO + 2 | yes | as BCO |
| BCLP 1810.01858; KPBC 2003.13753 | "very large", unoptimized | tracks + markers | yes | coupling phase via phase estimation, or length |

Sources: [read, AGIK §2, §4], [read, HNN §1, §3–6], [read, NW abstract, §I–II], [read, GI via BCO §2.1, §9], [read, BCO §7.2, §8.1, §8.2, §8.7], [read, WBG Cor. 27], [read, BCLP §V, §VII; KPBC §3.3–3.4, §6]. AGIK and HNN footnotes say a 13th state is needed for a 1/poly gap [read]. Caha–Landau–Nagaj (1712.07395) show idling clocks keep an Ω(N⁻²) gap [read]. No reduction below 42 (TI) or 8 (non-TI) after 2017 was found [searched].

Direct sums, with a qubit attached only to letters that hold data, account for AGIK's 12 → HNN's 8 and BCO's 48 → 39 letters. Length input (a counter) is cheaper than coupling input (phase estimation), and the selector already turns the coupling μ into a length. No paper gives an explicit small d for computations confined to segments.

### 2. One port-touching term suffices

Let p (port) and s (partner) be q-dimensional with Bell state |Ω⟩ and Π = |Ω⟩⟨Ω|, and let "rest" be the cell's other registers. A static term Δ(1 − Π) on the (p, s) bond acts at every clock time, so the history compression turns it into Σ_t f_t W_t†(1 − Π)W_t, with f_t ≥ 0 and W_t the cumulative circuit, which never acts on p.

The real span of {W†ΠW : W unitary on s ⊗ rest} is the set of Hermitian operators with scalar pure-port and pure-rest parts, of dimension (q² − 1)(4q² − 1) + 1 for one rest qubit [proved, Astra; ranks 46 and 946 at q = 2, 4 reproduced in `frame_checks.py`]. Positive weights cost only an identity offset [proved, Astra]. The route-B coupling g(B + B†) and route A's spectral terms lie in this span once their pure-rest parts, ordinary terms away from the port, are removed. At q = 2 the four Bell frames of (p, c), reached by a Pauli I, Z, X, Y on c followed by SWAP(s, c), with weights f, f, f₁, f₂ and 2f = f₁ + f₂, plus g(X_sX_c − Y_sY_c)/2, reproduce the dressed-Bell cell [checked: residual 6 × 10⁻¹⁶, code two-fold, compressions αX_c and αY_c with the predicted α = −0.511]. For q = 2^m Astra gives explicit recipes with q²(q² − 2)/2 frames per channel pair, each at most 13m + 2 gates from {H, S, Paulis, CNOT(s_l, c), SWAP(s_1, c)}, all acting on the (partner, neighbour) bond [proved, Astra; my replay with independent matrices gives residuals below 9 × 10⁻¹⁶ for every Q at m = 1, 2].

Pure Pauli frames on s pull the check back to Bell-diagonal (p, s) terms. Each logical state still spans a star block with a state-independent matrix, so the code stays exactly degenerate and only α becomes flavour-dependent [checked at q = 4, r = 8: sixteen-fold degeneracy, compressions α_aγ_a to 2 × 10⁻¹⁵]. Gates on s alone cannot relabel the channel, though. To first order in a rotation angle the induced ⟨Ω|·|Ω_a⟩ block is anti-Hermitian (F_a commutes with the generator), giving zero port compression [proved, mine], and Astra proves that no s-only conjugation of a realized coupling changes its channel, since tr_s is invariant. Labels need entangling frames under a program that knows the label (Astra conjugates every gate word by a controlled-Q Clifford on s), so route B needs q-dependent loops.

The interface is a port letter # and a partner letter S, each carrying only its q-dimensional register, with the static check on (#, S). Gates on S come from transitions on (S, w) in which only the bulk letter w changes. Cell boundaries read S # # S, with the inter-cell bond λΣF_a ⊗ F_a on (#, #). For route A (q = 2) the partner is an ordinary bulk qubit.

### 3. Program alphabets

Write d0 = Σ over letters of the dimension of the attached register.

Route B. d0 = 2q + R and c = ⌊q/2⌋⌈q/2⌉, so a violation needs q ≳ 4 + 2√(R + 3). The program is uniform. A conveyor of fSWAP sweeps brings each logical mode next to S, labelled frames couple it to its channel, and the right-end event is local because iγ_aP of the mode on the last site is a single-site Pauli. An inventory (mode letters with a phase bit, about 8 qubit-carrying and 8 empty head states, label counters) gives R ≈ 40 to 60 [heuristic].

| R | q | c | d0 = 2q + R |
|---|---|---|---|
| 30 | 16 | 64 | 62 |
| 40 | 18 | 81 | 76 |
| 50 | 19 | 90 | 88 |
| 60 | 20 | 100 | 100 |

Route A. q = 2, so every register is at most a qubit and d0 ≈ R_A. The G1/C9 gadget depends on N without structure I could make uniform, so a universal machine is needed. BCO's alphabet, with its two brackets replaced by the one 2-dimensional letter # (orientation tells the ends apart) and 2 to 4 event head states added, gives R_A ≈ 45, and its base-6 counter already decodes N from the cell length.

On the same machine route A is cheaper than route B by about 2q ≈ 36 and gives unbounded c. Route B wins only if its uniform program runs on a much smaller machine.

### 4. The selector

Tensored walker. C7's {L, X, R} on every bulk letter, # outside, so d0 → (port letters) + 3R. The computation track is frustration-free and blind to the walker, so C7's argument and the machine's illegal-sector bounds (BCO's 1/|M|³, or C10 §8) hold separately.

Merged clock. Run the history clock in hopping form, −JΣ(U|t+1⟩⟨t| + h.c.), so a cell whose legal history has K(b) steps gets −2J cos(π/(K(b) + 2)), and C7's convexity argument selects b with no walker letters. But in hopping form a penalty at one node barely raises a long path. Every locally legal configuration must therefore lie on a path no longer than the legal one, and internal lengths (number of modes, counter sizes) must be rigid or derived from b, or the cell picks whichever split has the longest path. Monotone counters started at zero, domain-wall letters for a unique head, and oblivious control flow give this by design. BCO's machine does not have it as written, since its illegal tapes are removed by penalties with a 1/|M|³ gap that a hopping clock does not respect.

### 5. Counts and what they rest on

| Route | Selector | d0 | c | Assumptions beyond proved inputs |
|---|---|---|---|---|
| A | tensored | ≈ 2 + 3 × 44 ≈ 130 | N − 1, unbounded | BCO-class machine adapts to marker-bounded cells with an idle port (a few letters) |
| A | merged | ≈ 45 | unbounded | a universal machine of about BCO's size with the longest-path property |
| B | merged | 76 to 88 | 81 to 90 | special-purpose machine with R ≈ 40 to 50 and the longest-path property |
| B | tensored | 174 to 208 | 182 to 210 | R ≈ 40 to 50, frustration-free |

The earlier estimates were 8k_c ≈ 10³ to 10⁴ (A) and about 64k_c² (B). All rows inherit the open infrared step (`contrib/ir_control/`, C11) and C10's finite-J caveat. Selected-time energies make the code exact only as J → ∞, with O(1/J) port errors, so route A needs the residual first-order compression below the λ² exchange. A violation at d0 propagates to larger d by REPORT §9.4 [heuristic].

### 6. Open directions

(i) A history machine with the longest-path property, written out with its legality proof, would take the merged rows from estimate to construction. (ii) A uniform route-A gadget, a second-order target whose cell program depends on N only through loop lengths, would remove the universal machine and put unbounded c at the scale of HNN-type alphabets, perhaps d0 ≈ 20 to 30 [speculative]. (iii) A cheaper channel relabelling than controlled-Q Cliffords would shorten route B's counters.

## Reproducing

- `frame_checks.py` (numpy, about a minute; set OPENBLAS_NUM_THREADS=1), output `frame_checks_stdout.txt`. Parts (a) span ranks at q = 2, 4, (b) the exact q = 2 synthesis, (c) the q = 4 degeneracy with Bell-diagonal junk.
- `replay_astra.py` (numpy, seconds) rebuilds every m = 1, 2 recipe in `astra_frames/explicit_recipes_m1_m2.json` with independent gate matrices, output `replay_astra_stdout.txt`.
- `spec_frames.md` is the Astra spec. `astra_frames/` holds RESULT.md, `verify.py`, `recipes.py` (about 10 seconds each) and their JSON outputs, and `astra_frames.out`/`.err` are the raw run.

## Relation to the rest of the repository

This addresses the "fixing d₀" step of "Where things stand" and REPORT §9.6 item 4. The status line could replace "roughly 10³ to 10⁴" by "about 130 with a tensored selector, possibly 45 to 90 if the history clock doubles as the selector".
