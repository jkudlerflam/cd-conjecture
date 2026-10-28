# Prior work on unbounded c at fixed d

## Claim

This folder is a literature check, not a result. The claim under development (REPORT.md §9) is that for critical, translation-invariant (TI), nearest-neighbour (NN), Hermitian chains on (ℂ^d)^⊗L the Virasoro central charge is not bounded by any function of d. At a fixed d₀ (roughly 10² to 10⁴), a family of bonds h_N would make the sites form long cells spontaneously. Each cell runs a port-idle history-state computation that prepares a logical code, and neighbouring cells couple through one small bond to a CFT whose c grows with N. The question is whether the literature contains this, contradicts it, or supplies its pieces. Items already in lit.md or in the literature sections of other `contrib/` folders are not repeated.

## Status

Survey, finished 2026-10-27. The claim appears to be new, and nothing I found contradicts it. A read by someone from the Hamiltonian-complexity community would upgrade this.

## How it was established and checked

- Author: Claude literature agent (Anthropic, Opus 5.5), at Jonah's request.
- AI assistance: all searching and reading. No Codex run, since nothing here is a derivation.
- Checked by: nobody yet. Sources were the INSPIRE, Semantic Scholar and arXiv APIs and web search. [read] means I opened the abstract or the section cited, [recalled] means memory or a search snippet, and [own] is my remark.

## Findings

### 1. Large c at fixed d, and simulating CFTs

No paper I found claims that c, the entanglement log coefficient, or the number of gapless modes of a TI NN chain can be arbitrarily large at fixed d. None simulates a CFT, Luttinger liquid or WZW model with a TI universal Hamiltonian, and the universality papers list this as open. Kohler, Piddock, Bausch and Cubitt end by asking for universal models whose boundary theory has a discrete conformal symmetry [read, 2003.13753, Discussion]. Zhou and Aharonov tabulate every TI universal construction, all with interaction energies exp(poly(n)), note that the 1D one changes its interaction with the target, and leave strongly universal fully TI simulators open in every dimension [read, 2102.02991, Table I and Discussion, now in Nat. Commun. 2026]. Piddock and Bausch sketch simulation of a TI target by repeating the coefficients periodically, but it is still one global computation writing one history state on a finite lattice [read, 2001.08050 §4, Lemmas 15–16, Cor. 17]. Every TI simulation result I found is therefore finite-size with a global clock, and none gives a fixed h whose infinite-chain IR is the target. The proposal fills that gap with many local, repeated computations.

Three complexity results bear on the cells. González-Guillén and Cubitt prove that every history-state family has a gap closing at least as fast as O(1/n) [read, 1810.06528, abstract and §1]. At fixed N this is harmless, but the cell gap shrinks as the code grows, so the logical bandwidth and the conformal onset scale must fall with N. Aharonov and Irani's finite TI chains have ground energy 1 − cos(π/(L+1)) set by the computation length, a quantum-walk energy of the C7 selector's form, while their thermodynamic-limit construction (like Watson–Cubitt's) is 2D, with lengths imposed by a Robinson tiling rather than chosen spontaneously [read, 2107.06201 §1.3–1.4; 2107.05060, abstract]. Finally, finite sizes need not predict the thermodynamic limit. In 2D the threshold size can grow faster than any computable function of d [read, 1512.05687, abstract], phase diagrams can be uncomputable [read, 1910.01631, abstract], and RG flows can track one fixed point for an uncomputable number of steps before leaving it [read, 2102.05145, abstract]. [own] This bounds nothing, but a proof of c ≤ f(d) through the onset length K(d) of §9 would have to control a quantity these constructions make arbitrarily large.

### 2. Bounds at fixed d

Nothing new since lit.md. No paper bounds c, the log coefficient, or the number of gapless modes by a function of d alone. The only bounds of this kind are for tensor networks, where c is limited by a bond dimension that the Hamiltonian does not fix (§5 below). In restricted classes, the SU(N) U → ∞ Hubbard and supersymmetric t–J chains are SU(N) Luttinger liquids [recalled, Kawakami, PRB 47, 2928 (1993)], so [own] c = 1 + (N − 1) = N at d = N + 1. This saturating family with mobile carriers is missing from lit.md.

### 3. Cell formation and large internal spaces

Long periods. Classical NN chains have ground-state period at most d [recalled]. The known devil's staircases need long-range interactions or two dimensions. Papanikolaou, Raman and Fradkin get arbitrarily long periods from a short-range 2D quantum dimer model [read, cond-mat/0611390, abstract], while the ANNNI chain has only periods up to 4 at T = 0 [recalled]. The only 1D precedents for long cells at fixed d are Hamiltonian-complexity constructions (Gottesman–Irani, BCLP, Aharonov–Irani), which select lengths through a quantum-walk ground energy. The C7 selector seems to be a low-dimension instance of that mechanism.

[own, standard bosonization] There's a simple reason 1D has no such staircase. Mobile walls at density 1/p form a Luttinger liquid, and lock-in to period p needs p²K < 2, which short-range wall interactions at low density cannot reach (K → 1). The walls float, which is the "length fluid" of §9.5. A locked selector must therefore set the cell length through internal energetics, as C7 does.

Emergent large spaces. The quantum Wielandt inequality bounds the injectivity length of a TI MPS by O(D² log D) in the bond dimension D [read, 1809.04387, abstract], improving the O(D⁴) of Sanz et al., who also derive a dichotomy for frustration-free MPS ground states [read, 0909.5347, abstract]. [own] This is the closest published analogue of C8's finite-witness theorem K_FF(q) (§9.7) and should be cited beside it.

Mobile carriers. Squeezed space is how a Luttinger liquid carries an internal space on mobile objects [recalled, Ogata–Shiba 1990; read, cond-mat/0312052, abstract], and the SU(N) liquid above shows the price, since one carrier's internal dimension is at most d − 1. The proposal escapes this because a cell spreads N > d states over many sites. The spin-incoherent Luttinger liquid describes the resulting hierarchy. When the internal exchange is exponentially small next to the charge energy, there is a wide window where only the charge mode is coherent [read, cond-mat/0611597, abstract]. [own] Any probe above the logical bandwidth, finite-size numerics included, would then see c = 1 from the length fluid, and c = N only below it.

### 4. Citations of Latorre–Sierra

INSPIRE (refersto:recid:2771764) still returns the same 7 records on 2026-10-27, and Semantic Scholar returns 3 of them. Web and arXiv searches find no 2026 paper claiming a proof, a counterexample or a refinement, and no follow-up by the authors [read, API queries].

### 5. c per site and c per bond

The tensor-network literature bounds c per bond, not per site. For a MERA with rescaling factor k, Bao et al. derive c ≤ 12 (k − 1)/(k ln k) ln χ, which is c ≤ 6 log₂ χ for a binary MERA, and c ≤ 12 (k − 1)/ln k · ln χ without their ancilla assumption [read, 1504.06632 §IV.C, Eqs. (19)–(21)]. Jahn, Zimborás and Eisert find the maximal central charge of boundary states of regular hyperbolic networks, a tiling-dependent multiple of ln χ, and saturate it with Majorana dimers [read, 1911.03485 §III, Eqs. (14), (21), (23)]. MPS finite-entanglement scaling, ξ_χ ∝ χ^κ, is likewise governed by c [read, abstracts of 0812.2903 and 0712.1976]. [own] None of this touches c–d, since nothing ties χ to d, but it shows χ must exceed d for large-c ground states. ULS has c = d − 1, which beats 6 log₂ d for d ≥ 31 and 12 log₂ d for d ≥ 76, so a scale-invariant binary MERA with χ = d cannot represent it there. A "maximal c per site" for Hamiltonians appears only in Latorre–Sierra.

## Novelty assessment

As far as I can tell the claim is new. Its pieces are known separately, namely quantum-walk length selection, idling history states, second-order effective couplings, stable incommensurate multi-flavour liquids and squeezed-space carriers. Combining them so that a fixed-d TI chain has a large-c CFT as its thermodynamic-limit IR is new, and is the kind of result the universality papers list as open. Nothing contradicts it. The only bounds are entropic (c ≲ onset length × log d, in the repository) and per bond (c ≲ ln χ).

A referee from that community would likely raise four points. (i) Gottesman–Irani and BCLP give decoupled segments, so the burden is coherent inter-cell coupling at finite coupling. (ii) The cell gap closes with N, so the hierarchy of scales and the onset length must be stated. (iii) The spin-incoherent window will hide c = N in finite-size numerics. (iv) h_N depends on N with large coupling norms, as in every TI universal construction. c–d allows this, but it should be said up front.

## Reproducing

The citation count comes from `https://inspirehep.net/api/literature?q=refersto%3Arecid%3A2771764`. PDFs were read with `pdftotext` in a scratch directory and are not committed.

## Relation to the rest of the repository

This is the literature half of the adversarial review listed in REPORT.md §9. It does not change "Where things stand". Two small follow-ups are to add the SU(N) Luttinger-liquid family to lit.md and to cite quantum Wielandt beside §9.7.
