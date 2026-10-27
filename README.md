# The c–d conjecture

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242), J. Stat. Mech. (2024) 113103) conjecture that a critical, translation-invariant, nearest-neighbour Hermitian chain on (ℂ^d)^⊗L flows to unitary CFT(s) with total central charge c ≤ d − 1. The bound is saturated by the SU(d) Uimin–Lai–Sutherland chain h = P (the swap), whose IR is SU(d)₁ with c = d − 1, so the conjecture is equivalent to c_max(d) = d − 1. This repository collects what we know about proving or refuting it. It is meant as a shared working space that collaborators push to and pull from.

Maintained by Jonah Kudler-Flam.

## Where things stand

The conjecture is open. We know of no proof and no counterexample within its scope, and it appears to be true.

- Entropic arguments cannot prove it. A monotone lattice c-function bounded by K log d at the lattice scale would give c ≤ K log d, which SU(d)₁ violates at large d. Thermodynamic, Lieb–Robinson and sum-rule bounds only control c times a power of the conformal window.
- Proved in restricted classes: quasi-free nearest-neighbour chains of any period, where c ≤ log₂ d; SU(d)-symmetric bonds; Gaussian IR theories whose modes are resolved by commuting on-site charges, where c equals the rank of their susceptibility matrix and is at most d − 1; WZW IR theories obeying a level ceiling.
- Interacting invariants bounded by d − 1 exist (the affine dimension of a parton occupation code, the susceptibility rank), but they equal c only when the IR has no neutral sector. Every route so far stops at the same missing lemma, a locality-based bound on neutral (Majorana, parafermion, coset) stress-tensor sectors.
- Numerical and analytic searches found nothing above d − 1. These covered all nearest-neighbour qubit chains, Sp(4), spin-3/2 SU(2), spin-orbital and spin-1 families, and integrable and non-compact candidates. The non-compact candidates turned out to be non-Hermitian.

The full account, with what was checked and how, is in [REPORT.md](REPORT.md).

## Open problems

1. A locality-based budget for neutral sectors: from a fixed two-site bond on ℂ^d ⊗ ℂ^d, bound all emergent stress-tensor sectors by d − 1.
2. An improvement-independent version of the velocity-free stress response F_h(q) = ∫ S_h(q,ω)/ω² dω → C|q|/(24π), and a microscopic bound on it.
3. A nearest-neighbour level ceiling k ≤ ind(R) for chains with an on-site symmetry flowing to WZW models, or a chain that violates it (for example Sp(4)₂ on ℂ⁴, or SO(N)₃ on an SO(N) vector).
4. Larger-size numerics on the two loose ends: qubit chains with entropy slopes 1.35 to 1.7, and small-angle Sp(4) (see REPORT.md §7). Cluster-ready scripts are in `codex/work_C4_D2/` and `codex/work_C4_D3/`.
5. Bethe-ansatz completeness for fundamental gl(d) chains together with control of massless strings, which would give c ≤ d − 1 for that integrable class.

## Layout

- `REPORT.md` is the working report: statement, literature status, proof routes, verdict, and the results of every model run, with provenance.
- `lit.md` is the literature map, with each claim tagged [read] or [recalled].
- `routes_draft.md` holds the short derivations behind REPORT.md §3.
- `specs/` contains the self-contained problem statements given to the model runs (`_background.md` is shared by the campaign specs).
- `codex/` contains the outputs of those runs (`*.out`), with code, data and full write-ups in `codex/work_*/`. `codex/opus_check_C3/` holds independent check scripts.

## Provenance

This material was produced with AI assistance and has not been refereed. Claude (Anthropic, Opus 5.5) planned the work, wrote REPORT.md and checked results, so "I" in REPORT.md refers to Claude. Closed-book derivations and numerics were run on GPT through the Codex CLI. Almost all used gpt-6-astra, and two early runs on gpt-5.6-sol are kept for comparison. Claims are tagged by how they were checked, and contributions should keep that convention.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).
