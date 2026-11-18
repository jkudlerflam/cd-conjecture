# Uniform local matching and frozen sector selection at d=9

2026-11-18. **[Proved local and sector bounds; numerical coefficient
checks; no independent review.]** This is stability-round assignment C.

## Claim

For the exact modified nine-state Hamiltonian of
[markerless_all_n_20261118](../markerless_all_n_20261118/README.md), the
elimination of clock and charged excitations has a connected fourth-order
local normal form with an explicitly bounded remainder independent of
the number of cells. All nonquadratic operators and the nonconserved L1
readout are retained. The model remains Hermitian, microscopically
nearest-neighbour and one-site translation invariant, at physical d=9,
with every coefficient fixed as physical length L grows.

There is also an exact nonzero-port ground-sector theorem on compatible
rings: with gamma=E_(n+1), B=sqrt(2)|lambda| and gamma>5B, each
nonselected boundary cell costs at least gamma-5B relative to the exact
dressed selected ground energy. This includes different head counts,
zero-data cells and inert forbidden/headless cells, and permits K=0.

## New results

- **[Proved]** The order-lambda^2 g^4 coefficient is at most quartic in
  logical fermions, with full-Fock norm bound zeta/32. Its complete finite-g
  tail is bounded locally, not only in the one-particle block.
- **[Numerical]** The first quartic coefficient is nonzero on literal n=1,2
  two-cell clusters; full tensors and Fock reconstruction checks are saved.
- **[Proved]** Clock distortion from the L1-only readout generates
  -2 h_L^2 C_ij n_i n_j, with C_ij>0 for every pair and every n>=2. A
  one-dimensional current formula gives C exactly. Endpoint visibility
  distortion is included explicitly.
- **[Proved / numerical checks]** The complete Hermitian fourth-order
  expression includes connected three-cell port processes, field-port
  terms and field-square terms. Independent wave-operator recursion and
  direct finite-lambda diagonalization check the normalization.
- **[Proved]** A full local interaction norm, an explicit finite-order
  majorant and an exponentially weighted remainder bound uniform in cell
  count. The unprojected remainder is O(lambda^5); its code projection is
  even and starts at O(lambda^6). These are not exact block-diagonalization
  or infinite-order convergence statements.
- **[Proved]** Exact frozen ground-sector selection for nonzero ports on
  L divisible by n+1, using the exact dressed bulk energy and complete
  two-site-boundary cell decomposition.

The finite coupling schedule can make these bounds small compared with
the smallest physical mode velocity and all leading-band momentum margins.
It remains severely hierarchical. An all-scale interacting matching/RG
theorem is still missing, so **no microscopic unitary compact Virasoro
central charge or unbounded-c CFT theorem is claimed**.

## Details

- [COEFFICIENTS.md](COEFFICIENTS.md): connected coefficient formulas,
  actual readout distortion, fermionic normalization and literal examples.
- [UNIFORM.md](UNIFORM.md): the local interaction norm, constructive
  normal form, explicit constants and global spectral-overlap limitation.
- [DEFECTS.md](DEFECTS.md): exact nonzero-port frozen sector selection,
  including wrong head counts and zero-data cells.
- [SCALES_AND_SCOPE.md](SCALES_AND_SCOPE.md): physical velocity and
  momentum units, forward scattering, an explicit smallness schedule and
  the remaining infrared hypothesis.

## Verification and reproduction

Author and AI assistance: OpenAI Codex, user-requested GPT-6 Astra / Max.
No additional agent or chat was used. The literal model is imported
read-only from the preceding contribution. The connected coefficients are
checked by a separate Bloch wave-operator recursion and direct low-energy
diagonalization in the positive-overlap frame. There has been no independent
human review. Local SW background was read in Bravyi, DiVincenzo and Loss;
the explicit majorants and applications are derived in this folder.

With numpy and scipy, from the repository root:

```sh
python contrib/markerless_uniform_matching_20261118/checks.py
```

The run used `/tmp/cd-small-d-venv/bin/python`, one numerical worker, and
writes [checks.json](checks.json). It took about 0.22 seconds. The largest
literal cluster has 3,136 states in an exact even-token sector of three
n=1 cells. The conversion checks use all logical Fock sectors on n=1,2
two-cell clusters; the sector check enumerates every physical word on
L=2,3 rings. A saved n=2 band point supplies the velocity audit without
a new root search. No DMRG, large scan or background production was run.

The example parameters make coefficients measurable and are explicitly
diagnostic, not newly certified all-n band members. Full matrices are
saved rather than discarding unresolved small higher-body entries.

## Relation to the project

This advances the modified frozen d=9 family beyond finite-cluster error
budgets and closes its nonzero-port compatible-ring sector comparison.
It does not transfer a conserved-charge field argument from the distinct
original d=10 mobile model. That model's vacancy response and selector
remain separate. The existing incompatible-ring issue of frozen sectors
also remains.

The next decisive calculation is an all-scale, volume-independent
integration/RG argument controlling the residual auxiliary coupling,
marginal forward scattering and dressed momentum cutoffs. The present
finite-order normal form is a quantitative input to that calculation.
Proposed coordinator text is in [PROPOSED_SHARED_SUMMARY.md](PROPOSED_SHARED_SUMMARY.md).
