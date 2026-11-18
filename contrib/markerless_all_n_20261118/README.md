# An arbitrary-n band mechanism at physical d=9

2026-11-18. **[Proved finite-dimensional statements; internal algebra and
numerical checks; no independent human review.]**

## Claim

For every finite n>=1, an explicitly modified version of the markerless
nine-state express-rail bond has a leading quadratic logical symbol with
exactly 2n simple Fermi points and strictly different sums for all different
unordered pairs of points, including repeated-point pairs. The microscopic
Hamiltonian is Hermitian, nearest-neighbour and one-site translation
invariant, with nine physical states and coefficients fixed as L grows.

The modification gives the vacuum head a boundary bias and hopping q,
the active data rail a hopping w, and uses a physical L1-only field. No
marker or other state is added. The source is X_n, the inactive endpoint is
D_n, and the actual carrier has 2n+2 vertices. This proves an analytic
arbitrary-n **logical-band family**, not the original uniform-field schedule
at all n and not a microscopic unbounded-central-charge theorem.

The central new theorem is **L_j<0 for every n>=2 and 0<=j<=n-2** in the
rational small-loop moment criterion left open by the eleven-state analysis.
In fact every summand is strictly negative. The proof uses ordinary
orthogonal polynomials after the substitution x=t^(-2), a signed resolvent
identity, and positive-zero interlacing. It does not extrapolate finite
root scans.

## Results and status

- **[Proved]** [All-n moment sign theorem](MOMENT_SIGN.md), including a
  no-cancellation formula for the derivative form.
- **[Proved]** [Literal nine-letter controls, clock and readout](DERIVATION.md):
  exact legal-sector invariance, an exactly degenerate zero-field code,
  clock gap >=1/2, and a zero-port length selector including head-number
  competitors and zero-data cells.
- **[Proved]** [Markerless small-loop response](DERIVATION.md): the express
  delay is 2(n+1)/D; the loop variable is eta=t_d D^n; normalized energy is
  D^(2n+2). The mixed pair response has coefficient -16(n+1), rather than
  the eleven-state value -16(n+2).
- **[Proved]** A finite parameter schedule for arbitrary n, with 2n simple
  roots and separation of every pair sum in the leading quadratic symbol.
- **[Numerical / exact arithmetic]** Literal two-cell Fock reconstruction,
  exact rational moment identities through n=6, finite-delay asymptotics at
  n=2,3, and independent 70-digit roots at n=2. These check the derivation;
  they are not interval certificates for all physical corrections.
- **[Open]** Nonzero-port thermodynamic sector selection, uniform interacting
  matching and the unitary compact Virasoro CFT. No new value of microscopic
  c or c_eff is established.

## Verification and reproduction

Author and AI assistance: OpenAI Codex, user-requested GPT-6 Astra / Max.
No additional agents were used. The exact sign proof and its finite
arithmetic implementation were developed here. The scalar resolvent method
is credited to the previous eleven-state loop analysis and the original
markerless source; the literal new bond and direct Fock checks were rebuilt
for this contribution. There has been no independent human check.

From the repository root, in a Python environment with numpy, scipy and
mpmath:

```sh
python contrib/markerless_all_n_20261118/checks.py
```

The script forces one numerical worker and writes [checks.json](checks.json).
This run used `/tmp/cd-small-d-venv/bin/python`, took about 0.43 seconds,
and launched no background job or DMRG. It checks 6,642 local transition
contexts, the complete n=1,2 two-cell Fock coefficient (400 and 3,600 literal
states), and the physical clock/readout. Largest relative hopping error was
6.9e-16; largest full Fock reconstruction error was 1.4e-10, limited by
subtracting large scalar energies in double precision. The finite-delay
checks converge to the predicted coefficient. At the saved n=2 high-precision
point the full pair-sum margin is about 2.98e-7; it is a validation point,
not an interval-certified microscopic member.

## Relation to the project

This addresses focused-round assignment 3. It replaces the unresolved
all-n moment condition by a proof and transfers the resulting band
mechanism to the smallest currently sourced express-rail alphabet, d=9.
The existing original certificates remain exactly n=9,10,11,12,14,16,20.
The published source in `contrib/markerless_d9_d10/` is present; older notes
calling it unlocated are superseded. This is not the separate unstable
single-rail nine-state seed.

[SCOPE_AND_SIMPLICITY.md](SCOPE_AND_SIMPLICITY.md) separates finite band
algebra from CFT claims and counts the physical controls and literal terms.
The existence proof does not improve the severe coefficient hierarchy.
[PROPOSED_SHARED_SUMMARY.md](PROPOSED_SHARED_SUMMARY.md) contains proposed
coordinator text; no shared root or coordination file was edited.
