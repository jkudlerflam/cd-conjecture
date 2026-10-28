# Red team: the unbounded-c block construction

## Claim

The construction under review (REPORT §9.1 to §9.8, `contrib/unbounded_c/`, `contrib/programmable_cell/`, `codex/work_C7/` Theorem 4, `codex/work_C8/`, `codex/work_C9/`) claims, conditionally, that one-site translation-invariant, nearest-neighbour, Hermitian chains on (ℂ^{d0})^{⊗L} with periodic boundary conditions can have a low-energy spectrum on rings that is a unitary CFT of Virasoro central charge c with

- route A: c unbounded at fixed d0 = 8k_c (k_c the alphabet of a segment-confined computation track, guessed at 10² to 10³), by second-order exchange (G1, C9 §2) between spontaneously selected, port-idle history-state cells with qubit ports;
- route B: c = q₊q₋ > d0 − 1 at d0 = 4k_c q, by first-order Majorana hopping between dressed Bell-pair cells.

Route A would refute every bound c ≤ f(d), and both would refute c ≤ d − 1.

## Status

Refereed, heuristic. The finite-dimensional ingredients survive every recheck. The chain-level argument has one gap that is fatal under the strict reading of the conjecture (frozen defects on rings with b ∤ L), one new problem found here (route B is not stable "untuned", Section 2), one target that should be abandoned (SU(N)₁ in route A), and the two known unwritten steps (the computation track, and a thermodynamic Schrieffer–Wolff reduction). None of these is a proof that the construction cannot work. What would upgrade it: a ruling on the ring-length reading, the computation track written out, and DMRG on the second-order logical chain at the intended parameters.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), red-team sub-session, 2026-10-27.
- AI assistance: none beyond this session. No Codex run, since every check here is a short script whose specification would be as long as the code.
- Checked by: me, by reimplementing (not rerunning) the decisive identities, and by a new second-order calculation for route B validated against exact diagonalization. Literature is tagged [read] (abstract or search summary seen) or [recalled].

## Argument

### 1. What survives

- G1 (Eq. 20 of `G1_result.md`). I built both blocks from the realization lemma with my own frame completion, diagonalized them, and summed the second-order bond over physical eigenvectors and denominators. H⁽¹⁾ = 0 to 10⁻³³ and ‖H⁽²⁾ − (−7177/8192 + 27/8192 P₃)‖_F = 2.5 × 10⁻¹⁵.
- Dressed Bell-pair closure. The 2 × 2 block [[0, g√r], [g√r, Δ]] is exact. Exact 2ⁿ degeneracy (splitting ≤ 7 × 10⁻¹⁵) and the analytic ground energy hold for six (q, n, r) with r up to 8, three couplings, and randomly rotated channel bases. The joint graded cell's compressions αγ_b and αη_b at both ends have residuals ≤ 3 × 10⁻¹⁵ (q = 3, 4).
- Port-idle identity. Fη = ηF holds exactly when no gate touches the port. One port gate, undone a step later, leaves ‖η†Fη − F‖ = 0.77.
- Second-order formula with the adjacent-bond term. For an open chain of three q = 2 cells, ‖H_eff − H⁽¹⁾ − λ²H⁽²⁾‖ falls as λ³ (4.1 × 10⁻⁵, 5.1 × 10⁻⁶, 6.4 × 10⁻⁷ at λ = 0.04, 0.02, 0.01). Dropping the three-cell term leaves a λ² residual. The Section 2 machinery is right.
- C7's selector was already checked twice, so I did not redo it. C8 §12's frozen-defect theorem is a three-line orbit argument, and I checked it by hand.

### 2. Route B at second order (new)

`route_b_second_order.py` computes the complete second-order effective chain of the joint graded dressed Bell-pair cell for (q, q₊) = (3, 1), (4, 1), (4, 2), expands it in global Jordan–Wigner Majoranas, and projects the quartic part on continuum currents J^{ab} = iψ_aψ_b, with c_{j,a} ~ ψ_{R,a} + (−1)^j ψ_{L,a}.

- Translation does forbid every mass. The only quadratic terms are on-site chemical potentials with h₂(0) = h₂(π), which move Fermi points and cannot gap chiral modes. The singlet Gross–Neveu coefficient is repulsive (+8.5α⁴ to +10.9α⁴), so C8 §11's attractive sign does not recur.
- The chemical potentials are degenerate, one value for every complex mode, because the graded chain has a global U(q₊) × U(q₋) symmetry (on ports, partners and logical spinors). At q₊ = q₋ they vanish altogether, since exchanging the two grading sectors reverses them. Each mode's Fermi points also sit at δk and π − δk, since the first-order hopping is purely imaginary.
- At q₊ = q₋ = 2 nothing cuts the flow off, and the one-loop SO(8) current-current flow from the computed couplings runs away (l ≈ 1/|G₀|) even with the gadget's exact symmetry. The free point is marginally unstable along a non-singlet direction.
- At q₊ ≠ q₋ the shift removes the charged channels, but two families stay at zero momentum, the U(n) flavour currents (n = q₊q₋) and "quartets" ψ_{aR}ψ_{bR}ψ_{aL}ψ_{bL}, whose momentum 2δk + 2(π − δk) vanishes for any shift. With exact gadget symmetry the restricted flow stays bounded. Flavour noise of relative size 10⁻² makes it run away for either sign (the isotropic point is a separatrix, and without an exact Cartan U(1) the easy-plane side is destabilized by XYZ-type anisotropy). The forward coupling controlling the quartets is attractive at q = 3 (G_QQ = −0.022) and repulsive at q = 4 (+0.059).

So route B is not "untuned at this order". Compiled history-state cells cannot be expected to keep U(q₊) × U(q₋) exactly, and any breaking opens a gap at an exponentially long but finite scale. The repair is design, not a no-go: distinct logical chemical potentials, real logical hoppings that unpair (δk, π − δk), and controlled signs of the remaining forward couplings, after which the chain is C8 §6's incommensurate Luttinger liquid.

### 3. Failure modes, ranked

| Failure mode | Severity | Status | Deciding test |
|---|---|---|---|
| Frozen markers give ≥ L degenerate ground states on prime rings (C8 §12). Repairs (mobile cells with payload transport, or global vacancy doping, C8 §17.5) are unbuilt. Mobile transport is hard here because the register does not move with the marker, and C8's vacancy transport needs a frustration-free core while the port-idle cells use signed energy terms. | Fatal under the all-ring reading, minor under a subsequence reading | Open | Jonah's ruling on scope. Then a toy: vacancy-doped frozen selector with a trivial c = 1 payload, ED or DMRG on prime L |
| Route A with the SU(N)₁ target. The swap point is a separatrix for SU(N)-breaking current-current terms (O(λ³/Δ²), O(ε_H)), as Heisenberg is for XXZ, and without an exact Cartan U(1) the easy-plane side is unstable too. | Fatal for that target | Resolved by switching | Use C8 §6's incommensurate free-fermion target, c = n at logical dimension 2ⁿ, still unbounded (c ~ log₂N) |
| Route B: no or degenerate Fermi-point shifts, runaway at q₊ = q₋, flavour and quartet channels otherwise (Section 2) | Serious | Open, repairable by design | Redo `route_b_second_order.py` with the repairs, then DMRG on the logical chain |
| Computation track unwritten: illegal-sector gaps, k_c, and zero energy on every cell length. Signed terms K⊗Y could give negative energy on wrong lengths and beat the selection margin t/b⁴ ≪ Δ. Terms K⊗(1+Y)/2 with K ≥ 0 and K\|Ω⟩ = 0 are positive and still exact on the code, which fixes this if they suffice. | Serious | Open | Write it out |
| Thermodynamic Schrieffer–Wolff. Bravyi–DiVincenzo–Loss give extensive errors [read, 1105.0675]. Lie–Schwinger block diagonalization (Del Vecchio, Fröhlich, Pizzo, Rossi) needs unique block ground states [recalled]. Nothing rigorous covers degenerate blocks with a critical effective model. | Fatal for a theorem, minor for physics | Open | None within reach, so accept SW plus RG |
| Route B dimension count needs k_c independent of q (finite gate menu on q-dimensional registers via Solovay–Kitaev). | Minor | Plausible | Part of writing the track |
| b-fold spontaneous translation breaking on L ∈ bℤ. The spectrum is b exactly degenerate copies of one CFT, which is a direct sum, so c is not multiplied by b. | Minor | Resolved | None needed |

### 4. Points examined and cleared

- Hierarchy. For route A, with Δ = 1, the ordering J = ρ²λ² ≪ λ ≪ ρ², ε_H ≪ ρ², J_c ≥ 2T²(2M + D_req), walker gap 3π²t/b² ≫ λ‖v‖, selection margin t/b⁴ ≫ λ²‖v‖² (the scalar second-order energy can differ between cell lengths), K > 4t is consistent. Choose λ, then t, then J_c. Every scale is fixed by N and L → ∞ is taken last, as the conjecture requires. Local SW needs only local gaps, so the extensive logical bandwidth does not enter.
- Resonances. A cell excitation hops through the ports at first order, so excited levels form bands of width O(λ) around Δ. They stay above the logical band. No low-energy resonance appears in the frozen route.
- Junk sectors. On L ∈ bℤ the frozen route has no gapless sector besides the payload. Illegal selector sectors cost ≥ K − 4t, clock towers ≥ J_c/2T², defect pairs 2γ_b with a flat (immobile) band above that gap. The mobile route adds the c = 1 length fluid, which is harmless as an extra CFT. Its optimal density ν* = 1/2 − 1/[π²(b + 1/2)] (C8 §13) sits close to the Peierls-resonant ν = 1/2, which gives a marginal coupling over about b cells and needs checking once a payload rides on it.
- Route A IR with the incommensurate target. Masses and inter-flavour backscattering carry momenta 2k_a, k_a ± k_b ≠ 0. The quartets of Section 2 carry 2θ_a + 2θ_b ≠ 0 because the synthesized hoppings have phases, and C8's condition (29) excludes all integer resonances. Twist fields have dimension n/4, irrelevant for n ≥ 9. At the RG level this target is stable without symmetry.
- Consistency with general bounds. Thermal entropy density πcT/3v ≤ log d inside the conformal window, with v = b·v_cell and window top ~ v_cell, requires b ≳ πc/(3 log d). The construction has b ≥ T_clock ≫ N ≥ c, so it passes. The area law for degenerate ground spaces [recalled] and C8's clustering tradeoff bound code size and end-to-end compressions only at fixed cell gap, and the cell gaps here shrink with N. Gosset–Mozgunov [recalled] and the claim that a CFT cannot emerge from a frustration-free chain [read, cited in 2610.24154] do not apply, since the walker selector and the payload are frustrated. A search found no theorem bounding c or edge degeneracy by a function of d for frustrated critical chains [read, search results only], and no counterexample in the literature [read, 2403.17242, 2506.16424].

### 5. Verdict

No fatal flaw in the finite-dimensional mathematics. The construction is a coherent physics-level route to unbounded c, and every specific obstruction tested here is either a design obligation or a matter of how the conjecture is read. The weakest links, in order, are the ring-length requirement, the unwritten computation track, and infrared control of whichever target is chosen. Route A with the incommensurate free-fermion target is the version to pursue, since it violates c–d by itself once n > 8k_c and avoids both the separatrix of SU(N)₁ and the degeneracies of route B.

## Reproducing

- `decisive_identities.py` (numpy, seconds): G1 at (d, N) = (2, 3) from physical eigenvectors, Bell-pair closure, port-idle identity. Output `decisive_identities_stdout.txt`.
- `route_b_second_order.py` (numpy, scipy, about 4 minutes, `OPENBLAS_NUM_THREADS=1`). Three-cell ED validation, then the second-order chain, its quadratic and quartic parts, and one-loop flows. Output `route_b_stdout.txt`. `--big` adds (q, q₊) = (4, 2) at cell dimension 4096, about 10 minutes, output `route_b_big_stdout.txt`. `--ed` adds a slow two-cell q = 3 check.

## Relation to the rest of the repository

This reviews REPORT §9.1 to §9.8 and does not change "Where things stand". Three edits seem warranted. The programmable-cell README should drop "untuned at this order" for route B, the SU(N)₁ target should give way to the incommensurate one, and the ring-length ruling (`contrib/unbounded_c/` §4(v), C8 §12) should come first, since under the all-ring reading the frozen-selector chains are outside the conjecture's hypothesis.
