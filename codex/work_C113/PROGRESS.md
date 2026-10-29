# C113 progress

2026-10-29

[Proved, work record] Read the project memory, the complete bound background, REPORT sections 9.68 to 9.77, the modular-route README, C105's memory, result and relevant proof sections, and the smooth-lemma specification. All edits are confined to work_C113.

[Assumed, task] Determine the precise continuum reduction, its converse, and its scope for gapped and nonconformal chains. The proposed universal lattice inequality is not assumed true.

[Proved, preliminary obstruction] The uniform Hamiltonian and gauge-optimized M are invariant under a bond divergence. The weighted Hamiltonian is not. A lemma quantified over all representatives cannot have the proposed right side. The effect survives smooth modulation and requires an explicit gauge prescription.

[Heuristic, being checked] Infrared irrelevant operators have suppressed long-distance contributions, but integrating out short-distance modes can also produce an analytic two-gradient vacuum-energy term. The requested assertion that every correction is higher order cannot be inferred from the uniform CFT spectrum alone.

[Assumed, next steps] Derive the exact stress-tensor formula including the circle term, give explicit lattice counterexamples where needed, distinguish a sufficient spectral-response inequality from an equivalent continuum statement, and analyze the XXZ ferromagnetic and free quadratic-band edges.

[Proved] Completed the chiral stress calculation, kappa = 1/(24 pi), with the finite-circle harmonic-length term retained. A thermodynamic limit over repeated modulation periods isolates the anomaly. The exact continuum formula covers neutral sectors and chiral asymmetry without on-site symmetry assumptions.

[Proved] Found an exact independent-spin gapped example with Def/G -> lambda^2/(2g) at fixed M=g and gap=2g. Stacking it with XX preserves A=1 while adding an arbitrary leading response. This disproves unqualified leading-coefficient matching and the converse.

[Proved] Found a fixed-bond ferromagnetic z=2 counterexample to the literal lemma. A coherent product-state witness gives a positive limiting deficit while G vanishes as 1/R. The mechanism is a telescoping longitudinal field. The canonical divergence-free ferromagnet and the free quadratic-band vacuum edge instead have zero deficit.

[Proved, conditional] Derived a gap and current-correlation bound for a gapped interpolation, and a weaker sufficient reduction using the positive inverse-cube current spectral moment. Explained why the latter implies the A bound under spectral matching without assuming that all microscopic response equals the anomaly.

[Numerical] Created and ran checks.py with python3. All small algebra and quadrature assertions passed. The first invocation with python failed because that executable is absent. No heavy local computation or Typhon job was needed.

[Proved, work record] Wrote RESULT.md, PROOF.md and local STATE.md. The parent project STATE.md and all shared files are unchanged under the task-specific restriction.

[Proved, audit strengthening] The obstruction survives a unique minimum-spread gauge. Constructed a d=4 direct-sum model with active AA ferromagnetic block, BB block 16 S.S+8I and mixed blocks 12I. BB extremal states enforce M=8; saturation proves uniqueness against every on-site Hermitian gauge operator, including off-diagonal label blocks. Domain-energy bounds establish the active uniform infrared. The weighted active coherent-state witness still has a positive limiting deficit.

[Proved] Added a unique-gapped-ground-state variant at the same d=4,M=8, with gap2g and Def/G>=1/(2g), and a critical XXZ variant with susceptibility divergence under the standard infrared dictionary. Minimum-spread selection alone is therefore disproved as a repair.

[Numerical] Extended the small checks to the 16 by 16 embedded bond and the 8 by 8 three-spin BB positivity input. Both passed. No new large computation was required.

[Proved, final audit] Checked the four-state construction's gauge uniqueness against off-diagonal as well as label-preserving gauges, the domain-energy bounds, and the gapped variant's exact gap. The work leaves no claim of unconditional equivalence. A different infrared-adapted representative with norm controlled by M remains a possible sufficient target, whereas minimum-spread selection is explicitly ruled out. Final files are RESULT.md, PROOF.md, PROGRESS.md, STATE.md, checks.py and checks.json.
