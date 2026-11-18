# What the earlier obstructions establish

The following repository reports were read, including their proofs and
implementation conventions. The new seed does not claim to refute their
stated theorems.

| Source | Established scope | What must not be inferred |
|---|---|---|
| [C84](../../codex/work_C84/RESULT.md) | An eleven-state frozen alphabet, an n=11 finite band certificate, and a zero-port selector with explicit finite-parameter budgets | An all-n CFT theorem, nonzero-port ground-sector selection, or a universal lower bound on physical d |
| [C86](../../codex/work_C86/PROOF.md), sections 3--4 | Signed-Gram hopping and reciprocal spectra in the common-carrier/reflection class; failure of the positive-metric and pair-sum tests | That indefinite Hermitian hopping or a resonant pair transfer necessarily gaps the spectrum |
| C86, section 7 | Eleven-letter lower bound **within** the product-token, matched-species, point-port class **and the retained nonresonant certificate** | A universal d>=11 result for interacting liquids or for designs with additional projected physical terms |
| [C87](../../codex/work_C87/PROOF.md) | Restricted token-phase, storage, port-dressing and marker-charge obstructions | That every data-dependent interaction, contextual encoding or resonant low-energy model is excluded |
| C87, section 6 | Abstract complementary channels can cancel an interacting second-order coefficient | An actual <=10-state nearest-neighbour realization of those four channels |

For a real common carrier with point sources, C86's coefficient is
`h(k)=epsilon I+T exp(ik)+T^T exp(-ik)` with T real and no generated
onsite M at leading conversion order. It is reciprocal. With four
simple roots the opposite-momentum pairs necessarily collide. That
collision is a **kinematic permission** for scattering. Its amplitude,
sign, forward couplings and RG flow require a separate calculation.

The unmodified nine-state single rail is not automatically a multichannel
metal either. At mass 2 the two-flavour symbol's two band ranges on the
saved grid are approximately `[-.0135421,-.00180516]` and their positive
reflection. The three-flavour symbol has two separated outer bands and
a narrow central band; its ranges are saved in `audit.json`. Grid
inspection sees only one Fermi pair at a time under a common field.
This is a finite diagnostic, not an all-mass or all-n theorem.

The new ingredient is consequently explicit: weak physical fields and
background-bit swaps generate a non-scalar onsite matrix in logical
space. They add no states, but they leave the **common-field-only**
Bloch ansatz above. The real seed deliberately retains reciprocity,
and its resonant quartic coefficient is computed rather than removed.

The labels *frozen*, *finite n*, and *conditional infrared interpretation*
are retained throughout. This audit does not identify the collaborator's
reported d=10 unbounded-c source. The source-verification task is separate.
