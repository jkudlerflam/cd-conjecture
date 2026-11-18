# Smaller-d priority and the reported ten-state family

Date: 2026-11-17. OpenAI Codex, GPT-6 family.
Status: source audit and research assignments, not a new construction,
independent certificate verification, or impossibility theorem.

The user reports a collaborator's possible d=10 unbounded-c example and
asks that it take priority over d=20/22 if confirmed. The user also
explicitly prioritizes the smallest possible physical on-site dimension.
This preference is now recorded in the root AGENTS.md alongside the
existing preference for fewer terms, moderate coefficient ratios and
less delicate tuning.

## What is and is not currently identified

**The reported d=10 arbitrary-c example has not yet been located.** A
source was requested from the user without suspending independent work.
Searches covered the current checkout, its fetched alternate collaboration
branch, the available project-chat inventory, and targeted public searches.
Failure to locate it is not evidence that no such construction exists.
No external collaborator was contacted.

The following existing contributions must not be conflated:

| Physical d | Construction and recorded evidence | Scope limitation |
| --- | --- | --- |
| 2 | [Four-term V=4 candidate](../v4_bulk_diagnosis_20261117/README.md), provisional finite-size entropy slope near 2.57 | Its thermodynamic c is unknown; largest accuracy checks are unfinished |
| 10 | [Nearest-neighbour demonstrator](../routeA_demonstrator/README.md), one logical flavour, numerical c approximately one | A mechanism demonstration, not an arbitrary-c family |
| 10 | [Mobile two-flavour toy](../mobile_cells/README.md), proposed total c=3 | Even that c is not established; no arbitrary-flavour construction at fixed d=10 |
| 10 | [C53 exact completion](../../codex/work_C53/RESULT.md), exact logical XXZ in the stated frozen sector and a free c=1 member | Uses four-site physical interactions and does not give unbounded c |
| 11 | [C84 frozen bond](../../codex/work_C84/RESULT.md), certified finite n=11 band member, predicted c=11 | The certificate does not establish arbitrary-n margins or the microscopic interacting CFT |
| 12 | [Single-vacancy mobile candidate](../mobile_below_13/README.md), proposed c=12 using C84's n=11 member | Additional slow-vacancy and interaction arguments are heuristic; no measured c |
| 13 | [C88 two-colour mobile bond](../../codex/work_C88/RESULT.md), certified n=13 material member, predicted total c=14 | Mobile sector selection, matching and infrared hypotheses remain |

The eleven- and thirteen-state constructions are also written out in
[the small-d note](../small_d_note/note.tex), which explicitly calls
their CFT conclusions conditional. The twelve-state mobile proposal is
smaller than thirteen but carries additional dynamical assumptions: its
vacancy motion depends on the clock and is very slow, with an attractive
induced interaction requiring a compensating contact term. Small d and
a simple or well-controlled Hamiltonian are separate comparison axes.

The recent d=20/22 work was aimed at an all-n construction and stronger
clock/defect control. Those dimensions must not be presented as the
smallest existing finite-member candidates. Conversely, replacing them
by a finite n=11 or n=13 certificate does not establish unbounded c.
The useful next step is to test whether the all-n mechanism and relevant
proofs can be transferred to the smaller alphabets.

## Known restrictions do not rule out d=10

[C86](../../codex/work_C86/RESULT.md) and
[C87](../../codex/work_C87/RESULT.md) found no certified design with at
most ten letters. Their lower bounds assume particular encodings, token
motion and port rules; they are not lower bounds on d for all physical
constructions. Several failures concern the existing sufficient band
certificate, notably reciprocal spectra and coincident Fermi-point pair
sums. These do not alone prove that the interacting system is massive.

One concrete remaining possibility is a smaller resonant construction
whose allowed interactions are irrelevant or symmetry-forbidden. Testing
it requires the actual generated interactions of a literal bond, not just
counting free-band crossings. Another possibility is a new contextual
encoding or port mechanism outside the restricted counting assumptions.

## Fourth-round assignments, delivered and active

1. **Establish the qubit counterexample:** analyze the V=4 saved states
   and accuracy stages to distinguish numerical error, density/parity
   drift and boundary corrections from an extra gapless sector. Test
   whether the observed correlation momenta support a one- or two-mode
   description under the actual microscopic symmetries. Keep the previous
   battery-saving production stop in force; use saved-data analysis and
   short bounded checks, and specify the cheapest decisive continuation.
   Output: `contrib/v4_interpretation_audit_20261117/`.

2. **Bounds/mechanism chat:** audit the resonant <=10-state loophole in
   C86/C87. Start with the smallest explicit retained model at two or
   three logical flavours, its full Bloch symbol and the leading actual
   interaction channels. Determine whether resonance produces an
   instability or merely falls outside the sufficient non-resonance
   certificate. A narrow no-go result, an explicit viable seed or a
   precise obstruction are all useful; no generic d>=11 bound is assumed.
   Output: `contrib/ten_state_resonance_audit_20261117/`.

3. **Stabilize the fixed-d unbounded-c family:** verify the reported
   d=10 source if it becomes available and make it primary if credible.
   Meanwhile audit the smaller existing bonds and test transfer of the
   all-n band mechanism to C84's eleven-state shifted-bit alphabet.
   Trace the previously reported large-n certificate deterioration and
   distinguish a bad coefficient prescription from an architectural
   obstruction. Include the inactive endpoint and finite-K splitting;
   changing clock bias must not conceal new physical letters. Retain the
   twelve-state heuristic mobile route in the inventory. No further
   d=20/22-specific extension is assigned unless needed for a smaller bond.
   Output: `contrib/minimal_d_unbounded_audit_20261117/`.

All three existing chats were observed in progress after receiving these
prompts. They retain their settings and use separate contribution folders.
No new chats or external messages were created. No new long numerical
campaign was started by this coordinating chat.

## Review boundary

This update reads and reconciles source statements and does not rerun
their finite certificates or claim independent review of their proofs.
Documentation links and the patch were checked. The top-level overview
is clarified where it previously promoted conditional predictions to
established counterexamples. The current plan explicitly distinguishes
microscopic d, finite-member evidence, all-n control, and CFT identification.
