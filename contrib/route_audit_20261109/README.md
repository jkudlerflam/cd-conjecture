# Independent route audit and improved contact bound

Date: 2026-11-09. AI assistance: OpenAI Codex (GPT-6 family). No human
review, new CFT measurement, or unconditional counterexample is claimed.
The contact proof was independently reviewed by the coordinating Codex
agent (AI only).
The audit started at `a45ecf5` and incorporated the new C175 result in
`1a4015f` while work was in progress.

## Result

**The d=13 mobile construction remains the strongest mature route to a
counterexample.** It has a literal nearest-neighbour, one-site invariant
bond and certified finite material ingredients, with predicted total
Virasoro central charge c=14. Its remaining microscopic matching and
relative defect estimates are substantive hypotheses, even when the
working standard is theoretical physics rather than a complete
constructive-RG theorem.

The new contribution here is a sharper theorem for one of those inputs:
under precisely C92's legal-space, positive-length-tilt assumptions,

\[
 \langle A_j\rangle>1/264,
 \qquad
 F_{\rm pinned}-F_{\rm unpinned}
 \ge k[J/135-T\log24576].
\]

C92's corresponding conservative bounds were 1/4096 and J/4096.
The improvement comes from retaining the binary entropy of the
vacancy-removal channel's environment, rather than bounding it by
log 2. See [the proof and full scope](CONTACT_BOUND.md) and the
[exact arithmetic checks](contact_checks.json). This strengthens the
projected pin estimate; it does not close the whole construction.

## Route ranking and decisive next step

| Route | Evidence currently available | Decisive unresolved step | Priority |
|---|---|---|---|
| d=13 mobile, predicted c=14 | Explicit bond, certified frozen material member, independent local grammar and band checks, legal projected pin theorem | A dressed gap/contact representation at nonzero J and a common-bulk-subtracted wrong-length/finite-K comparison, followed by the multichannel infrared identification | Main analytic route |
| d=11 frozen, predicted c=11 | Smaller explicit bond, certified 22 Fermi points, frozen sector comparisons on compatible rings | Many-cell infrared transfer; qualification to lengths divisible by 13 and the translation-breaking sectors must remain explicit | Simpler analytic control, less clean as a full-ring counterexample |
| d=2 Fibonacci cusp, proposed c=2 | Reproducible effective Hamiltonian and suggestive unconverged entropy fits | C175 now finds a bound pair below the atomic continuum; finite-density pairing gap and full-qubit scaling decide whether c=2 survives | Bounded finite-density diagnostic, no longer the leading dilute mechanism |
| d=4464, unbounded proposed c | Finite construction and low-n projected numerical controls | The stated microscopic-to-infrared input for every chosen n | Long-term extension; does not shorten the current missing matching step |

The original conjecture explicitly allows a collection of unitary CFT
sectors with different velocities. Velocity inequality alone is not an
exclusion. [Latorre and Sierra, section 2](https://arxiv.org/html/2403.17242v2#S0.SS2).

For the d=2 route, [C175](../../codex/work_C175/RESULT.md) distinguishes
its negative dilute result from the still unresolved density near 0.17.
Its binding and cluster data prevent using the atom-only two-valley XXZ
mapping as a ground-state argument there. A bound pair does not by
itself rule out an intermediate-density two-component phase.

## The d=13 scope ledger

| Ingredient | What is actually established | What cannot be inferred |
|---|---|---|
| Literal bond and legal coding | Independent reconstruction of the d=13 bond, charge/pin conservation, unique vacancy colouring and finite-order grammar checks | A thermodynamic CFT or convergence of all perturbative orders |
| Frozen band certificate | 26 simple Fermi points and low-order nonresonance for the certified n=13 material member | The same bands at arbitrary fixed J/t after adding vacancies |
| [C91](../../codex/work_C91/RESULT.md) ultraviolet filter | A volume-independent first-moment bound for the legal compressed Hamiltonian | Replacing filtered contact by its mean, or a dressed massive auxiliary gap |
| [C93](../../codex/work_C93/RESULT.md) | Block parity representation; a separately specified strong-spatial-norm ultraviolet schedule; leading colour-field budget | A bound on the retained dynamical contact-clock coupling or temporal influence kernel |
| [C92](../../codex/work_C92/RESULT.md), strengthened here | Uniform contact and activated pin estimates at infinite grammar penalty in a positive-length-tilt ensemble | A canonical fixed-vacancy contact bound, finite-K pin theorem or wrong-length comparison |
| [C94](../../codex/work_C94/RESULT.md) | Counterexample to a proposed temperature-uniform finite-volume response bound; conditional quadratic comparison | A refutation of the thermodynamic construction, which permits a global subextensive term |
| [C95](../../codex/work_C95/RESULT.md), [C96](../../codex/work_C96/RESULT.md) | Exact free-charge response and tagged-kernel formula; weaker spatial criterion; quadratic resummation controls | The full relative low-temperature interacting remainder, uniformly in wrong lengths and nonzero port |
| Multichannel Luttinger argument | Conditional operator and resonance analysis | Verification that the microscopic mobile chain satisfies its dynamical matching assumptions |

This ledger distinguishes a missing proof of a standard infrared
principle from a missing verification of the model-specific hypotheses
needed to apply that principle. Both are present in the source record.

## Parameter prescription: a concrete correction to carry forward

The [red-team report](../redteam_d13/README.md) found two fixable failures
of the wording "fixed positive J,t" in the original displayed member:

1. At the fixed marker chemical potential, vacancy response favours
   length 14 over 15 above J/t about 0.0056 at second order.
2. The mobile rail traversal agrees with the frozen contact-average
   dictionary only when D J/t is sufficiently small. The toy traversal
   suggests J/t below roughly 2^-175 for the raw certificate tolerance.
   This is a heuristic extrapolation of a checked toy calculation, not
   a uniform proof for the entire chain.

The leading colour field simultaneously requires

\[
 K\ge\max(2^{2560},2^{2243}(t/J)^2).
\]

There is no algebraic conflict among these finite requirements. C91
already supplies J/t=2^-5000 and K=2^14000 as a restricted example for
its first-moment ultraviolet and leading colour budgets. That example
satisfies the displayed C93 colour condition with 1757 binary orders
of margin (14000-(2243+10000)), and is far below both red-team
thresholds. It is **not** a certificate for C93's much stronger
spatial-weighted schedule or for the missing infrared hypotheses.

The published `codex/work_C88/MEMBER.json` and
`contrib/small_d_note/note.tex` still use an implicit finite-K bound with
fixed positive J,t. Any future finite-member presentation should state
the fast-vacancy condition explicitly and then distinguish the proved
finite budgets from the still conditional mobile dictionary. Making
J/t small does not on its own remove bare clock leakage: C93 equation
(10) gives a leakage matrix element 0.15071494 J, while its reference
clock gap is itself proportional to J. Their ratio stays finite as
J/t tends to zero. A dressed low-energy bath estimate is needed.

## Reproduction and review

`python3 -B contrib/route_audit_20261109/check_contact_bound.py` uses only
the standard library and regenerates `contact_checks.json`. All
certificate assertions use integer or rational arithmetic. Printed
floating-point roots are diagnostic only. The proof's entropy argument
and inherited operator identities require mathematical review in
addition to this arithmetic check.

`check_entropy_channel.py` (NumPy) checks 72 finite complex density
matrices and seven examples saturating the entropy bound; all pass,
with subadditivity slack no smaller than -2.3e-16 from roundoff.
The script checks the entropy-channel step only, without making a
Hamiltonian truncation or a thermodynamic claim.

The existing `contrib/small_d_note/regen_bonds.py` was also rerun
read-only. It independently reconstructs the displayed d=11 and d=13
move tables and compares them to the repository bond implementations
at diagnostic parameter points. Both maximum entry differences are
1.43e-14 or smaller; the pin identity residual is zero. The output is retained in
`bond_rebuild.log`. This is corroboration of the finite bond, not a
new microscopic central-charge measurement.

Only this contribution directory was modified by this audit. The main
coordinator handles repository synchronization and any broader status
changes after review.
