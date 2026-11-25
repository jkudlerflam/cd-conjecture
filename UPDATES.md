# Dated updates

Announcements of each research round, newest first. They were moved here from the README and STATUS.md on 2026-11-24 without changes. Several are superseded by later entries, so read the summary in the README's "Where things stand" first. The older material at the bottom is not in date order.

## 2026-11-25: quartet estimator calibration and IAS pilot

[Numerical.] The [quartet study](contrib/quartet_typhon_20261125/README.md)
reproduces the six-file table exactly but finds a large filtering bias:
the published 48-site procedure fits an exact asymptotic p=4 two-boson
control as p=5.381 over r=32–56. Forward-filtered fitting repairs that
calibration; operator/envelope dependence still prevents a replacement
charge-four dimension. Thus the earlier x4≈2.05±0.2 is not a calibrated
estimate. Dense Hamiltonian, field, gauge and extended-correlator checks
pass. IAS pilot **17764244** completed in 50:58 and saved an unconverged
L128 chi256 checkpoint; chi384/512 were not reached. The h=1.772 control
matched same-chi entropy references through chi256. Target quartet signals
still change by order one between chi128 and 256, and centred enlarged
matrices do not certify the neutral-rank assumption. The conditional
L256/384 expansion was not launched because these checks failed. Recipes,
small data, source hashes and convergence limits are published (REPORT §9.219).

## 2026-11-25: proposed IAS numerical campaign

The [cluster plan](contrib/cluster_numerics_plan_20261125/README.md) prioritizes
converged quartet measurements and operator separation on the simple d=2
chain, then bulk modes and velocities. A smaller original-d=10 contact
response campaign targets the missing mobile-sector bounds. It incorporates
the new pi-pair and direct-quartet data of §§9.216–9.217, includes CPU
pilots and resource gates, and records the live Typhon capacity check.
This is a proposal; no research jobs were submitted and no CFT verdict changed.

## 2026-11-24: the quartet correlator on the four-term chain

The k = 0 quartet correlator, measured directly on the stored open chains at δ = 3/16, h = 1.783 (L = 256 and 384), gives x4 ≈ 2.05 ± 0.2. That is one error bar above the vertex inference 1.84 and within one of marginal, with the L drift pointing toward irrelevant. The signal is 7 to 23 times weaker in the c = 1 control. At the phase-1.53 bond the quartet stays momentum-mismatched (period about 430 sites). The independent code gives c(128 → 256) = 1.83 in even parity, against TeNPy's 1.76, and 2.03 in odd parity (REPORT §9.217).

## 2026-11-24: three checks on the four-term chain

An independent DMRG code reproduces the c ≈ 2 open-chain states at δ = 3/16, h = 1.783 to L = 128. c(h) between L = 128 and 256 steps from 1 to about 2 between h = 1.778 and 1.780 and stays at 1.7 to 2.7 across the wedge. The pair correlator in the wedge has a q = π component that the c = 1 neighbour lacks, so the quartet is commensurate, with dimension about 1.84 if the pair is a vertex operator. That is the one result so far that argues against an asymptotic c = 2. Koo–Saleur on rings (ED to L = 24, puMPS to L = 64) passes single-chain calibrations but cannot resolve this chain or a known c = 2 ladder at reachable bond dimension (REPORT §9.216).

## Moved from STATUS.md on 2026-11-24

These sat at the top of STATUS.md and inside its constructions section. Their content is folded into the current STATUS.md.

**First CFT audit completed (2026-11-23):** The managed investigation of
the fixed delta=3/16,h=1.783 qubit Hamiltonian is published, including an
independent AI audit. The finite-scale entropy evidence remains near two;
the asymptotic CFT remains unidentified. Under a matched generic two-boson
operator dictionary, the quoted exponents make the quartet the only
relevant generic locking harmonic, even with mixing. That dictionary and
the target coupling are unmeasured; a different assignment reproduces the
same two exponents without relevant generic locking. See the
[audited findings](contrib/qubit_ir_cft_team_20261123/README.md) and
[prioritized missing measurements](contrib/qubit_ir_cft_team_20261123/TODO.md):
pair momentum and primary identity, independent bulk modes and velocities,
normalized response, compactification, and longer-scale stability.

**New field-diagram investigation (2026-11-18):** One new GPT-6 Astra
Max chat is active on the h-dependent phase diagram of
H=sum[(9/16)XY+(7/16)YX+ZZ]-h sum Z. The
[assignment](contrib/research_coordination_20261116/CHIRALITY_FIELD_DIAGRAM_20261118.md)
covers the full field range and emphasizes the proposed c=3/2 transition
and unresolved nearby c about 2 readings. This is a research deployment,
not a new phase result; long production remains stopped.

**Stability round completed (2026-11-18):** All four GPT-6 Astra Max
chats completed. The modified d=9 bond gains volume-independent local
bounds through fourth port order and exact nonzero-port ground-sector
selection on compatible rings. Qubit work computes a connected
finite-density quartet response, exposes a bare-versus-dressed momentum
discrepancy, and specifies the missing full response measurement. The
original d=10 mobile model gains a uniform theorem for the canonical
zero-momentum contact response at half occupancy. All-scale CFT matching
and the full mobile selector remain open; no new microscopic central
charge is established. See the
[results](contrib/research_coordination_20261116/STABILITY_RESULTS_20261118.md)
and [assignments](contrib/research_coordination_20261116/STABILITY_ROUND_20261118.md).

**Current research priority and completed round (2026-11-18):** Large c
at microscopic d=2 and unbounded c at fixed small d. All four GPT-6 Astra
Max chats completed their assignments. The strongest new result is an
arbitrary-n leading-band theorem for a modified markerless d=9 bond.
The d=10 mobile work certifies selected finite sector comparisons and an
exact charge-field improvement. Qubit work computes finite quartet mixing
and proves crystal preemption in the tested three-valley leading slice.
No new microscopic central charge is established. See the
[completed results](contrib/research_coordination_20261116/FOCUSED_RESULTS_20261118.md)
and [assignment record](contrib/research_coordination_20261116/FOCUSED_ROUND_20261118.md).
The earlier battery-saving stop on long qubit production remains in force.

**Priority set on 2026-11-17:** Minimize physical on-site dimension while
also tracking term counts, coefficient hierarchy and strength of evidence.
A collaborator's reported d=10 arbitrary-c example had not yet been
identified then; the publication below resolves that source search.
The known d=10 demonstrator and two-flavour toy did not prove that claim.
That round audited <=10-state mechanisms and transfer
to the existing eleven-state alphabet, before further d=20/22-specific
work. A finite n=11 or n=13 band certificate must not be promoted to an
all-n result. See the [source ledger and assignments](contrib/research_coordination_20261116/SMALLEST_D_FOLLOWUP_20261117.md).

**Earlier follow-up (2026-11-18, before publication below):** Live GitHub
issues, comments, both branches and targeted history searches had not
located the reported d=10 arbitrary-c model. Results gave a separate
single-rail nine-state two-flavour resonant seed and an
eleven-state all-n clock/zero-port-selector/simple-tree reduction; neither
establishes arbitrary-c physics. The next calculations test three flavours
in the nine-state seed and all-n loop splitting in the eleven-state model.
The qubit chat studies symmetry-allowed locking operators while production
remains stopped. [Audit and scope](contrib/research_coordination_20261116/GITHUB_D10_AUDIT_20261118.md).

**Located (2026-11-18):** The reported d=10 candidate family is
[contrib/markerless_d9_d10/](contrib/markerless_d9_d10/README.md), written
on 2026-11-12 and kept in a local checkout until today. Nine material
states without markers have predicted c=n with frozen cells; one
uncoloured vacancy gives predicted c=n+1 at d=10. Its status is that of
the d=11 and d=13 bonds: certified logical band data at the listed n,
conditional infrared step, and no all-n theorem for that original
parameter schedule. The newer modified all-n band proof and finite mobile
comparisons are distinguished above. Thermodynamic mobile selection and
the dressed-band J/t budget remain open. See REPORT for the referee outcome.

**That earlier round is now complete:** The tested single-rail nine-state three-flavour point
has the desired bands but a coupled leading RG instability. The eleven-state
loop had analytically separated pair sums for a new n<=4 family; its all-n
moment condition was left open then and is resolved by the new sign proof
above. The qubit analysis allows quartet and mixed
locking and corrects the drift covariance used in C180; C191's chiral
implementation agrees with the correction. No c value was established and
no DMRG was restarted. [Results and qualifications](contrib/research_coordination_20261116/ROUND5_RESULTS_20261118.md).

**Finite-density follow-up (2026-11-18):** The
[connected quartet calculation](contrib/qubit_detuned_quartet_20261118/README.md)
gives nonzero finite-ring response coefficients at the fixed asymmetric
point, but they are not yet a continuum coupling. The
[dressed-mode audit](contrib/qubit_finite_density_modes_20261118/README.md)
finds a finite-window carrier shift of opposite sign to the bare-band
estimate and proves that the usual leading parameter probes have
insufficient response rank. A normalized four-channel response, target
velocities and a controlled thermodynamic mismatch are still needed to
apply the new coupled locking-versus-detuning criterion. These results
do not validate the order-of-magnitude stability estimate in item 5.

## Moved from the README on 2026-11-24

**Qubit update (2026-11-23):** The earlier c=3/2 interpretation of the
simple chirality chain was withdrawn on November 19. Its proposed c
approximately two region now has both infinite-chain and open-chain
entropy evidence, with the cleanest new point at delta=3/16,h=1.783.
The latest correlator estimates support the two-valley interpretation
there, but do not settle long-distance locking or establish a microscopic
CFT. The older large-blockade and phase-1.53 candidates remain conditional.
The [current review](contrib/qubit_c2_status_20261123/README.md) compares
the evidence, convergence limits and coefficient simplicity, and
reproduces the central-window fits from the committed tables.

The subsequent [independently audited team investigation](contrib/qubit_ir_cft_team_20261123/README.md)
identifies missing operator assignments and conditional locking tests;
it does not identify the asymptotic CFT. The two reported exponents can
admit different stability conclusions under different dictionaries.
The [prioritized to-do list](contrib/qubit_ir_cft_team_20261123/TODO.md)
specifies the missing pair momenta and quadratures, bulk modes and
velocities, normalized response, compactification and locking amplitudes,
with minimum measurements and data/resource requirements. See REPORT
§9.215 for corrections to the earlier interpretation.

**Stability round completed (2026-11-18):** The modified d=9 family now
has finite-order local interaction bounds independent of chain length
and exact nonzero-port ground-sector selection on compatible rings.
Qubit work computes a connected finite-density quartet response and
finds that saved carrier shifts disagree in sign with bare-band estimates;
the continuum matching needed to settle c=2 remains missing. The original
d=10 mobile model gains a uniform bound on one half-occupancy response
channel, with full thermodynamic selection still open. All four GPT-6 Astra
Max chats completed. No new microscopic c is established and the
coefficient hierarchies are not improved. The
[results and remaining gaps](contrib/research_coordination_20261116/STABILITY_RESULTS_20261118.md)
also distinguish the collaborator's separate, simpler c=3/2 numerical
qubit candidate. Long production remains stopped.

**Focused round completed (2026-11-18):** A new moment-sign proof gives
an arbitrary-n leading quadratic band family at physical d=9, using a
modified biased clock and L1 readout. The original d=10 mobile model gains
integer-certified finite sector comparisons and an exact conserved-charge
field improvement. Qubit work computes finite interacting quartet mixing
and proves crystal preemption in the tested number-conserving three-valley
leading model. No new microscopic c or unbounded-c CFT is established, and
the small-d coefficient hierarchy remains severe. All four GPT-6 Astra Max
chats have completed; [results and scope](contrib/research_coordination_20261116/FOCUSED_RESULTS_20261118.md).

**Ten-state family located (2026-11-18):** The reported d=10 candidate
family is [contrib/markerless_d9_d10/](contrib/markerless_d9_d10/README.md),
dated 2026-11-12 and published today. It removes the marker letters of the
d=11 and d=13 bonds: the predictions are c=n on nine-state frozen cells
and c=n+1 after adding one uncoloured vacancy, at d=10. The legal sector
is exactly invariant, eliminating finite-K leakage. Logical band data
are certified at n=9,10,11,12,14,16,20, not for all n. The infrared step
remains conditional, and the corrected mobile bond still needs a
quantitative J/t budget and multi-cell sector checks. This is now the
priority fixed-d family. An AI referee report is in the folder.

**Earlier completed round (2026-11-18):** The separate single-rail nine-state three-flavour
bond has six certified simple Fermi roots, but its generated interactions
destabilize the tested three-mode Gaussian state in the coupled leading
RG flow. The eleven-state loop gives full pair-sum separation for a new
analytic n=1,...,4 family and an explicit remaining all-n nonvanishing
criterion, without adding states. The qubit audit permits quartet and
mixed locking, corrects C180's drift covariance, and identifies the
momentum-pi pair exponent as a useful future diagnostic. No microscopic
central charge or unbounded-c theorem follows. See the
[completed-round summary](contrib/research_coordination_20261116/ROUND5_RESULTS_20261118.md).

**Earlier source audit (2026-11-18, before publication above):** Direct
GitHub branch, issue, comment and history checks found no d=10 arbitrary-c
source then; that source search is now superseded.
New work supplies a literal nine-state two-flavour resonant seed with a
favorable leading RG flow, and transfers clock, zero-port selector and
simple-band arguments to eleven states for arbitrary finite n. The first
does not establish arbitrary flavours; the second retains exact Cooper
resonances in its tree and an open all-n loop-splitting problem. Neither
establishes a microscopic CFT or unbounded c. The new tasks address those
gaps, plus qubit mode-locking criteria without restarting DMRG. See the
[source check and follow-ups](contrib/research_coordination_20261116/GITHUB_D10_AUDIT_20261118.md).

**Smaller-d priority and source audit (2026-11-17, historical):** The
reported d=10 unbounded-c construction had not yet been located; the d=10
models identified then were one- or two-flavour demonstrations. The
markerless contribution above supersedes this source-search conclusion.
Existing smaller candidates include d=11 frozen and d=13 mobile bonds
with certified finite members, and a d=12 mobile proposal with additional
heuristic dynamics. These are distinct from the recent d=20/22 all-n
program. The new round prioritizes <=10-state mechanisms and transfer
to the smallest concrete alphabet while preserving all CFT assumptions.
See the [source ledger and active assignments](contrib/research_coordination_20261116/SMALLEST_D_FOLLOWUP_20261117.md).

**V=4 size-scaling follow-up (2026-11-17):** The
[bulk-phase and entropy study](contrib/v4_bulk_diagnosis_20261117/README.md)
reaches L=144 in both parities and reproduces the L=96 state independently.
The latest interrupted chi=320 entropy slope is 2.57 for L=96 to 144,
versus 2.82 for L=48 to 96; size/window drift and incomplete accuracy
checks prevent identifying a thermodynamic c. Density/bond and secondary
valley-channel diagnostics do not establish persistent bulk order.
The user requested a battery-saving stop; all jobs are stopped and
checkpoints are preserved.

**Mobile length-selection test (2026-11-17):** A
[tail-sensitive comparison](contrib/mobile_length_selection_20261117/README.md)
reduces legal vacuum length words to explicit scalar vacancy problems at
J/t<=1/80. For every length, the second-order clock correction favors
the intended word in the weak marker-changing test. The scalar marker
response does not cancel; uniform selection and retuning remain unresolved.
Bounded checks retain the frozen selector and need no new penalty or
smaller J/t. The full construction's coefficient hierarchy remains.

**Two-term qubit diagnostic (2026-11-17):** The
[line `2XY+gZZ`](contrib/two_term_qubit_phase_20261117/README.md)
has an exact obstruction to the standard boost third charge, including
all range-two repairs, and a controlled large-|g| bulk-gap exclusion
with a non-explicit threshold. Near the flat point, a slow hopping is
generated at first order; this is not the tuned cancellation line.
Sector-resolved rings through L=20 and an audit of existing L=128 data
favor c=1 at g=1/2,1; stronger-coupling controls favor density order.
No all-g central-charge bound or new counterexample is established.

**Simpler qubit candidates (2026-11-17):** A
[46-run calibrated campaign](contrib/simple_qubit_candidates_20261117/README.md)
closes the missing local chi check for the old V=4 point at L=96.
Its common-basis form `(9/16)XY+(7/16)YX+ZZ-(7/4)Z` has four
monomials and coefficient spread 4 (six monomials and spread 28 in
the original basis). It is the simplest retained follow-up lead;
strong density modulation and unresolved size scaling prevent a CFT
claim. V=8 reaches L=128,chi=256 but still fails the chi-step gate;
small parameter probes show substantial density/sector changes.
No robust c>1 phase or new counterexample is established.

**Mobile vacuum-clock control (2026-11-17):** A
[coherent-gradient comparison](contrib/mobile_vacuum_clock_20261117/README.md)
proves a positive compressed clock-pin cost and bounds vacuum-clock
excitation weight with the terminal selector and vacancy dynamics retained.
Its sufficient J/t<=1/80 allowance is independent of the charged mass D
and cell length, in the stated tilted ensemble and its canonical
minimizing blocks. A parameter audit distinguishes genuine simplifications
from proof-induced hierarchies. Dressed dynamical matching, wrong-length
and illegal-sector selection, and the CFT remain open.

**Minimal qubit terms (2026-11-17):** The
[four-face classification](contrib/qubit_minimal_terms_20261117/README.md)
of `a(XX-YY)+b(XY-YX)+cZZ+hZ` excludes central charge above one when
pairing or interaction is absent, using the XXZ and free classifications.
The staggered-field and zero-field interacting faces retain explicit
unresolved regions. The latter is an asymmetric Ashkin–Teller model,
with a controlled C=1 neighbourhood of its XXZ limits. A common rotation
reduces six Pauli monomials to four but leaves all four independent
couplings; it cannot remove the reference candidate's invariant hierarchy.
The two-term subfamily `2a XY+cZZ` is a simple residual test, not a newly
established counterexample.

**First d=22 mobile auxiliary estimate (2026-11-16):** The
[mobile matching contribution](contrib/small_d_mobile_matching_20261116/README.md)
derives the literal contact-weighted representation and proves a positive
charged-rail pin cost with vacancy dynamics retained, in a positive
length-tilt ensemble and its canonical minimizing blocks. It also bounds
charged admixture at nonzero port and handles physical cyclic gluing with
one vanishing global energy allowance. Vacuum-clock matching, arbitrary
prescribed densities, wrong-length selection, and the interacting CFT
remain open; a separate canonical scalar-impurity lemma is not full matching.

**Exact qubit resonance hierarchy (2026-11-16):** The
[all-period parent analysis](contrib/qubit_resonance_hierarchy_20261116/README.md)
constructs exact nearest-neighbour qubit parents for every resonance
period m. On compatible rings their ground space has dimension L+1,
its singular limit is a span of single-residue Dicke multiplets, and
their gap has an explicit O(L^-2) upper bound. Exact local calculations
and larger complete-space checks find eighth-order motion at m=5,6.
These exact parents do not supply m positive-velocity conformal modes;
the general detuned phase problem and uniform local bound remain open.

**Resonant qubit detunings (2026-11-16):** The
[many-body and scattering analysis](contrib/qubit_resonant_liquids_20261116/README.md)
shows that the exceptional leading models sit at ferromagnetic density
thresholds. Their natural pair-amplitude detuning also creates particles
and retains a positive parent structure; treating it as repulsion alone
is incorrect. A separate hopping detuning mixes residue sectors and
passes the stated few-particle tests. An explicit fixed microscopic test
point is supplied, but its interacting infrared phase remains unverified.

**Coordinated follow-up (2026-11-16):** An
[independent effective-model check](contrib/research_coordination_20261116/README.md)
shows that the new three-/four-valley qubit branches, at their specified
sixth-order threshold, are positive configuration-graph Laplacians with
ferromagnetic single-residue sectors and exactly L+1 ground states on
compatible rings. Their band minima alone do not supply three or four
conformal modes. Follow-up work studies admissible detunings, the general
resonance hierarchy, and mobile matching for the d=22 family. The original
large-chi qubit convergence test remains blocked on production resources.

**Qubit cancellation classification (2026-11-16):** The
[new perturbative classification](contrib/qubit_cancellation_classification_20261116/README.md)
identifies both dilute and dense classical cancellation branches, gives
their complete fourth-order operators, and finds a sixth-order obstruction
to scalar cancellation on the dense branch. Two exceptional dilute paths
produce distance-four and distance-three hopping; the earlier parent paths
instead extend the blockade through distance eight. Exact symbolic and
full-sector finite-ring checks are supplied. Neither a two-valley argument
nor parent tangency proves a finite bound on total interacting c; universal
local boundedness remains open.

**Frozen d=20 defect stability (2026-11-16):** An
[exact cutting comparison](contrib/small_d_defect_stability_20261116/README.md)
now selects the all-length-n ground sectors at strictly nonzero port
coupling on rings divisible by n+1. It gives a positive cost per wrong
cell and per forbidden bond, uniformly in volume and wrong-cell length,
with finite grammar penalty and the required positive readout. No
chemical-potential retuning is needed. The interacting CFT identification,
incompatible-ring qualification, and d=22 mobile matching and defect
questions remain; this is not an unbounded-c counterexample.

**Smaller fixed-d candidate (2026-11-16):** The
[new 22-state mobile family](contrib/smaller_d_unbounded_20261116/README.md)
has an explicit bond, an exact all-length cell selector and an arbitrary-n
leading-band argument, with proposed c=n+1. Its frozen precursor has d=20.
The infrared conclusion remains conditional on mobile ground-sector stability,
vacancy matching and multichannel conformal control; it has not had
independent review. The note distinguishes this unbounded-c target from
the smaller finite-c candidates and records the
[compactness argument and its missing hypothesis](contrib/smaller_d_unbounded_20261116/COMPACTNESS.md).

**Small-d continuation (2026-11-16):** The
[stacking and vacancy note](contrib/small_d_lifting_20261116/README.md)
proposes attainable c=2,3,4 at d=2,3,4 while leaving their maximal values
open. Two onsite copies give d=4, c=4 exactly if the qubit c=2 seed is
confirmed. An explicit d=3 vacancy lift targets c=3, with additional
Gaussian stability and finite-coupling matching hypotheses. Exact closure
proofs, an entropy identity, and finite-matrix checks are supplied; no new
thermodynamic central charge is measured.

**Proposed revised qubit conjecture (2026-11-16):** For a fixed Hermitian
nearest-neighbour qubit bond repeated with **one-site translation symmetry**,
the [new working conjecture](contrib/d2_updated_conjecture_20261116/README.md)
is total unitary conformal central charge c≤2, conjecturally sharp but now
explicitly assessed as a tentative upper-bound guess. This is
not a proved upper bound or a promotion of the numerical c≈2 candidate to
an established counterexample. The note gives the exact microscopic-U(1)
reduction to c≤1, checks the two-valley mechanism, and tests repeated
blockade tuning, including two exact nonconformal frustration-free families.
Additional composite or neutral critical sectors remain the principal gap
in a general argument.

**Independent continuation (2026-11-09):** [new checked contributions](contrib/counterexample_push_20261109/README.md)
give an exact binding boundary and pair wavefunction for the Fibonacci
qubit candidate, identify its leading omitted number-changing term, and
improve the d=13 legal projected pin bound from J/4096 to J/135.
The qubit result agrees with C175's new two-body obstruction; its
finite-density c remains open. The [route audit](contrib/route_audit_20261109/README.md)
ranks d=13 as the strongest mature construction and keeps its
model-specific matching and relative defect hypotheses explicit.
The same continuation supplies a [new transverse-DM qubit family](contrib/counterexample_push_20261109/FACILITATED_ROUTE.md)
whose complete effective Hamiltonian weakens binding perturbatively;
the apparent small-ring crossing is strongly size dependent and is not
a counterexample claim.

**d=2 update (2026-11-10):** More accurate full-spin reruns recover c≈2
at the phase-1.53 point through L=256. The earlier unconverged L=128
downturn does not reproduce. Six tuned bonds give ground-branch entropy
estimates 1.94–2.03 from 128 to 256 sites, with calibrated controls.
The L=256 entropy still changes too much with bond dimension, so no
registered pass or established d=2 counterexample is declared. L=512
runs are recorded as in progress; infinite-chain microscopic convergence
remains open. See [current status](STATUS.md) and the [ground-sector
selection audit](contrib/d2_status_audit_20261110/README.md). The transverse-DM
points previously tested read near c=1.

**Conditional large-d constructions (scope clarified 2026-11-17).** Explicit
chains predict c=14 at d=13 and unbounded c at fixed d=4464, subject to
their stated matching, sector-selection and infrared inputs. The finite
certificates do not by themselves establish those microscopic CFTs.
The proposed replacement, Σ_a c_a v_a ≤ π M log₂ d, survives the recorded
tests and is proved in the stated classes, but remains open in general.
See [STATUS.md](STATUS.md) and each construction's hypothesis ledger.

A d=24 conditional counterexample candidate is recorded with independent AI reviews and a [fresh adversarial audit](contrib/counterexample_audit_20261028/README.md). The earlier verified-counterexample wording is withdrawn: finite-J matching and a uniform dressed-defect comparison remain open. The working standard is theoretical physics (standard RG, controlled perturbation theory, calibrated numerics), not mathematical rigor.

- **Small-alphabet candidates with finite certificates (2026-10-29): d = 11 frozen with predicted c = 11, and d = 13 mobile with predicted c = 14.** Explicit route-A bonds with certified finite members, refereed, and red-teamed for d = 13; both rest on the same inherited infrared inputs as the d = 24 candidate. Full bonds and status in [contrib/small_d_note/](contrib/small_d_note/note.pdf) and STATUS.md. The logical vacancy model gives calibrated c = 2, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3 (REPORT §9.34, §9.55, §9.67). The smaller d=12 mobile proposal is separately recorded with additional heuristic dynamics.

- **Unbounded c at fixed d0 = 4464 (the C39/C44/C45 wrapper).** c = n + 1 for every n, conditional on the Wilsonian input S. Several pieces of S are now proved and refereed (massive matching C48, ensemble transfer C58, c from a Gaussian limit C51); the open core is the constructive RG (Lemma N and lemma TP). The full microscopic demonstrator flows to its matched XXZ Luttinger liquid numerically (REPORT §9.58). Expert summary: [contrib/expert_note/note.pdf](contrib/expert_note/note.pdf).

- **d = 2.** No c > 1 found in about 103,000 sampled qubit chains of the residual strata, on the staggered branch to L = 256, or inside floating phases (REPORT §9.37, §9.40, §9.62); the analytic covering argument is not complete.

- **d=24, predicted c=24: a one-vacancy candidate.** Add one mobile vacancy to the explicit 23-state chain, with its original bond acting on occupied pairs. The local phase-model representation, certified 23-band reference, nonzero-coupling sector bounds and physical momentum rules support 23 payload bosons plus one charge boson, conditional on finite-coupling matching and defect stability. Published local perturbation bounds control the original cell machinery. The conclusion uses conventional Gaussian phase stability; no microscopic c=24 measurement or rigorous CFT-limit theorem is claimed. Read [the construction and independent reviews](contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md).

- **Independent vacancy continuation.** A proved density window removes exact half-filling tuning from the proposed d=24 construction. New reviews check the local perturbation argument and C39's all-n finite input, including a new exact four-mode certificate. A full microscopic d=3 control has entropy and independently calibrated energy scaling consistent with c=2; this tests the vacancy mechanism, not the c=24 chain. See [the proofs, reviews and numerical data](contrib/vacancy_continuation/README.md).

- **Small d: no counterexample.** Nothing above d − 1 at d = 2, in symmetric families at d = 3 to 7, or at the best d = 4 target. The d = 3 hint failed its own follow-up (C12); an independent check remains pending in the recorded status.
- **New d = 2 proof reductions.** Arbitrary uniform fields in the symmetric ferromagnetic rank-one class give at most c = 1/2, using a recent rigorous gap theorem and the exact Ising solution. An exhaustive classification of additive on-site spin charges leaves one additional two-parameter staggered-field family; exact gap arguments exclude both outer regions. The full qubit conjecture remains open. See [the proofs and remaining region](contrib/d2_gap_reduction/README.md).
- **No general proof.** Every lattice quantity we know that bounds c grows like (a length scale) × log d. Proofs exist for quasi-free chains (c ≤ log₂ d), SU(d)-symmetric bonds, and some classes of IR theories, and there are conditional bounds.
- **Large d: a candidate physics-level refutation.** The new sparse species-graph family has predicted c=(p+1)(p²+p+1)+1 and microscopic d=O(p² log p), so c/d diverges. Exact finite constructions now supply the ports, reversible circuit, O(q log q) controller and initialized positive core; the existing mobile wrapper applies. The conformal conclusion uses stated standard massive-sector integration and Gaussian RG inputs. It is not a rigorous thermodynamic theorem or a measured numerical violation. Read [the complete candidate argument](contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md) and [the infrared audit](contrib/research_continuation/IR_AUDIT.md).
- **Unbounded c at one fixed d remains conditional.** C39/C42 and the C44 wrapper now supply a finite spatial construction at d₀=4464. C45 gives a coupling schedule conditional on a stated Wilsonian input S; C46 proves partial rigorous results but leaves the full infrared step open. Read [the referee](contrib/referee_C44_C45/README.md) and [C46](codex/work_C46/RESULT.md). The sparse first-order family with growing d is a separate route.

[STATUS.md](STATUS.md) has the list of where help is most useful, and [STATUS_DETAIL.md](STATUS_DETAIL.md) the details.

The [direct-attack continuation](contrib/direct_attack_20261028/README.md) supplies ultraviolet matching bounds, an exact exclusion of a hidden free-fermion qubit mechanism, and an independently AI-audited thermal massive-matching theorem allowing infinite excited cell spaces. The retarded multichannel RG and central-charge conclusion remain open.

The [physical-ensemble continuation](contrib/thermal_ensemble_20261028/README.md)
proves interacting physical-to-rooted pressure transfer, an all-length
pressure-zero theorem and a length-tilted thermal reference comparison.
Its [positive completion and matching application](contrib/thermal_ensemble_20261028/MOBILE_MATCHING_APPLICATION.md)
repair the mobile bond on illegal words, retain d0=4464 and the complete
legal dynamics, and supply the full-space massive and ensemble inputs.
The parallel [frozen-chain continuation](contrib/decisive_attack_20261028/README.md)
controls the frozen-defect pressure while retaining the incompatible-ring
qualification. The logical infrared conclusion remains open in both.

**Reversible digit compiler (2026-10-28):**
[digit_compiler](contrib/digit_compiler/README.md) gives fixed tables for all
generator arithmetic and selected-array operations, then compiles them to
one physical reversible head. Actual arithmetic runs restore every work and
history cell. Full-space local permutation gates also pass autonomous
compute/reset cycles. Whole-interpreter storage and dispatch, initialized
Hamiltonian composition, and the final violating dimension remain open.

**Digit compilation and infrared progress (2026-10-28):**
[digit_compilation](contrib/digit_compilation/README.md) supplies an explicit
96-entry parameter-reader table, eleven finite local arithmetic kernels,
and checked execution of eight arithmetic opcodes inside the complete
21-mode circuit generator. All 287,798 abstract gates match. A separate
matrix theorem gives the explicit Gaussian stability budget delta<5/13,
including density-current couplings and unequal velocities. The remaining
memory/controller composition and microscopic low-energy reduction are
still required; no new violating dimension or central-charge measurement
is claimed.

**Microscopic implementation progress (2026-10-28):**
[executable circuit and controller work](contrib/controller_instantiation/README.md)
now supplies a full sparse-encoding gate generator, explicit reversible
dispatcher and length-counter tables, a repaired input-format check, and
a saved d=240 one-mode test Hamiltonian. The autonomous program generator,
verified program initialization and endpoint readout still need to be
composed into the high-channel core. The saved fixture is not a c-d violation.

The [autonomous-controller continuation](contrib/autonomous_controller/README.md)
now composes blank program generation, coherent execution, endpoint readout
and complete reset in one explicit local machine for a demonstration program.
Its initialized clock is also instantiated. A separate fixed 3,115-instruction
word program reproduces the full sparse-circuit generator; digit-level local
compilation and high-channel composition remain. Two explicit parameter-query
gates remove the need to decode the chosen prime from cell length. None of
these demonstration dimensions is a newly claimed violating d.

Additional derivations from the chat "Prove c d conjecture" are collected in
[contrib/chat_followup/](contrib/chat_followup/README.md): a species-torus
selection rule for the encoded target, a massive-auxiliary construction,
exact extensions of controlled classes, and a detailed qubit catalog. The
imported proposal remains conditional; the [new continuation](contrib/research_continuation/README.md) gives a different sparse-graph construction. Also see [C14's progress
log](codex/work_C14/PROGRESS.md) for its reported physics-level counterexample
to the proposed velocity-weighted Casimir bound; that does not refute c-d.
