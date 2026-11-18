# Where the c–d conjecture stands

Last scope and priority update 2026-11-18; the older numerical overview below
is dated 2026-11-12. This page lists only the main results. [REPORT.md](REPORT.md)
is the chronological record, and [STATUS_DETAIL.md](STATUS_DETAIL.md) is an
earlier long version. The standard is theoretical physics, and nothing has
had independent human review. Conditional infrared predictions are not
established microscopic counterexamples.

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

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242)) conjecture c ≤ d − 1 for translation-invariant nearest-neighbour chains on (ℂ^d)^⊗L whose low-energy spectrum is a unitary CFT.

## In short

1. Explicit conditional candidates exist at d=9 (frozen cells, predicted
   c=n) and d=10 (mobile cells, predicted c=n+1), with 2n simple Fermi
   points certified at n=9,10,11,12,14,16,20, alongside the earlier d=11
   and d=13 bonds and proposed
   unbounded-c families at fixed d. Their finite certificates and partial
   matching results do not establish the full microscopic CFT conclusions.
   A modified d=9 clock/readout now has an arbitrary-n leading quadratic
   band proof, finite-order volume-independent local bounds and exact
   frozen ground-sector selection on compatible rings. All-scale
   interacting matching and a CFT theorem remain open.
2. At d = 2 there is a candidate counterexample, an explicit nearest-neighbour qubit chain that reads c = 2.00 on open chains up to L = 256, converged in bond dimension, and 1.95 to 1.98 from L = 256 to 512, passing the same test there, with a mechanism and an effective model that give the same answer. It remains subject to the convergence and infrared-stability checks listed below, including charge locking; these finite-size readings do not settle its thermodynamic c.
3. A replacement bound, Σ_a c_a v_a ≤ π M log₂ d, holds in every solvable family and is proved in several classes.

## Conditional large-d constructions

1. Proposed unbounded c at fixed local dimension d₀=4464, with predicted
   c=n+1 for every n (C39, C42, C44, C45). The finite inputs are reported
   as proved and refereed. The microscopic-to-Luttinger-liquid conclusion
   remains conditional on the stated Wilsonian input.
2. Small-alphabet candidates with certified finite members include d=11
   (frozen cells, predicted c=11) and d=13 (mobile cells, predicted c=14),
   written out in [contrib/small_d_note/](contrib/small_d_note/note.pdf).
   Their finite ingredients have AI referee checks, including a red-team
   review of d=13. They retain the stated infrared inputs; the
   [independent audit](contrib/counterexample_push_20261109/README.md)
   leaves the d=13 model-specific hypotheses open. A
   [d=12 mobile proposal](contrib/mobile_below_13/README.md) has additional
   heuristic dynamical steps and is not an established smaller counterexample.
   The [markerless chains](contrib/markerless_d9_d10/README.md) lower these
   to d=9 frozen and d=10 mobile. Every head moves by one rule, so the legal
   sector and the cell pattern are exact invariants. This removes finite-K
   leakage and the need for vacancy colours; a penalty can still raise
   competing illegal sectors. The band data are certified with C84's pipeline at
   n=9, 10, 11, 12, 14, 16 and 20 (predicted c up to 21 at d=10), not for
   all n. An AI referee confirmed the frozen chain and the certificate and
   re-checked the corrected mobile bond. The same infrared inputs apply,
   and J/t for thermodynamic mobile selection and dressed bands is not
   quantified. Subsequent [finite mobile certificates](contrib/markerless_mobile_control_20261118/README.md)
   compare (10,10) with (9,11), (20), and (6,7,7) at N=20,V=2 and
   J/t=1/25. Moving only the conserved charge field outside the contact
   embedding exactly removes its bit-position error; other mobile errors
   and the severe coefficient hierarchy remain.
   Separately, the [all-n modified family](contrib/markerless_all_n_20261118/README.md)
   proves strict moment negativity and 2n simple, pair-sum-separated roots
   for every finite n at physical d=9. It uses a biased clock and an L1
   field, so the mobile charge-field theorem cannot be transferred to it
   without a new argument. The [uniform follow-up](contrib/markerless_uniform_matching_20261118/README.md)
   supplies finite-order local interaction bounds and exact nonzero-port
   frozen ground-sector selection on compatible rings; it retains residual
   auxiliary coupling and does not establish the CFT. The distinct
   [mobile half-occupancy analysis](contrib/markerless_thermodynamic_selector_20261118/README.md)
   proves one uniform susceptibility channel, but not the full selector.
3. Numerics agree. The logical vacancy model gives c = 2.0, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3, and the projected d₀ = 4464 chain gives c = 3.00 ± 0.02 at n = 2.

Not established: a rigorous lattice-to-CFT theorem for any counterexample.

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

## d = 2: a candidate counterexample with c = 2

The chain is

H = Σ_j [ (X_j X_{j+1} − Y_j Y_{j+1})/2 + (X_j Y_{j+1} − Y_j X_{j+1})/160 + 10 Z_j Z_{j+1} ] + h_z Σ_j Z_j, with h_z = −19.976047,

together with a family of neighbours of it. One of them, proposed by Pace, adds a small XX + YY term (hopping phase 1.53). Off the symmetric line a momentum mismatch may cut off charge locking in the estimate of open item 5, making the phase-1.53 bond the stronger candidate in that analysis. This is not a matched thermodynamic stability result.

The strong ZZ term, with the field near its blockade value, leaves a dilute gas of blockaded particles. Their direct hop (the XY − YX term) is tuned down to the size of the next-nearest hop that the pairing term XX − YY generates at second order. The particles then move on a zigzag ladder whose band has two valleys, and the proposed infrared theory is a two-component Luttinger liquid with c = 2.

The spread of couplings is one small parameter ε = Δ/V (pair amplitude Δ = J_x − J_y, blockade energy V = 4J_z), not three tunings. The next-nearest hop exists only at second order, |t₂| = Δ²/2V, and the direct hop 2D is set equal to it (D = 1/160 above), giving the ladder εΔ : Δ : Δ/ε. The field sits at the blockade value plus a remainder of order t₂ that sets the density. As ε → 0 only the hopping phase, the pairing strength and the density matter, and the readings agree at V/Δ = 40 and 80. How large ε can be is open.

Evidence, all numerical:

- The qubit chain itself. Parity-resolved DMRG on open chains of L = 64, 128 and 256, with no projection, converged in bond dimension up to χ = 1024, reads c = 2.00 for the bond above and for the phase-1.53 bond from the entropy difference between L = 128 and 256, which passes the pre-registered test at both. Four further tuned bonds (two densities, V/Δ = 40 and 80) read 1.94 to 2.03 at χ = 768, where the reading is still rising with χ. The c = 1 and c = 2 controls calibrate both size pairs, and independent CPU and GPU codes agree (REPORT §9.178, §9.186, §9.192).
- The effective model. The second-order model of the blockaded particles reads c = 1.99 to 2.00 on the infinite chain at three densities, out to correlation lengths of 500 to 600 sites. Its open chains read 1.82, 1.88 and 1.92 over the same three size doublings, approaching 2 from below while the qubit chain approaches from above (REPORT §9.163, §9.179, §9.204).
- The mechanism. One particle from each valley binds into a pair in the dilute limit, which alone would give a pair liquid with c = 1. At the candidate density the inter-valley backscattering responsible is irrelevant if the pseudospin Luttinger parameter K_s exceeds 1. K_s reads 1.2 to 1.5 in the effective model and 1.4 to 1.6 on the qubit chain at L = 256, where the pair exponent and K_c also match the effective model (REPORT §9.160, §9.161, §9.182).

Open, before this is called established:

1. Convergence at L = 512. The third size reads c = 1.95 and 1.98 from 256 to 512 at bond dimension 2048, where an exact c = 2 control reads 1.95 at the same sizes, and passes its registered test at both bonds. The entropy at L = 512 still moves by 0.010 to 0.012 from bond dimension 1536 to 2048, just above the registered 0.01 (REPORT §9.199, §9.203).
2. Entanglement scaling on the infinite qubit chain. Two attempts did not converge at the candidate, because the finite-χ state breaks particle number (REPORT §9.183, §9.187).
3. The pseudospin gap. In the effective model the pair gap is seen directly at low density, on finite and converged infinite chains, and closes near n ≈ 0.07, against the candidate's 0.19. The registered test puts the candidate on the gapless side in substance, and falls short of the letter only because no single density shows the gap above three errors on both chains at once (REPORT §9.188, §9.204).
4. One unexplained reading. An older-solver run at V/Δ = 8 reads about 1.4 on the infinite chain. (The c = 1 control of the infinite-chain calculation, the other loose end, now reads 0.99, REPORT §9.201.)
5. Charge locking. The qubit chain conserves only parity, so particle number is emergent. A number-changing term O₄ (two inter-valley pairs, scaling dimension about 1.25) is relevant when its momentum matches, and on the symmetric line it does. Its coupling is fourth order in the pairing, about (16/9) C₄ ε² with C₄ not computed and not forbidden by any known symmetry. For C₄ of order one this gives a locking length of several thousand sites, beyond which the charge mode on the symmetric line would gap and leave c ≤ 1. At the phase-1.53 bond the coupling oscillates with a period of about 420 sites, more than ten times shorter, and the liquid survives in this estimate. This is an order-of-magnitude estimate, not a calculation, and chains of L ≤ 512 cannot test it (REPORT §9.161, §9.162; [collaborator note](contrib/collaborator_notes/qubit_c2/notes_qubit_c2.pdf), Sec. 6).

Pace's earlier runs of the phase-1.53 bond read a downturn from L = 64 to 128 on states that were not fully optimized ([record](contrib/d2_tenpy_20261109/RESULT.md), [rational bond](contrib/d2_dmrg_priority_20261109/SIMPLE_PARITY_CANDIDATE.md)). The converged rerun lies lower in energy and reads 2.2, so the downturn does not reproduce (REPORT §9.167). Their [transverse-DM deformation](contrib/counterexample_push_20261109/FACILITATED_ROUTE.md) reads c ≈ 1, and their [symmetry audit](contrib/d2_priority_audit_20261109/UNIFORM_Y_SOURCE.md) allows a small uniform Y source.

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

## Elsewhere at d = 2

Away from this corner no chain with c > 1 was found in about 103,000 sampled chains, targeted searches, or the literature (REPORT §9.110 to §9.150, refereed only where noted). A recalibrated screen of the uniform sample bounds any open c > 1 region at about 0.1 percent of bond space.

- c ≤ 1 is derived at the physics standard at generic points of every symmetric class: U(1), free fermion, XYZ, cyclic, parity (outside the corner above), staggered charge, and the singular strata. On the time-reversal strata T and K_I it holds given a one-species assumption, and the sheet K_II is open. Isolated multicritical points on the strata have not been inventoried.
- The XZ point, the only bond with a second decoupled gapless species, gives no c > 1. Its induced interaction has the wrong sign for the c = 3/2 phase. Its neighbourhood, which reads c ≈ 1.3 on open chains to L = 256, is a gapped Néel bulk seen through a slowly relaxing cat, confirmed at L = 512.
- Close to the ferromagnetic Heisenberg point every conformal bond has c ≤ 1, given standard semiclassics (refereed, and checked on 48 lattice bonds).
- Where a symmetry makes two soliton valleys degenerate, the two-valley gas does not mix, so c = 1.

Open here:

1. Generic bonds with no symmetry, and attractive liquids whose U(1)-breaking is as large as the exchange. The symmetric results do not propagate into this region (REPORT §9.127). It includes single-species c > 1 theories such as Z₅ parafermions and three-state Potts order next to a floating liquid (c = 9/5). A structured search for exact maps into the nearest-neighbour class found no c > 1 model (REPORT §9.134).
2. Two mechanisms for a second gapless sector at strong U(1)-breaking, a staggered Ising transition inside an incommensurate liquid and accidental two-valley soliton onsets. Searches found neither. A numerical example of the first now exists: a parity-only qubit chain whose staggered chirality orders at h_I = 1.7466 inside an incommensurate liquid, reading c = 1.53 to 1.59 against 1.42 to 1.52 for c = 3/2 controls (REPORT §9.208, [contrib](contrib/d2_xy_yx_chirality_imps_20261118/README.md)). The c ≈ 2 readings just below h_I are unresolved.

## A replacement bound

The proposed replacement is A = Σ_a c_a v_a ≤ π M log₂ d, with M the gauge-optimized bond half-spread. It is proved for quasi-free chains (A ≤ 2M, sharp), for free-fermion and classical chains, and for WZW sectors coupled to on-site symmetries. Every solvable family satisfies it, and an adversarial search found no violation. The general case reduces to C(d) = sup L(E_{2L} − 2E_L)/M < ∞.

At d = 2 on the U(1) class (XXZ plus DM along z plus field) C(2) is now finite on every bond, every row proved and refereed by AI agents (2026-11-18). The ferromagnetic strip −1 < J_z/J ≤ −cos(1/64) is closed by a uniform counting bound on the Bethe roots, with a universal constant that is not yet explicit, and the easy-axis side (C = 10781) passed its referee. The constants are far from sharp, about 10⁴⁸ on −1 + 1.2·10⁻⁴ ≤ J_z/J ≤ 0 (REPORT §9.94 to §9.109, §9.207).

At general d, three theorems hold (REPORT §9.207). For every chain, (L e₀ − E_L)/M → 0 uniformly in the bond at fixed d, a qualitative result with no content at reachable sizes that does not bound A. Two-species exchange chains with any internal register obey A ≤ 2M at zero interaction, and critical ferromagnetic Potts chains at Q = 3, 4 obey a crossover-free bound inherited from XXZ, which covers the neutral c = 4/5 sector. A general finite C(d) is still open.

## Small d other than 2

At d = 3 to 7, nothing above d − 1 in symmetric families. A Casimir-type bound Σ c v ≤ 2π‖h‖ is false (SU(N) chains with rank-2 symmetric representations, from N = 7).

## Where help is most useful

1. An independent reproduction of the d = 2 candidate. The Hamiltonian above, any DMRG code, L = 128 and 256 at bond dimension 768 or more.
2. Constructive RG for the Wilsonian step behind the large-d counterexamples.
3. Independent review of the markerless d = 9 and d = 10 bonds, especially
   the corrected mobile selector, its quantitative J/t budget and all-n
   band separation.
4. A general proof of the replacement bound (a finite C(d) at any d ≥ 2 beyond the classes above), or an explicit constant on the d = 2 strip.

All of this was produced with AI assistance (Claude Opus 5.5 and Claude Fable 5 coordinating, with gpt-6-astra runs).
