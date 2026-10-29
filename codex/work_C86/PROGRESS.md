# C86 progress

[Read] Read the simplest background, project memory, and the C75, C80,
and C84 construction notes, alphabet arguments, rejected alternatives,
and certificate source. C80 has no separate ALPHABET.md; its accounting
is in PROOF Section 8. Reused source was compared by full source diffs.

[Scope] All new scripts and outputs stay in work_C86. The project-level
STATE is outside the permitted write scope. No network or cluster work.

[Completed search] Tested a single spatial rail, reflection-related
endpoint sources, a port supported on two carrier vertices, and the
charge constraints on a single species with a data-bit conversion.
No alphabet reduction or new certified member was obtained. The
candidates fail before C84's positive-metric root stage.

[Proved] Reflection replaces the desired Gram matrix by F=GJ, with
J the flavor-reflection involution. Its Hermitian part has blocks
G_++ and -G_--. The derivation includes complex scalar carriers and
extended sources. Reflection also gives h(-k)=J h(k) J, so four or
more distinct roots force a zero unordered pair-sum margin. A common
port phase only shifts the reflection center.

[Proved, restricted] The alphabet bound a>=11 now covers the
single-token product class with matched common carrier species,
frozen on-site markers, nearest-neighbor token motion, and local
payload-independent pair ports. At most ten letters forces a scalar
path with point sources. Its spectrum is reciprocal. PROOF Section 7
states every assumption and explains why this is not a universal bound.

[Proved, computer-assisted] The C84-style 256-bit resolvent enclosures
reject 20 single-rail cases and 98 complex two-vertex-source cases.
All have negative principal determinants, with rational sign-gap
bounds. The largest product graph has dimension 121. An all-parameter
two-vertex formula is verified at 12 rational points. Four direct
Fraction inverses independently check selected Jacobi enclosures.

[Rejected checker, corrected] The first extended-source run failed
an entrywise-positivity assertion. Signed sources need not meet it.
The corrected test bounds the determinant with all interval endpoint
signs retained. rejected_checks.json records the failure.

[Numerical consistency checks] Independent heat-kernel integration
on complex graphs of sizes 4,5,6 agrees with the new Gram identity
to relative error below 1.7e-15. Reflection and onsite identities
agree below 2.5e-15. These are separate from the rational acceptances.

[Completed] RESULT.md, PROOF.md, ALPHABET.md and STATE.md contain the
outcome and scope. All four Python scripts pass, and the final audit
verifies all 118 negative determinant bounds and positive sign gaps.
No new root certificate, finite-K schedule, central-charge inference,
or independent referee is claimed. All work remains in work_C86.
