# C25 progress

2026-10-28. G3 only. No earlier local progress existed.

[Assumed inputs] Read project STATE, route-A background, assembly README,
STATUS, C18 Section 8 and its clock code, C8 tagged-vacancy, helix and
selector notes. Read C21, C22 and C8 progress. C26, C27 and the three
Claude contribution progress files were absent. No network or delegation.
The project STATE is outside the permitted write folder and will not be
modified. This file carries the local handoff.

[Proved, derivation being written] The C18 clock can be replaced termwise
by positive propagation Laplacians. The uniform history, Bell pairs and
right input pin then give the same two-dimensional exact code. Wrapping
its four occupied sites gives physical lengths 5 and 6, not 4 and 5.
Vacancies gate every original two-site term, including Bell checks.
The half-onsite turn in the original repeated bond is gated by the
average of its two occupied-neighbor projectors.

[Proved, derivation being written] Bell syndromes and gauged input bits
reduce each port-accessible charge-clock block to 150 or 450 dimensions.
The actual helix compensator creates only zero and one flipped spin in
the untwisted ferromagnetic frame, so its full accessible spectrum fits
in 6 or 7 dimensions per active track. This is an invariant reduction,
not a projection that discards virtual excitations.

Next: assemble the exact sectors, contract their two-cell resolvents,
and compare all four adjacent-length combinations with m_p=1/r_p.
Also evaluate the mixed marker-port coefficient using the syndrome
selection rule and verify positivity of the proposed G1 check.

## Exact model assembled

[Proved] In the legal clock and ordered-tag sector the wrapped Hamiltonian
is H_ch-e_p plus the core's three positive bond terms gated by adjacency
projectors D_j. Bell checks are gated too. The uniform clock code remains
exact for every charge configuration. The occupation filter stays in its
exact product ground state throughout a port excursion and can be omitted
from the reachable sector without approximation.

[Proved] The coefficient at order u lambda is zero. Marker transport
preserves every ordered Bell syndrome, whereas a single port exchange
creates a nonzero syndrome on two cells. The unperturbed resolvent preserves
these labels. This applies to both orderings and all marker destinations.
It does not prove vanishing for a syndrome-mixing G1 implementation.

[Numerical, ongoing] The charge matrices have dimensions 5 and 15.
At w/U=0.16 the two-vacancy tunneling gap is about 1e-7 U. Raw diagonalization
mixed reflection eigenvectors and triggered the reflection check. The
implementation now diagonalizes reflection sectors separately and solves
the exact charge-ground / charge-excited block decoupling before source
spectral decomposition. The correction retains both blocks and all weights.
No failed numerical result is used in the report.

[Proved] For fixed w>0 and source scale s tending to zero, the effective
core is sum_j <D_j> h_j, not the original ungated core. The compensator
removes endpoint factors but cannot remove this finite-w change of spectrum.
The scan compares the exact result to both that limit and the original core.

## Completed calculation and audit

[Proved] RESULT.md now gives the full source matrices, both endpoint
measures, the exact two-cell resolvent contraction, a finite Schur bound,
the missing-bond estimate and the mixed-term selection rule. The limits
are ordered: first small nonzero w/U, then small s/Delta_wrap. The latter
scale includes the two-vacancy tunneling gap of order w^8/U^7. At fixed
w, the limiting core has mean-adjacency factors <D_j>, so the existing
compensator is not an exact second-order equalizer.

[Numerical] composition.py completed 27 parameter points and 108 ordered
length-pair bonds in about 36 seconds. Source blocks have dimensions 900
and 3150 including the accessible helix sectors. At w/U=0.06 and
s/Delta_wrap=0.001 the relative traceless length mismatch is 0.0013266;
its limiting value at that w is 0.0013359. The full scalar corrections
are larger and are explicitly retained in the output and report.

[Numerical] checks.py passed independent raw source inversions at
dimension 450, full helix checks at dimensions 32 and 64, literal Bell
source checks at dimension 16, and the G1 demonstrator at dimension 32.
Largest scaled raw-inverse discrepancy is 1.72e-9. All claims from the
smallest charge gaps are floating-point results, not certified intervals.

[Proved] The syndrome-mixing G1 local check is positive for D0,D1>0,
kappa^2<=D0 D1, and has the desired complete local kernel under strict
inequality. It preserves termwise annihilation of idle-Bell work histories.
This confirms compatibility with the wrapper hypothesis, not a compiler.
The u lambda coefficient is zero only for the syndrome-conserving C18
cell calculated here.

[Assumed remaining work] No hopping-producing G1 cell was substituted,
no whole repeated bond or all-ring finite-coupling spectrum was computed,
and no infrared conclusion was established. The obstruction to exact
finite-parameter equalization restricts this implementation and does not
kill route A. The named G3 next calculation is finished.

Read the newly available C26 and C27 progress at the final audit. Their
fourth-order and gluing work is disjoint from this calculation. The three
Claude PROGRESS files remain absent. No external file was modified.
RESULT.md, composition.py, checks.py and their local JSON outputs are the
handoff. Local STATE.md records the result because project STATE is
outside the user-authorized write directory.
