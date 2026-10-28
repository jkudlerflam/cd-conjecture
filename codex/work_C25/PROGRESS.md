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
