# C81 handoff

2026-10-28. The task is complete as a conditional candidate and finite audit.
Read RESULT.md, then PROOF.md. No microscopic c=16 result is established.

The literal 16-state bond is in bond16.py, independently constructed from
C75's move table. The whole occupied bond, including half-site terms, is
embedded in occupied pairs. Exact N,Q,M survive. Physical half occupancy
is the proposed density. One vacancy means one extra local state.

Inherited C75 inputs, not recomputed: n=15 certificate, thirty simple
roots, delta2>0.0145540, delta4>1.1444e-7, uniform finite-K whole-cell
energy bound, selected-cell code/gap. Parameters and refinements are in
PROOF.md Section 2. If g decreases, K must increase to control 1/(K g^2).

New result: finite-K nonzero-port FROZEN sector selection follows by using
full cells of dimension 3*12^(l-1), including illegal strings. Marker
positions remain conserved in the frozen chain. Only boundary ports
couple full cells, with norm sqrt(2)|lambda|. The wrong-cell density
margin is gamma=2^-16-1/(K-4). A ring with wrong material length D obeys
E_wrong-E_sel >=(gamma-5sqrt(2)|lambda|)D. The thermal counterpart is
0<=f_cell/16-f15<=480 T exp[-(gamma-3sqrt(2)|lambda|)/T] at low T.
These proofs use C75's all-length bound and have not been independently
refereed. They are not mobile dressed-defect bounds.

Transferred reference results: phase-fluid K_gap=4, physical K_particle=1,
cyclic projection coefficient alpha=1/2-2/pi^2, physical dangerous-vertex
margin >3.57625e-9, free nonforward charge threshold >4.5534e-9 t,
marker exponent 1/128. The full harmonic lattice is retained.

Important matching qualification: exp(i pi Q/2) is logical parity on
the code but has eigenvalues +/-i on odd-Q auxiliary states. The d24
microscopic Z2 grading cannot be copied verbatim. A suitable full-space
representation or preliminary elimination must be checked before citing
the N23 massive-matching application for C75.

Checks passed: independent bond identity; literal/squeezed matrices;
330-state cyclic sea projection inside full 1350-state L4,N2 space;
dense ED of all 120 lowest levels with J^2 errors; two n1 cell sectors
at L5,6 of dimensions 2320,6960; scalar phase-model ED through L12;
full frozen cells of dimensions 1296 at resolved parameters. No n15 ED
or central-charge extraction. See checks.json and auxiliary_checks.json.

Remaining for the physics claim: strong-norm massive matching, uniform
low-frequency vacancy matching, mobile defect comparison, density tuning,
retarded multichannel RG, and conformal/thermal identification. The supplied
spatially weighted ultraviolet lemma applies at fixed finite K but does
not close these steps. The project STATE and all external files were left
untouched, as required.
