# A finite endpoint and color shuttle

Status: proved finite interface, with exhaustive local checks and finite
history checks. Conditional on prepared marked output records. It is not a
complete initialized high-channel core or a conformal-chain construction.

## 1. What this adds

The run-7 clock and clean-readout contract are inherited. The new sparse
construction in contrib/controller_instantiation supplies a circuit and a
preloaded dispatcher, but its README explicitly leaves endpoint transport
and elementary-clock readout unfinished. The following machine supplies
that interface for any four prepared output records. It uses no travel
history, no logical-value-dependent control, and no growing controller.
The implementation is independent of that dispatcher and can be called by
one. Its integration with a generator and verifier is still required.

Let n>=8. The ordinary tape has static marks

    L LI ... LP LC ... RC RP ... RI R,

with LP at a, LC at a+1, RC at b, RP at b+1, where
2<=a, a+1<b, b+1<=n-3. Unlisted sites have mark B. Each site has one
READY bit and one quantum record of dimension D. The latter can be a
q+1-state graded port register. Other inert work factors can be tensored
on, but their dimensions must then be counted.

The intended source records are LP, LC, RC and RP. These may be entangled
with each other and with all the other work registers. The machine starts
at L with READY=0 everywhere. It does not assume the source code is a
product state or retain a copy of a quantum record in its control.

## 2. The literal collision table

endpoint_shuttle.py defines 21 direction-fixed controls, one for each
segment in the following itinerary. The letter S means SWAP the quantum
records at every elementary move. The letter I means identity on records.
A flag toggle occurs on arrival at the end of a starred segment.

| Phase | Segment | Gate |
| :-- | :-- | :-- |
| 0 | L to LP | I |
| 1 | LP to L | S |
| 2 | L to LC | I |
| 3 | LC to LI | S |
| 4 | LI to RP | I |
| 5 | RP to R | S |
| 6 | R to RC | I |
| 7 | RC to RI | S |
| 8 | RI to R, toggle right READY | I |
| 9 | R to L, toggle left READY | I |
| 10 | L to R, clear right READY | I |
| 11 | R to L, clear left READY | I |
| 12 | L to RI | I |
| 13 to 20 | Reverse phases 7 to 0 | Same gate, reverse direction |

The ordinary neutral symbol is (mark,READY), giving B=18 states. At the
end of a segment its control advances modulo 21 and its specified flag
is toggled. At a mark different from its start and end, the control and
flag are unchanged. No semantic rule is specified at the starting mark.
Those local pairs never occur as arrivals in the promised trajectory.

These are 336 rules. Their domains and ranges are disjoint. For any
control k, its scan-loop range excludes its starting mark; the sole rule
entering control k from the preceding phase arrives at that mark. Its
scan-loop domain excludes its endpoint; the rule leaving that phase uses
that endpoint. Flag toggles are permutations and the static marks are
unchanged. Thus this is an injection on local pairs, not merely an
invertible trajectory.

For each static mark, fix every unused domain/range pair common to both
sets, then pair the remaining domain/range lists in increasing order.
There are equal numbers, since all prescribed entries preserve that mark.
This specifies the complete permutation C on 378 states, exported as
endpoint_shuttle_table.json. No behavior on malformed local inputs is
left undefined.

## 3. The repeated-bond implementation and its count

Use the single-head ordered-sweep interface of
contrib/exact_route/SWEEP_INTERFACE.md. This is an explicit alternative to
the two-lane packet layout of run 7. A site is either an ordinary record
or a record with one of 21 controls and a fresh/updated head bit. Its
dimension is

    w = 18(1+2*21)D = 774D.

A control has a fixed direction, so no additional direction factor is
needed. For a right-moving fresh head, G moves it from the left site to
the right, changes fresh to updated, and applies its I or SWAP record
gate. Its reverse block applies the inverse gate. For a left-moving head,
use the opposite orientation. All other two-site sectors are identity,
including the zero-head and two-head sectors. These subspaces are
orthogonal, so G is a Hermitian unitary on the entire two-site space.

The onsite F maps an updated head through C to a fresh head, and maps a
fresh head through C^{-1} to an updated head. It fixes ordinary sites.
Thus F is also a Hermitian unitary on the entire onsite space.

In an ordered left-to-right sweep, a fresh head makes only one move. A
right-moving updated head cannot match the forward rule at the next bond
or its inverse, which requires it in the right slot. A left-moving head
moves behind all remaining bonds. The subsequent onsite return sweep
performs one collision. Hence a complete clock sweep executes one of the
machine's elementary move-and-collision steps.

Decorate the 15-letter clock outward moves with G and return moves with
F. Its complete clock graph and conserved legality penalties are unchanged.
This gives a uniform positive NN propagation Hamiltonian of dimension

    d_shuttle = 15w = 11610D.

For D=q+1, 11610 is the concrete linear prefactor for this module only.
It is not the K of the completed machine in readout_machine_contract.md.
In particular it excludes input encoding, generation, verification and
any other work factors. The formula K=30 s0 r^2 applies to the different
two-lane architecture; substituting this module's 18 and 21 there would
not describe this machine.

All controls and marks are even. Give fresh and updated overlays the
grade of their underlying record. Identity and ordinary SWAP preserve
the total product grading on the two sites. The onsite collision does
not change records. The entire circuit therefore preserves grading on
all states, including malformed ones. No fermionic SWAP convention is
being assumed for the microscopic tensor-product sites.

## 4. Clean ports at every clock vertex

The first eight segments move LP to L, LC to LI, RP to R and RC to RI.
The paths do not disturb an already placed record. LP travels left from
a, leaving LC at a+1 unchanged. LC then stops at LI rather than entering
L. RP travels right from b+1, leaving RC at b unchanged; RC stops at RI.
The ordering of source positions keeps the two pairs separated.

No flag is enabled during any of these record moves. After phase 7,
the next five phases use identity quantum gates only. The right flag is
enabled first, the left second, the right cleared third and the left
cleared fourth. Only after both are clear does the inverse transport
start. It returns the whole record permutation to identity, not only the
four designated records. The neutral configuration also returns to its
start. The cycle length is

    T(n,a,b) = 10n+8a-4b-12.

Indeed the first eight distances sum to 4a+3n-2b-4, the readout itinerary
to 4n-4, and inverse transport has the same length as forward transport.
This is less than 14n for all the allowed positions, so one full cycle
fits among 2^n sweeps for all n>=8. There is no absorbing halt.

Define the endpoint readout to be READY times the projector onto no head
at that endpoint times the desired record operator. Define a colored
readout on an endpoint and its adjacent site with no-head projectors at
both sites. The head overlay is not silently identified with an ordinary
record by either physical observable.

At every active elementary clock vertex, the selected record is already
at its intended site. During a partial outward sweep the one head may
be fresh or updated, and during a partial return it may or may not have
collided. Those differences change only neutral control. If a head sits
on a readout support, the physical operator vanishes by definition. If it
does not, its pullback is exactly the intended source-record operator.
This proves the claim at all partial-sweep vertices, including flag
activation and removal, not just at completed sweeps.

The scheduling and all active-time sets are independent of quantum data.
Therefore if an input isometry V has source port maps A_L,A_R and source
colored maps D_L,D_R, its normalized clock history has maps

    A'_u = (N_port,u/C_n) A_u,
    D'_u = (N_density,u/C_n) D_u.

Here C_n=(2n+4)2^n+n-4. The counts are positive but need not agree with
each other or between ports. An incomplete last cycle adds only clean or
inactive vertices and does not spoil either identity. Input-independent
counts, rather than equal dwell lengths, are the necessary condition.

These statements preserve the correct C8 Clifford parity dressing and
its colored densities when applied to that prepared code. They also
apply to arbitrary prepared entangled codes. They do not construct V.

## 5. Finite checks and limitations

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 endpoint_shuttle.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 endpoint_shuttle_checks.py

The first writes the complete collision table. The second checks:

* all 378 collision entries and their inverse, preserving every static mark;
* 12096 fresh/updated bond cases, including inverse branches and grading;
* 75 complete cycles, n=8,...,32, totaling 14100 machine steps;
* every clock vertex at n=8,9,10, respectively 5124,11269,24582 vertices,
  including exact comparison with the machine after every complete sweep;
* source-code port and colored-density matrices at (q,m)=(4,2),(8,4),(8,8),
  with logical dimensions 4,16,256 and up to 4096 code entries.

All pass. Port and density matrix discrepancies are zero in the
implemented arithmetic. The clock active counts (left,right) are
(1028,1040), (2300,2300), (5276,5356) for ports, and (884,880),
(2020,2020), (4704,4628) for densities. They are not CFT estimators.

The source code in the algebraic test has an abstract residual-environment
label. The test does not assert that this label has already been compiled
into this module's D-dimensional records. The finite-history check proves
transport as an operator identity for arbitrary records, while the
separate source-code test verifies the inherited Clifford formulas. Their
composition is exact, conditional on the actual prepared input code.

Remaining tasks are an autonomous input encoder/verifier and a counted
composition with this module. For C19, the excited-state spectral measures
are another requirement. A theorem about ground-code port compression
does not identify the source resolvents of the enlarged history cell.
