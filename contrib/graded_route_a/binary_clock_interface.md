# Interface and limitations of the initialized binary clock

The existing C8 fifteen-letter full-sweep extension already remedies
the sparse schedule analyzed here. See the accompanying
[full-sweep interface](sweeping_clock_interface.md) for the stronger
available clock and its remaining interpreter obligations.

## Claim and status

[proved] The nine-letter clock of
`codex/work_C8/initialized_binary_clock.md` admits arbitrary commuting
diagonal nearest-neighbour input-format checks during its initialization
sweep, with the same complete-history-space argument and gap bound.
Its work-gate schedule has an exact nonuniform edge-use budget: edge
j|j+1 in an n-site cell is crossed 2^(n-j) times in each direction
after initialization.

This gives an interface for a possible compiler and an obstruction to
a naive repeated full-tape-sweep implementation. It does not establish
universality of this schedule or rule out a more careful compiler.
A reversible cellular automaton or Turing machine has not been embedded
into this clock here.

Author: OpenAI Codex, GPT-6 family, 2026-10-28. The existing clock source
was read; the extension and counting argument were independently derived.
`clock_schedule_checks.py` checks n=1,...,10. No human review is claimed.

## 1. Commuting formatting checks preserve the history theorem

Write C=C_n=n+4(2^n-1). During the first n configurations the I head
visits every site and the work circuit is identity. Let Q_a be mutually
commuting work projectors, each supported on one site, two neighbouring
sites, or an allowed endpoint. Assign each Q_a to an initialization
time t_a at which a repeated local term can check its support. For
example, a two-site check Q can be coupled to the clock pattern IZ on
the same two sites; that pattern occurs once per internal bond during
initialization. Boundary terms can check the endpoints.

Add positive terms eta_a times that clock projector tensor Q_a,
with eta_a>=eta_min>0. Diagonal checks can enforce adjacent tag patterns,
blank workspace and a single domain-wall head. They need not select a
tensor-product input space. Define A=intersection_a ker Q_a. The complete
legal zero space is exactly

    V psi = C^(-1/2) sum_t |t> W_t psi, psi in A.        (1)

Proof: gauge away the post-initialization circuit with controlled W_t.
The pins remain unchanged since W_t=I at every pin. In the joint
eigenspaces of the commuting Q_a, accepted work inputs see the path
Laplacian; every rejected input has a potential >=eta_min at some
vertex. The original telescoping bound therefore gives

    gap >= min{K,J[1-cos(pi/C)],eta_min/(2C),J/(4C^2)}.    (2)

Illegal clock sectors still cost K because the added terms are positive
and clock-defect projectors remain conserved. This proves completeness,
not just the existence of the intended histories.

The dimension of A must be counted. If the format admits arbitrary
programs or delimiter locations, all of them survive in (1). The clock
does not select a preferred program sector implicitly.

## 2. Exact edge-use budget

Number sites from the least significant clock digit 1 to digit n.
Carry crosses edge j|j+1 when the current counter value has its lowest
j bits equal to one. Among all 2^n values, including terminal overflow,
exactly 2^(n-j) satisfy this condition. Thus

    number of carry gates on edge j = 2^(n-j).           (3)

Every such carry except terminal overflow returns. The initial all-zero
return supplies the missing one, so

    number of return gates on edge j = 2^(n-j).          (4)

Initialization crosses once more with its assigned identity. The total
number of potentially nontrivial two-site gates on edge j is

    B_j = 2^(n-j+1).                                    (5)

Each direction splits equally between its b=0 and b=1 rules. A one-site
C0->R1 turn occurs at site i exactly 2^(n-i) times.

The rightmost edge is consequently used only four nontrivial times,
regardless of n. Exponential total clock duration does not provide
exponentially many full-tape sweeps. A conventional repeated sweeping
machine cannot simply be attached to carry/return transitions and
inherit its intended update schedule.

This bounds the available interactions, not general computational power.
A circuit could compute in a prefix and transport a small final output
to the right once. The sparse port isometry might admit such a design.

## 3. Padding helps but must itself be initialized

If an active prefix ends at site a>=2, its last internal edge is used
2^(n-a+2) times. With p=n-a padding sites, a direct implementation
needing T gates on that edge therefore requires

    T <= 2^(p+2), equivalently p >= ceil(log_2 T)-2.      (6)

This is necessary, not sufficient. The gates still appear in the nested
carry/return order. The computation/padding split must also be selected
without introducing extra ground sectors.

Purely diagonal local formatting has a precise limitation. Fixed
allowed adjacent symbols and endpoint symbols on a finite alphabet
define a regular language: the allowed words are paths in a finite
directed graph. Forgetting local tags preserves regularity. Such
constraints cannot enforce the nonregular format

    data^m separator padding^m

for unbounded m. This does not cover quantum length filters,
computational verification, or a target-specific finite length list.
It excludes only an attempted all-length equality check using fixed
diagonal nearest-neighbour zero-energy constraints.

A compiler admitting arbitrary splits must prove all of them harmless,
derive the split during computation, or computationally reject the
wrong ones. Selecting the intended sector by hand would leave the
complete-ground-space obligation unsatisfied.

## 4. A circuit identity the schedule does supply

Assign a forward gate U_j to both carry labels on edge j and U_j^dagger
to both return labels. An increment whose first zero lies at site i
has net work action

    U_1^dagger ... U_(i-1)^dagger G_i U_(i-1) ... U_1,    (7)

where G_i is the turn gate. In the repeated-bond model the U_j are
copies of one two-site gate, possibly controlled by local initialized
tags, and the G_i are copies of one one-site gate. If the turn is
inactive, the partial sweep cancels exactly.

This suggests reserving complete active-prefix sweeps for designated
events and cancelling shorter ones. It is not a universal interpreter:
local reversible tags must recognize those events, carry a controller,
initialize any delimiter, and handle the initial return and terminal
unreturned carry. No complete rule set doing this is claimed here.

## 5. Halting, idling and endpoint flags

The unique terminal clock word is L^(n-1) C1. Its right endpoint
identifies termination locally, but its left endpoint is L, which
occurs at many earlier times. The terminal clock alone therefore
does not isolate a simultaneous unique readout at both endpoints.
There is no later return sweep to broadcast terminal overflow leftward.
A valid output protocol must establish clean flags beforehand or use
a different locally checked readout arrangement.

A finite-dimensional unitary acting as identity on an entire READY
subspace cannot map an orthogonal non-READY state into that subspace.
This follows immediately from preservation of inner products. More
generally a finite-dimensional invariant READY subspace has no incoming
transition from its orthogonal complement. A repeated reversible update
therefore cannot implement an absorbing halt state. Reversible idling
can advance auxiliary counters while preserving outputs, but the
counter, initialization and lack of premature wraparound need proof.
Different successive gates may also prepare then preserve an output;
the statement is about one repeated update, not all halting protocols.

Left and right readouts need not have equal weights. Known nonzero
weights can be absorbed into couplings. The obligation is that a port
operator vanish at every contaminated time and realize the same desired
logical operator at every time when it is active.

## 6. The precise remaining computational lemma

A sufficient completion using this clock would supply fixed rule
unitaries and a fixed-overhead work alphabet such that:

1. Commuting input checks select exactly the intended logical inputs
   and initialized auxiliary registers, including format registers.
2. The actual binary-clock word computes the sparse encoding or clean
   readouts, with operation counts fitting the edge budgets.
3. Every active endpoint readout has its claimed compression and an
   exactly known nonzero total weight.
4. Every local update preserves the required parity or species charges.
5. Every length-dependent program is generated from blank inputs and
   length, with no unverified program or split sector.

Once these are supplied, the existing clock and Section 1 provide
legality, exact histories, initialization and spectral control. They
do not currently provide the computational lemma itself.

## Reproducing

Run `python3 clock_schedule_checks.py`. It uses the standard library
and writes `clock_schedule_results.json`. It takes under one second.
