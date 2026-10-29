# C84 alphabet accounting

[Proved] The certified alphabet has eleven letters. The reduction from
C80 is E_1 deletion, with one unused boundary bit and one added forbidden
pair type L_1 followed by a marker. It retains both charged species and
the folded data/express graph. The member has eleven active flavors on
twelve data sites plus one marker.

| Function | Letters in C80 | Letters here | Status |
| :-- | :-- | :-- | :-- |
| Binary L/R backgrounds and vacuum head | 4 | 4 | Retained |
| Charged data and express carriers | 4 | 4 | Retained |
| Markers, including displaced-bit storage | 4 | 3 | E_1 removed |
| Total | 12 | 11 | Certified at n=11 |

## What is necessary and what is a design choice

[Proved, restricted] The code needs enough independent payload states,
the propagation needs independent endpoint filters for its active
flavors, the conversion needs a carrier charge change, and the selected
frozen cells need their boundary observable. These are necessary
functions in the retained mechanism. They do not each require a
separate letter. The vacuum head already shares R_0, and the empty and
vacuum markers already share E.

[Proved] The marker need not store a data bit. Fixing that bit costs one
flavor per cell, whereas the number of flavors still grows without bound
with the cell length. This is the overlooked saving in the dense shifted
encoding. Omitting E_1 without accounting for the lost flavor would be
wrong. The scalar graph must retain its extra inactive endpoint.

[Proved, restricted] In a faithful product encoding of the folded graph,
D_s and X_s at the same position and payload differ only on that one
site. Distinct masses and rail amplitudes require orthogonal letters.
The binary payload needs L_0,L_1 and R_0,R_1; the orientation distinguishes
the clock positions and enforces its single wall. A neutral port with
two charged outputs needs three marker states. These disjoint sets give
4+4+3=11. This bound has hypotheses about where the information is
stored. Reencoding the rail or orientation in neighboring sites evades
those hypotheses and requires another local construction.

[Proved] The full additive-charge solution has three free parameters u,
v,delta. Background charges are u,u+delta. Carrier charges are
u-delta/2,u+delta/2, and marker charges v,v-delta/2,v+delta/2.
The logical U(1) needs delta nonzero. Hence marker states cannot merge,
A and B cannot merge, and a carrier cannot merge with a background
without destroying that logical charge. Frozen on-site marker positions
exclude marker/data identifications. With no rail or orientation count,
these rules yield only a>=7. Keeping orientations gives a>=9; keeping
the faithful fold gives a>=11. PROOF Section 8 derives the equations.

[Proved, weak] Code dimension alone gives a>=3 for a dense unbounded
family with one frozen on-site marker per cell and b=n+O(1). One
nonmarker state and a fixed finite marker space cannot encode 2^n states.
There is no stronger universal lower bound established here.

## Rejected or unresolved reductions

[Proved, restricted] Renaming A,B as bit values of a single token does
not save a state: a charged binary token still has two orthogonal local
letters per rail. Sharing those states with L/R instead is excluded by
the charge equations above for the literal conversion and port. A new
port that also changes payload occupation is outside this argument.

[Proved, restricted; exact finite word counts] Identifying D_s=X_s
collapses 72 legal words to 60 at p=3 and identifies vertices with
different diagonal energies. This destroys the faithful folded graph.
It does not prove that every single-rail mechanism is impossible.

[Numerical] L_0=R_0 collapses those words to 64 and gives vacuum levels
-3,-2.61803,-2,-1.80194 at K=100. L_1=R_1 is injective on the old words,
but the enlarged zero-penalty grammar gives levels
-2.020008,-1.624405,-1.618034,-1.549903. The intended code level is
fourfold -1.618034. These are failures at order one. Integer Rayleigh witnesses in the
zero-penalty subspaces give energies strictly below -2 in both cases,
whereas the target -2cos(pi/5) is greater than -2. Larger grammar
penalties do not remove paths already in its zero-penalty space.

[Proved, restricted] The other six equal-charge single identifications
merge a marker with data. They no longer admit the required on-site
projector whose eigenvalues distinguish frozen marker positions. More
strongly, each quotient has a connected on-site transition graph, so
every diagonal on-site projector conserved by the rules is trivial.
All ten cases, their exact word counts, and the numerical spectra are
in `alphabet_checks.json`.

[Inherited rejected alternatives, not rerun] C75 and C80 already reject
the literal identical-real-carrier fold deletion at n=2. Its two spectral
intervals are disjoint, so a common field cannot cross both. C80's fast
B replacement also fails at n=2,3 for its five tested speeds, through
loss of onsite cancellation. Neither test excludes arbitrary redesigned
carriers or a controlled counterterm construction. These failures were
read before this search; their negative results were not recast as
universal lower bounds.

[Proved, restricted] A tree of scalar hopping amplitudes has no
gauge-invariant flux. If the entire effective quadratic Hamiltonian can
be made real by an on-cell gauge, its spectrum is even in k and the
sea centers coincide modulo pi when each crossed band has one connected
sea, as in the retained reference. A nontrivial loop
or another source of time-reversal breaking is needed for the retained
pair-sum strategy. A loop in configuration space does not by itself
require an express letter. Consequently this observation cannot justify
a universal two-rail alphabet count.

[Heuristic] Below eleven, the promising unresolved issue is how to
store rail or orientation information in a context without adding
zero-penalty branches or losing the adjacent port. No ten-letter bond
with those properties and a finite certificate was found here.
