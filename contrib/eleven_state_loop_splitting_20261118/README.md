# Eleven-state loop: linear centre splitting needs finite express delay

2026-11-18. Work in progress; no microscopic CFT claim.
Author and AI assistance: OpenAI Codex, GPT-6 family.
No independent human or separate-agent review is claimed.

The same eleven physical states, nearest-neighbor rules and physical
C84 phases are retained. Only the already allowed closing amplitude
t_d is turned on. The inactive D_n and X_n endpoints remain in the
2n+3-vertex carrier. The previous contribution's biased clock and L1
readout are unchanged.

The exact derivative of the otherwise zero opposite-root pair sums is
implemented in loop.py. Hopping is even in t_d, while the imaginary
onsite term is odd. The first derivative vanishes if the tree hopping
is replaced by its strict symmetric Gram limit. At finite express
energy D the tree's antisymmetric delay term combines with the loop
onsite term to give a nonzero linear response in explicit families.
For eta=t_d D^(n+1), its leading size is eta/D.

The all-n question reduces to the explicit rational moment coefficients
in moments.py. The last coefficient is exactly -1/(16n-12), giving two
distinct leading centre responses for every n>=2. Full distinctness
for all n has not yet been proved. Small exact cases are retained as
finite results, not extrapolated to an unbounded-c family.

The literal modification adds no states and no new rule type beyond
restoring the two existing loop transitions. Finite-g, finite-K, port,
ground-sector and interacting infrared qualifications will be stated
against the new centre-splitting scale in the completed notes.

The source search belongs to the coordinating chat. The GitHub clue
has not identified a ten-state arbitrary-c construction; this folder
does not repeat that search. All edits are confined to this new folder.
