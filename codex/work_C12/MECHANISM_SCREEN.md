# Screening simple mechanisms

These are exclusions of particular implementations, not a proof of the conjecture.

## Fixed MPS and frustration-free cores

For a fixed injective MPS tensor of bond dimension D, an open parent-Hamiltonian interval longer than the injectivity scale has at most D^2 boundary ground states. Increasing its length does not increase the logical dimension. A family whose bond dimension itself increases with the chosen cell length would require a new construction, not the usual fixed-tensor edge-state argument.

A ferromagnetic SU(q) cell of b fundamentals is more promising. Its frustration-free ground space is Sym^b(C^q), dimension binomial(b+q-1,q-1). If t_i^A is a constituent generator and T^A=sum_i t_i^A, permutation symmetry gives P t_i^A P = T^A/b. A weak single boundary exchange therefore projects to (lambda/b^2) sum_A T_left^A T_right^A. The large representation is real, but it does not amount to b gapless copies.

A deliberately tuned SU(q)_b WZW chain would have c=b(q^2-1)/(b+q). This is enough to make the idea worth retaining. The missing step is producing its special representation-dependent interaction. The first-order boundary projection gives only the bilinear, not the required higher-degree interaction. Higher-order gadgets must also suppress lower-order terms and furnish the relevant tunings. No explicit economical bond achieving this was found.

Colored-path or stack constraints can yield a ground space that grows with length. Their obvious boundary operators push or pop a stack. Counting stack words as independent CFT flavors is unjustified. An explicit dispersing, relativistic logical Hamiltonian remains necessary.

## Squeezed space and composite objects

For the supplied spin-charge construction, if c_i<=d_i-1 then c_1+c_2+1<=d_1+d_2-1. Repetition preserves the bound. It cannot be a counterexample by itself.

An extended molecule containing b q-state constituents has q^b internal states, but an ordinary collision exchanges or couples its exposed constituents first. A whole-molecule SU(q^b) permutation does not follow from its internal dimension. If microscopic exchange is forbidden entirely, the internal labels become frozen zero modes, not a collection of unitary CFTs. Protecting whole-molecule exchange while canceling partial exchange is the additional engineering problem.

The odd-conversion qutrit model in DERIVATION.md is a different squeezed-space possibility. It seeks extra momentum valleys of mobile bound objects rather than counting their internal labels. It has an explicit nearest-neighbor bond and passes the two-body test, but its many-body CFT is unresolved.

## Gauge constraints

A simple local Z2 gauge constraint does not compress k on-site fermions to k-1 fermions for free. Matter has dimension 2^k per site and a link gauge qubit contributes a factor two. Imposing one Gauss constraint per site removes that factor again, leaving the original extensive Hilbert-space count, up to global sectors. Eliminating the link qubits instead makes matter hopping nonlocal or constrains it. This count does not exclude more elaborate emergent gauge mechanisms; it excludes the naive local-parity compression argument.

Ground-state sector labels and nondispersing gauge modes must not be added as central charges. Symmetry can protect gapless modes but their existence still needs a Hamiltonian and a linear spectrum.

## Ladders and spontaneous periods

Putting k independent q-state critical chains on one rung gives d=q^k and c=k(q-1). Since q^k >= 1+k(q-1), this obeys c<=d-1. Unrolling a rung into a long spontaneously formed cell saves local dimension only if the cell can transmit k independent logical channels through a single physical cut. A static sequence of nearest-neighbor pairings is still a one-dimensional path; drawing it as a ladder does not create the missing leg couplings.

The supplied operator-rank observation is the useful distinction. Large logical dimension alone is insufficient, but low-rank fermion-like couplings could transmit many modes through a modest port. A simple frustration-free encoding with this boundary operator system was not found. The failure of the elementary constructions does not prove that a computation track is necessary.

## Nonrelativistic and sector-counting pitfalls

An SU(q) ferromagnet can have large logarithmic entanglement and a large finite-size ground multiplet. Its quadratic magnons fail the z=1 condition. The same warning applies to the vacuum thresholds of the qutrit candidate.

Translated copies of a spontaneously broken ground state are alternative vacua of the same critical theory, not independent tensor-factor CFTs whose central charges should be added. No counterexample based on multiplying the number of such sectors is claimed.
