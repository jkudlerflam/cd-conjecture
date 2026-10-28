# C22 progress

2026-10-28. New run. Read project STATE, route-A background, STATUS, C7 selector, C19 spectral/pole result, ir_control, C20 and the collaborator's exact/graded/sparse-route summaries. C8 PROGRESS read; C21 and C23 PROGRESS were absent. All writes stay here. No network or subagents.

[Proved, derivation to be recorded] A one-particle bus on a tree with unitary data-valued hopping can be gauged to a scalar graph Hamiltonian. Same-source return resolvents are scalar, regardless of the number of poles. This excludes simple outward-and-return unitary wires, not looped buses or general spatial codes.

[Candidate, not yet validated] A looped three-species bus, with vacuum plus neutral/positive/negative excitations, alongside n free data qubits. Neutral hopping carries the local data Z; charged hopping does not. Local neutral-to-charged conversion flips the data spin. In the one-particle sector a gauge transform converts the local flips into Jordan-Wigner fermions. Empty bus is an exact 2^n ground code if its chemical potential exceeds hopping and conversion norms. Endpoint excitation sources have zero first-order compression. This is computation-free and admits a 28-letter C7 tensor alphabet (4 marker bus letters plus 3*4*2 bulk letters), pending complete bond/sector accounting. No successful hopping/IR claim yet.

Next: calculate the full boundary resolvent, including self-energy terms and the left/right parity strings, before interpreting pole count or claiming n gapless bands. Exact small-sector numerical checks will distinguish this candidate from an arbitrary matrix gadget.

## Exact local architecture and obstruction

[Proved] Merge the bus with the selector token. There are 17 states: L/R carrying a data bit (4), four X-token species carrying a data bit (8), one empty marker, and four occupied marker species (5). A cell has n data sites and a terminating marker, one token with n+1 positions, and vacuum-token ground wavefunction sin(pi(j+1)/(n+2)). The code has exactly 2^n states. Endpoint source amplitude is a_n=sqrt(2/(n+2))*sin(pi/(n+2)); all active first-order compressions vanish. Exact charge is data occupation plus bus charge, and its mod-two reduction is product logical parity. This is a C7-type spatial walker, not a computation history.

[Proved] A neutral-token Z string gauge gives local fermion conversion in the excited sector. Spatial normal modes already supply n+1 distinct positive poles at conversion g=0. Both ports and the complete finite cell have been constructed in spatial_bus.py. The legal-cell dimension is 4(n+1)2^n, with every vacuum-walker excitation retained.

[Proved and numerical] For equal bus hopping magnitudes and equal charged conversion amplitudes, the order-lambda^2 g^2 one-cell terms are scalar. Flux-free hopping cancels. Uniform flux gives hopping coefficients proportional to
exp(i chi(n-i+j))*sin(kappa(n-i-j))*F_ij,
where F_ij>0 is an explicit spatial heat-kernel integral. This matrix is generally nonnormal. At n=2 its bottom-right entry vanishes and the two off-diagonal magnitudes differ. Its half-filled Bloch determinant is strictly negative, with a nonzero band gap. Thus this concrete 17-state design fails to give two independent gapless flavors. Uniform phases are also gauge-equivalent to a common overall phase times a real hopping matrix; they do not supply independently displaced flavor seas.

[Numerical] Exact g^2 resolvent expansion agrees with the independent heat-kernel formula below 1e-18 for n=1,2,3. Full legal-cell ED on two and three cells is under way, using exact total-charge sectors. First completed (n,M)=(1,2),(1,3),(2,2),(2,3) checks all include every legal cell state, not a one-excitation truncation. Final parameters and errors will be in RESULT.md. No central charge is measured. The construction does not prove a no-go for nonuniform flux mechanisms, more rails, noncommuting codes, or higher orders.
