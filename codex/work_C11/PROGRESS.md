# C11 progress

Closed-book calculation begun 2026-10-27. Only files created in this session will be read. The supplied context is treated as input, with any additional modeling assumptions recorded explicitly.

The first issue is whether the one-port Hamiltonian plus the stated mirror compression fixes a two-port cell Hamiltonian and its excited-state resolvent. First-order port compressions do not in general determine the second-order interaction. I will separate conclusions requiring only the given data from explicit one-port calculations and any proposed two-port completion.

The exact one-port resolvent has two excited gaps, D = sqrt(Delta^2 + 4 r g^2) and d = (Delta + D)/2. Its spectral port tensors were derived analytically and checked by dense diagonalization for q=2 and q=3. For q=2 the isolated-bond four-Majorana coefficient is positive at every nonzero finite g/Delta and starts at (83/2) lambda^2 g^4 / Delta^5. This statement is conditional on each endpoint having the supplied one-port resolvent. It is not yet the full two-port-chain answer.

An explicit additional model, Hcell = Hleft + Hright on two independent port-partner pairs sharing the logical register, has been diagonalized. At g/Delta=0.1 its ground degeneracy and Clifford port compressions hold numerically for r=2,4, but its compression coefficient differs from the supplied single-port alpha. Its second-order bond coefficients also differ. Thus importing the single-port resolvent into this completion is unjustified. The q=3, r=4 bond correction is anisotropic, with nonzero mixing between distinct flavor pairs. Orthogonality of port labels does not establish O(r) symmetry of higher-order processes.
