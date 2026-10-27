# C10 progress

The task is closed book. I will not inspect existing files or use outside references. All work and generated artifacts stay in the current directory.

Initial plan: analyze the compatibility of port-idle propagation, initialization, and energetic Bell preparation; give a finite explicit short-program track where possible; locate any missing uniform compiler and illegal-sector estimates; verify a small example with exact diagonalization. No universal construction or fixed local dimension will be claimed without those ingredients.

## First structural findings

A hard-coded finite program admits a transparent nearest-neighbor track. Give each physical position its own background labels A_j,R_j and each clock time its own head label H_t, whose prescribed position is j_t. A legal cell is A_1...A_(j_t-1) H_t R_(j_t+1)...R_m. Neighboring clock times differ by one head step, with the gate carried by that bond. Marker boundaries and forbidden A-R interfaces enforce exactly one head. This proves a path clock and an illegal-sector penalty, but the alphabet size is 2m+T+2 and therefore is not a fixed-alphabet compiler.

There is an independent exactness problem. For a port-idle circuit U_t, the propagation kernel is V psi = (T+1)^(-1/2) sum_t |t> U_t psi. A time-local energy D = sum_t |t><t| h_t preserves this entire kernel only if U_t^dagger h_t U_t is the same operator for every t. On an individual history eigenvector psi, every such operator must have the same eigenvalue and eigenvector. An energy present at one selected time generally fails this condition. Consequently an arbitrary prescribed program with selected-time energies does not have its uniform history space as its exact ground space at finite clock coupling. This will be proved and tested rather than assumed away.

No existing files have been read or listed.

## Explicit small Route B program

For q=2, choose F_1=X, F_2=Y and gamma_1=X, gamma_2=Y on one logical qubit. Direct Pauli algebra gives

B+B^dagger = (X_p X_l + Y_p Y_l + X_s X_l - Y_s Y_l)/2.

A four-site cell has port p at site 1, partner s at 2, logical l at 3, and an idle right port at 4. Three clock states put the head at positions (2,3,2). The two gates are SWAP_(2,3), then identity. At time 0 put 3 Delta(1-P_Omega) on bond (1,2) and (3g/2)(XX-YY) on bond (2,3). At time 2 put (3g/2)(XX+YY) on bond (1,2). Their uniform-history compression is exactly H_L. Both gates leave both ports idle; the Bell projector appears only as energy. The alphabet has 12 letters, so this short example has d0=96.

This construction realizes the requested operator after history compression, not as an exact finite-J ground-space restriction. Exact diagonalization will quantify that distinction. The full prescribed-program version uses k_c=2m+T+2 and therefore fails the fixed-alphabet requirement as N grows.
