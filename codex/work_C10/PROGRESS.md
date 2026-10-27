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

## Numerical results

`python` was unavailable; `python3 verify_c10.py` succeeded in under one second. The script and `verification.json` were created in this directory.

The exhaustive classical check examined all 12^4 = 20,736 computation strings in a four-site cell, including misplaced M symbols. Exactly three are legal, and their transition graph is the intended three-vertex path. A transfer-matrix check finds no legal cell of any tested length 1 through 8 except length 4. The general proof covers all lengths.

The 48-dimensional legal-sector ED gives the target compressed Route B Hamiltonian to 3.8e-16 and the ideal port compression to 5.6e-17. Port operators preserve the unrestricted propagation-history space exactly. The purported uniform Route B ground histories have eigenvector residual 0.328264, independent of J, so they are not exact eigenvectors. At J=100 the true ground energy is -0.07658245, compared with the compressed value -0.07445626; the port-compression error is 0.00696736. At J=1000 and 10000 it falls to 0.000708510 and 0.0000709677.

Bell preparation alone, as a positive energy term at time zero, has the exact history ground space. For J=7 its measured gap is 0.635209, above the proved 0.269231 bound. The clock tower is (0,3.5,10.5).

An explicit illegal no-head block has energy kappa=20. A uniform per-head shift C=3|g|=0.6 makes all clock-time energies nonnegative. The analytic bound for any illegal cell, allowing shared endpoint synchronization penalties, is kappa/2=10 before subtracting the valid-cell energy. An exact local marker counterterm subtracts that energy and leaves the bound at least kappa/2-C=9.4, preserving the supplied selector ordering and its margin. This is proved only for the hard-coded finite track.

## Full bond and exhaustive cell-sector verification

`build_bond_c10.py` now builds the actual translation-invariant two-site bond in local order selector / computation / data. For d0=96 the sparse bond is 9216 by 9216 with 10,470 nonzero entries. Its Hermiticity error and same-bond selector commutator are zero. The overlap-bond commutation is proved algebraically from dependence only on the marker projector.

`verify_full_cell_c10.py` diagonalized all quantum blocks for all 20,736 four-site computation strings with correct endpoint markers, total dimension 331,776. There are 19,876 blocks, the largest of dimension 144. The illegal-sector minimum is exactly 20 before the marker counterterm, and 19.47658245 after it. Independently compressing the built physical bond onto the legal cell matches the expected matrix to 1.42e-14. Marker mismatches and all other lengths are covered by the analytic kappa/2 bound, not this finite enumeration.

At J=100, the actual compressed port spectrum is +/-0.355122674, whereas the ideal spectrum is +/-0.348155312. This discrepancy cannot be removed by changing the logical basis. The second-order ground-energy coefficient is -0.2118188985/J, consistent with the finite-J scan.

The final report will distinguish the explicit hard-coded construction and its proved bounds from the missing fixed-alphabet compiler and the false exact-history claim for arbitrary selected-time energies. No outside references or pre-existing files have been accessed.

## Improved selector-preservation argument

The final bond uses an explicit counterterm rho=C=3|g| per marker, rather than a computed ground energy. Let e0 be the valid-cell ground energy before this counterterm. The history variational state gives 0 <= e0 <= C+E_minus <= C. Thus a=(e0-rho)/b <= 0. A good cell contributes at least a b, while a bad cell of any marker spacing B contributes at least kappa/2-rho >= 0 >= a B. Summing over cells proves H_comp >= a L, attained in the desired all-good b-periodic sector. Adding the computation therefore cannot reduce the selector energy margin. This argument uses only rho >= e0 and kappa/2 >= rho, and works generally with rho=C+v, v=max_t ||h_t||.

The default built bond now has rho=0.6, bare marker chemical potential mu-rho=0.3, and the conservative illegal-cell bound 9.4. The optional exact-counterterm mode remains available. Earlier reported 19.47658245 values used the exact counterterm; the final default exhaustive test instead has illegal minimum 19.4.
