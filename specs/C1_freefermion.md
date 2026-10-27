# C1: Central-charge bound for nearest-neighbour quasi-free fermion chains

Closed book: do not read or list any files. Work from this statement alone.

## Setup
A 1D lattice of sites x in Z_L (periodic, L -> infinity). Each site carries n Majorana operators
gamma_{x,a}, a = 1..n, with {gamma_{x,a}, gamma_{y,b}} = 2 delta_{xy} delta_{ab}, gamma^dagger = gamma.
For n even the local Fock dimension is d = 2^{n/2}; n odd is allowed (fermionic chain with "dimension"
2^{n/2} per site, e.g. one Majorana per site).
Quadratic Hamiltonian H = (i/4) sum_{x,y} sum_{a,b} A_{ab}(x,y) gamma_{x,a} gamma_{y,b}, where A is real and
antisymmetric as a matrix indexed by (x,a),(y,b), and nearest-neighbour: A(x,y) = 0 unless |x - y| <= 1.

Case (a), translation invariant: A(x,y) = A(y - x), so A(0) = -A(0)^T and A(-1) = -A(1)^T.
Case (b), period p: A(x+p, y+p) = A(x,y), still nearest-neighbour on the original lattice.

Criticality assumption (CFT): the single-particle (BdG) spectrum has finitely many zero-energy points; at each,
every band that vanishes there does so linearly (nonzero group velocity); there is no band identically zero.
You may use the standard fact that the scaling limit is then a theory of free chiral Majorana fermions, with
each chiral Majorana contributing 1/2 to c_L or to c_R. Define c := (c_L + c_R)/2.

## Tasks
1. Define the Bloch Hamiltonian h(k), state its particle-hole symmetry, and prove carefully that
   c_L + c_R = (1/2) x (number of zeros of det h(k) for k in the Brillouin zone, counted with multiplicity),
   under the linear-crossing assumption. Be careful with the double counting from particle-hole symmetry and with
   k = 0, pi. Check the formula on: critical transverse-field Ising (n = 2, c = 1/2), XX chain at generic filling
   (n = 2, c = 1), the one-Majorana-per-site chain H = i t sum_x gamma_x gamma_{x+1} (n = 1, c = 1/2).
2. Prove that in case (a) the number of zeros is at most 2 rank(A(1)) <= 2n, hence c <= rank(A(1))/2 <= n/2
   (= log_2 d for n even). Show c_L = c_R.
3. Case (b): show the same bound c <= n/2 holds for every period p, using that only the last site of a unit cell
   couples to the first site of the next, so the inter-cell coupling has rank <= n. Make the degree count of
   det of the (np x np) Bloch matrix, as a Laurent polynomial in z = e^{ik}, precise.
4. Is the bound attained for every n? Give explicit examples (n even and n odd).
5. Say what changes for (i) charge-conserving fermions with n_c complex orbitals per site (d = 2^{n_c}),
   (ii) non-Hermitian quadratic models with c replaced by c_eff (only a short discussion), (iii) Dirac points where
   two bands cross zero at the same k with linear slopes versus quadratic band touchings (not conformal).
6. Compare log_2 d with d - 1 for d = 2^m and state where equality holds.

Deliverable: theorem statements with complete proofs, then the checks of item 1 done explicitly. Be explicit at
every step and flag anything you assume rather than prove.
