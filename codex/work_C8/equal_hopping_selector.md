# An exactly soluble mobile selector on every ring

Status: the Hamiltonian mappings, sector bound, entropy identity, and finite-size
energies below are proved. The identification of a free Fermi sea with its
compact c=1 CFT is the stated standard continuum input. This supplies a
microscopic c=1 model at d=4, not a counterexample to c <= d-1.

## 1. Hamiltonian and legal arrows

Use the four letters #,L,X,R with arrows

    #=(1,0), L=(0,0), X=(0,1), R=(1,1).

Write Q for the projector onto pairs whose touching arrows disagree and
M=|#><#|. Set a>0 and

    h = K Q + (mu/2)(M tensor I + I tensor M)
        -a(|LX><XR| + |XR><LX|)
        -a(|#L><R#| + |R#><#L|).                         (1)

This is the preceding mobile-selector Hamiltonian at t=u=a, with no small-u
projection. It is one-site translation invariant and nearest neighbor.
Every Q_i commutes with H because each allowed move preserves its three
affected legality indicators. Thus the full Hilbert space splits by the
binary violation pattern epsilon_i=Q_i.

In the legal sector let b_i be the right arrow on site i. Its left arrow is
b_(i-1). Every binary sequence b gives one legal basis vector, including
on odd rings. The legal Hilbert space has dimension 2^L.

A hop at bond i flips b_i precisely when b_(i-1) != b_(i+1). Both orientations
have the same amplitude -a. Let Z_i|b_i> = (-1)^b_i |b_i>. The restricted
Hamiltonian is therefore

    H_legal = -(a/2) sum_i X_i(1-Z_(i-1) Z_(i+1))
              +(mu/4) sum_i (1-Z_(i-1) Z_i).             (2)

The chemical term follows because # counts the 1-to-0 transitions. On a
ring they are half of all transitions. This step would fail on an open
interval without a boundary term.

## 2. All forbidden sectors stay uniformly above the ground state

In a fixed violation sector, the left arrow of site i is

    a_i = b_(i-1) xor epsilon_(i-1).

Consequently every violation sector also has dimension 2^L and the same
binary-coordinate Hilbert space. The comparison is not a variational
identification of a few basis states.

For a violation on bond i, that bond's hopping is absent. This changes the
operator by norm at most a. The hopping coefficient on bond i+1 can also
change, because its left outside arrow a_(i+1) has flipped. The two
conditional flips differ by norm at most a. All other kinetic terms agree,
apart from the analogous effects of other violations. The marker projector
on site i+1 changes by norm at most one. It costs at most |mu| in operator
norm. The triangle inequality gives, for k violations,

    ||H_epsilon - K k I - H_legal|| <= k(2a+|mu|).        (3)

This remains valid when affected neighborhoods overlap. Thus

    E_min(epsilon) >= E_min(legal) + k(K-2a-|mu|).        (4)

Choose K>2a+|mu|. The legal ground state is then a ground state of the full
d=4 chain, and every illegal sector is separated from it by a positive,
L-independent gap at least K-2a-|mu|. This includes all ring lengths L>=3.
The small L=2 convention, in which the same unordered edge occurs twice,
is irrelevant to the limit and is not included.

## 3. Free domain walls, including both ring twists

Define wall occupations w_i=b_(i-1) xor b_i. Their total W is even on a ring.
Every even wall pattern has two preimages, related by the global flip
G=product_i X_i. The move in (2) exchanges w_i,w_(i+1)=(1,0) and (0,1).
It is nearest-neighbor hard-core hopping with amplitude -a, and the chemical
potential is mu W/2.

Fix the orientation b_0 when reconstructing b from w. A hop across that
reference changes b_0; every other wall hop leaves it fixed. In a sector
G=g, g=+1 or -1, that crossing has amplitude -ag. This is the hard-core
particle boundary twist. After Jordan-Wigner, even W gives fermion
boundary phase -g. The complete spectrum is the union over

    W even, g=+1,-1,
    k_j = 2pi(j+phi_g)/L,
    phi_+=1/2, phi_-=0,
    E = sum_(occupied k) [-2a cos(k) + mu/2].             (5)

Both twists are necessary. Counting all these states gives
2 times 2^(L-1)=2^L, the whole legal Hilbert space. This also proves there
are no omitted extensive zero-energy sectors.

For |mu|<4a the dispersion has two simple Fermi points

    cos(k_F)=mu/(4a),       v=2a sin(k_F)>0.               (6)

Linearization at them gives a massless complex fermion. The even-W
projection and the two twists impose finite global sector restrictions;
they change the compactification and allowed towers but not the local
stress tensor or its c=1 central charge. This last free-fermion continuum
identification is the standard input. At |mu|=4a the slope vanishes and the
band edge is quadratic, so those endpoints are excluded.

On odd or otherwise incommensurate rings, W and the filled momenta are
rounded within O(1) of their thermodynamic values. Their energies differ
at order 1/L. Each Fermi level has at most two orbitals with a given energy,
so ground-state multiplicities stay bounded, independently of L. No
immobile defect with L possible positions occurs. Irrational k_F gives
the usual size-dependent charge-sector offsets of a lattice Fermi sea;
there is no claim that its vacuum offset is identical for all L.

## 4. Energies, towers, and entropy at mu=0

Take L divisible by four. The positive ground state is in g=+1 with
W=L/2. Summing the occupied cosines gives

    E0(L) = -2a/sin(pi/L)
          = -(2a/pi)L - (pi a)/(3L) + O(L^-3).           (7)

With v=2a, the Casimir coefficient is -pi c v/(6L), giving c=1.
The lowest state in g=-1 has energy -2a cot(pi/L). Its gap is

    2a tan(pi/(2L)) = pi a/L + O(L^-3),                  (8)

so it has scaling dimension x=1/4. Particle-hole excitations near either
Fermi edge have energy 2pi v/L times an integer, with the ordinary
fermionic partitions producing the oscillator towers. Changing even W
by 2 has leading energy 4pi a/L, or x=1. Thus both the first intersector
gap and the oscillator scale are linear in 1/L.

There is also an exact entanglement identity. Let |psi_w> be an even-wall
state, and take its g=+1 legal-arrow lift. For a nonempty proper interval
A of physical sites, write n_A for its wall parity and s for its initial
arrow. The basis change within A identifies its allowed letter strings
with (w_A,s), and the complement with (w_B,s xor n_A). The lifted state is

    (1/sqrt(2)) sum_(w_A,w_B,s) psi(w_A,w_B)
       |w_A,s>_A |w_B,s xor n_A>_B.                     (9)

Tracing B forces equal w_B. Global even parity then forces equal n_A
on the bra and ket, and hence equal s. Therefore

    rho_A^physical is unitarily equivalent to
       (rho_A^walls/2) direct-sum (rho_A^walls/2),
    S_A^physical = S_A^walls + log(2).                  (10)

Zero eigenvalues from illegal interval strings do not affect the entropy.
This proves equality of the logarithmic entropy coefficients, rather
than estimating the effect of the encoding by an unspecified constant.

## 5. Independent checks and limitations

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 equal_hopping_selector.py`.
The output is `equal_hopping_selector_checks.json`.

The script compares every eigenvalue of (2) at L=5,6,7,8,9 and mu=0.37 with
the independent occupation sum (5), with maximum error 1.96e-14. It tests
(3) for all 15,31,63 nonempty violation patterns at L=4,5,6 respectively. It constructs the many-body Slater wavefunction and
computes physical-site reduced density matrices, without storing a 4^L
state vector. The XX wall state is the known c=1 calibration for entropy.

At L=20, fitting the interval entropy to log[(L/pi)sin(pi ell/L)] for
ell=5,...,10 gives c=1.003981. The half-chain entropy slope from L=16 to 20
gives c=0.992700. Each equals its XX control to numerical precision, as
(10) predicts. The Casimir fit at L=40,48,64,80,96,128,160, with independently
known v=2 and a 1/L^4 correction in E0/L, gives c=0.999999704.

This resolves the finite-u CFT question only at the equal-hopping point.
It does not prove that the long-cell, small-u regime is in the same phase,
and it contains no high-dimensional logical code. Moving from this
soluble selector to a large-c payload still requires a microscopic
Hamiltonian and a separate stability analysis.
