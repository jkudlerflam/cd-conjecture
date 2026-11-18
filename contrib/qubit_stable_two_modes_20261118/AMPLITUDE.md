# Full interacting quartet matching

## 1. Physical and effective conventions

Let b=(X+iY)/2 annihilate a down spin and n=(1-Z)/2. A uniform spin
rotation puts the four-term displayed Hamiltonian in the convention

\[
H=V\sum n_jn_{j+1}-\mu\sum n_j
 +\sum[w b^\dagger_{j+1}b_j+\mathrm{h.c.}]
 +\Delta\sum[b^\dagger_jb^\dagger_{j+1}+\mathrm{h.c.}],
\]

up to a constant, with Delta=a+b, w=i(a-b), V=4u, mu=4u-2h.
We take Delta,V>0. For the controlled cusp set

\[
\epsilon=\Delta/V,\quad t=\Delta^2/(2V),\quad
w=i\kappa t,\quad \bar\mu=(\mu-3\Delta^2/V)/t.
\]

Both kappa and mu-bar are fixed as epsilon tends to zero. This is an
expansion parameter limit, not an L-dependent coefficient choice.
With P imposing n_j n_(j+1)=0, the leading dimensionless number-conserving
operator is

\[
K=P\sum_j[i\kappa b^\dagger_{j+1}b_j+\mathrm{h.c.}
 -(b^\dagger_{j+2}b_j+\mathrm{h.c.})
 -2n_jn_{j+2}-{2\over3}n_jn_{j+3}]P.
\]

The effective Hamiltonian divided by t is

\[
K-\bar\mu N+\gamma(B^\dagger+B)+O(\epsilon^2),\qquad
B^\dagger=P\sum_j b^\dagger_jb^\dagger_{j+3}P,\quad
\gamma={4\over3}\epsilon.
\]

The omitted order-epsilon-squared terms have no direct Delta N=4 matrix
element on the full blockaded space: the path obstruction in C179 still
applies. They do shift the diagonal energies and the physical crossing
field. This distinction is essential in comparing an avoided crossing.
The path amplitudes and the full leading interaction are from the
[amplitude audit](../v4_locking_operators_20261118/AMPLITUDES.md).

All matching states below are eigenstates of K with the attractions
included. Their charge is the dressed blockaded N. Raw physical spin
density is not substituted for that charge. The microscopic spectral
comparison automatically includes high-energy dressing, rather than
assuming that a physical spin-pair operator equals its projected version
to all orders.

## 2. A definite, on-shell finite amplitude

Fix an even periodic ring, physical momentum K_tot=0 or pi and a common
Sigma eigenvalue. Let |i> and |f> be specified normalized eigenstates of
K in sectors N and N+4, with energies e_i and e_f. Choose

\[
\bar\mu_*=(e_f-e_i)/4,\qquad \bar e=(e_i+e_f)/2.
\]

At this field their grand energies coincide. If no intermediate N+2
eigenvalue equals e-bar, degenerate perturbation theory gives

\[
\boxed{ A_{fi}=-\langle f|B^\dagger
       (K_{N+2}-\bar e)^{-1}B^\dagger|i\rangle. }\tag{1}
\]

The physical off-diagonal coupling is

\[
t\gamma^2 A_{fi}={8\over9}A_{fi}{\Delta^4\over V^3}.
\]

Equation (1) includes **every** interacting intermediate state in the
accessible momentum/Sigma sector. Excluded symmetry sectors have exactly
zero overlap with either source. Off shell, the Hermitian canonical
second-order coefficient instead uses the average of the two denominators
at the initial and final grand energies; (1) is where they coincide.

`quartet.py` fixes each external eigenvector's largest component to be
real and positive. This is reproducible but not a physically preferred
continuum field origin. We quote |A|, not a supposed invariant sign of a
cosine. A finite-state phase convention cannot determine the continuum
vertex normalization or its scaling dimension.

Representative results, all with the full generated attractions:

| L | N to N+4 | kappa | mu-bar at matching | min intermediate denominator | magnitude A |
|---:|---:|---:|---:|---:|---:|
| 14 | 0 to 4 | 2 | -3.1037088582 | positive | 0.06381310621 |
| 16 | 0 to 4 | 1 | -2.1236411285 | 0.5573150584 | 0.02948434173 |
| 20 | 0 to 4 | 1 | -2.2213776734 | 0.6403859237 | 0.17246173011 |
| 24 | 0 to 4 | 1 | -2.2492847376 | 0.6356218514 | 0.23611731886 |
| 20 | 2 to 6 | 1 | -2.1374431332 | 0.7703198250 | 0.32887616260 |
| 24 | 2 to 6 | 1 | -2.1926273631 | 0.6764184377 | 0.51834521853 |

The N=2 to 6 rows have K_tot=pi and Sigma=-1; N=0 to 4 has K_tot=0,
Sigma=+1. Initial and final states are the lowest eigenstates *in those
declared sectors*. They need not be the full grand-canonical ground
states. In fact a lower-momentum-sector bound pair can pre-empt a vacuum
crossing. No thermodynamic density extrapolation is made from this table.

At L=16 the absolute sum of the spectral contributions is 0.8649224 while
their final magnitude is only 0.0294843. Retaining only a free denominator,
one intermediate state, or the absolute values before summing would miss
a large cancellation. The N=2 background rows are also a direct check that
the calculation is not limited to propagation in an empty lattice.

## 3. Exact vacuum resolvent and a connected-source identity

For N=0 the source has total momentum zero. A relative-coordinate
wavefunction has r>=2. At K_tot=0 the two nearest-neighbour hopping
contributions add to w+w*=0. Thus the imaginary hopping cancels **exactly**
in the full two-particle intermediate Hamiltonian. Even and odd r decouple.
B^dagger|0> occupies only the odd-r sector, whose boundary interaction is
V3/t=-2/3; the V2/t=-2 interaction remains in the inaccessible even sector.
It has not been discarded as a weak correction.

On the infinite half-line, put r=2m+3, m>=0. For mu-bar<-2, define

\[
 z={-\bar\mu-\sqrt{\bar\mu^2-4}\over2}\in(0,1),\qquad
 d_m={z^{m+1}\over2(1-z/3)}.
\]

Then

\[
[-2(\mathsf S+\mathsf S^\dagger)-2\bar\mu
  -(2/3)|0\rangle\langle0|]d=|0\rangle.\tag{2}
\]

This is the exact interacting vacuum-pair resolvent, not a band
approximation. All d_m are positive. Positivity **does not** fix the
overlap with an interacting four-particle eigenstate. The explicit
finite-ring inverse used in (1) includes periodic images exactly.

For a connected interpretation, extend the vacuum solution to a projected
pair creator D^dagger, so that

\[
D^\dagger|0\rangle=(K_2-2\bar\mu_*)^{-1}B^\dagger|0\rangle.
\]

Its coefficients are solved in the full finite-ring relative basis;
the antipodal pair is counted correctly. Define

\[
|R\rangle=(K_4-4\bar\mu_*) (D^\dagger)^2|0\rangle
           -2B^\dagger D^\dagger|0\rangle.
\]

Since (K4-4mu-bar*)|f>=0 at the vacuum crossing,

\[
\boxed{A_{f0}={1\over2}\langle f|R\rangle.}\tag{3}
\]

For a quadratic canonical-particle Hamiltonian the Leibniz rule gives
R=0. With the blockade and interactions it contains the collision
corrections to two independently dressed pairs. For mu-bar<-2 the
geometric decay in (2) makes this a localized connected construction;
the localization statement is not used at the threshold mu-bar=-2.
Equations (1) and (3) agree to <6e-15 in the saved vacuum cases.
`checks.py` separately shows exact cancellation of the on-shell spectral
sum for a generic quadratic four-fermion control. This guards against
mistaking a disconnected B^dagger B^dagger term for a quartet interaction.

## 4. Literal microscopic validation of a nonzero leading coefficient

`microscopic.py` builds the full occupation Hamiltonian of section 1,
with no blockade projection, in the even-parity, K_tot=0, Sigma=+1 sector.
At L=14,kappa=2 this is a 330-dimensional matrix. For each epsilon it finds
the avoided crossing by the Hellmann--Feynman condition that the two
levels have equal <N>. The crossing field includes all finite-epsilon
diagonal shifts. The next level is separated by more than 2.1t, so the
two-level interpretation is resolved in this symmetry sector.

The prediction from (1) is

\[
\frac{\mathrm{gap}}{2t\epsilon^2}\longrightarrow
\frac{16}{9}|A|=0.1134455221430.
\]

| epsilon | literal qubit value | ratio to prediction |
|---:|---:|---:|
| 1/10 | 0.10629337618 | 0.93695524 |
| 1/16 | 0.11057778205 | 0.97472143 |
| 1/24 | 0.11215863750 | 0.98865636 |

An additional calculation in the leading projected model approaches the
same prediction more rapidly, as recorded in `microscopic.json`.
This is a microscopic **finite-level** confirmation of the power counting
that was left as an unknown matching coefficient in C179. It does not
assign the same coefficient to a thermodynamic cosine or extrapolate the
controlled expansion to epsilon=1/4 at V4.

## 5. Exact cancellation control and its next surviving order

On L=12, define the real, normalized state

\[
|D_4\rangle={1\over\sqrt{30}}
 \sum_{\substack{|S|=4\\S\subset\mathrm{even}\ \mathrm{or}\ S\subset\mathrm{odd}}}|S\rangle.
\]

The t2=-1, V2=-2 terms on one sublattice are the uniform ferromagnetic
spin-exchange combination, giving K|D4>=-8|D4>. The imaginary nearest hop
annihilates it by cancellation between the two ways of reaching each
allowed opposite-sublattice configuration. V3 vanishes on it. Consequently
this eigenstate exists for every kappa; it is the lowest target state in
the declared sector at kappa=1.

Furthermore **B|D4>=B^dagger|D4>=0**. Annihilation would require a particle
on each sublattice. Creation would require three vacancies on the occupied
sublattice: two around the newly occupied opposite-sublattice site and
one at the other endpoint of the separation-three pair. There are only
two vacancies. This proves a zero amplitude in the leading projected
model to every order in gamma, not just a numerical second-order zero.

It is not an exact dark state of the microscopic nearest-neighbour chain.
For V=1, Delta=epsilon, w=i epsilon^2/2, write the microscopic chemical
potential as

\[
\mu=(3+\bar\mu/2)\epsilon^2,\qquad
\bar\mu=-2+{17\over5}\epsilon^2+O(\epsilon^4).
\]

After this diagonal matching, the common energy begins
E=-12epsilon^2-12epsilon^4+O(epsilon^6), in the occupation Hamiltonian
with vacuum diagonal energy zero. The first off-diagonal coefficient in
the (vacuum,D4) dressed basis is exactly

\[
\boxed{ M_{D4,0}= {6\sqrt{30}\over5}\epsilon^6
       +O(\epsilon^8). }\tag{4}
\]

Thus the minimum gap is (12sqrt(30)/5)epsilon^6+O(epsilon^8), after also
aligning the order-epsilon-six diagonals. Restoring V gives the coefficient
in the README. Numerical gap/(2t epsilon^4) values are 11.82306, 12.62326,
12.93426, 13.05033 at epsilon=1/12,1/20,1/32,1/48, approaching the exact
12sqrt(30)/5=13.14534138. The truncated projected model instead has a true
crossing, resolved numerically to rounding error.

### Exact certificate for (4)

`exact_dark.py` constructs rational complex matrices in unnormalized
translation/Sigma orbit bases; no decimal reconstruction of matrix
entries or guessed rational coefficient is used. The full symmetry
sector has 106 states: 11 blockaded and 95 with positive interaction
energy. The orbit norm of D4 is exactly 30.

Write H=D+epsilon A+epsilon^2 B+epsilon^4 C, where D counts adjacent
occupied bonds, A is adjacent pair creation plus annihilation,
B=(i hop+h.c.)/2-2N, and C=-mu2 N/2. First eliminate the 95-state sector
using its full D inverse, with E=e0 epsilon^2+e2 epsilon^4+... . If
R(epsilon)=(E-H_QQ)^(-1)=sum R_j epsilon^j, its recurrence is

\[
R_0=-D_Q^{-1},\quad
R_j=-D_Q^{-1}\sum_{k=1,2,4}X_kR_{j-k},\quad
X_1=A_Q,\ X_2=B_Q-e_0,\ X_4=C_Q-e_2.
\]

Terms with negative indices are zero. This gives the 11-state,
energy-dependent coefficients F2,...,F6 exactly. Its leading two-state
kernel has e0=-12. Let T span (vacuum,D4), and let R_l be (e0-F2)^(-1)
on the other nine blockaded states. The order-epsilon-four matrix is

\[
M_2=T^\dagger(F_4+F_3R_lF_3)T.
\]

With C initially zero, mu2=(M2_ff-M2_ii)/2=17/5 and e2=M2_ii=-12.
Recompute F4,F5,F6 with those counterterms. Put X=F3,Y=F4,Z=F5,W=F6.
The order-epsilon-six matrix, before its diagonal chemical counterterm, is

\[
T^\dagger[W+ZR_lX+XR_lZ+YR_lY
+YR_lXR_lX+XR_lYR_lX+XR_lXR_lY
+XR_lXR_lXR_lX-e_2XR_l^2X]T.
\]

The unknown next chemical counterterm is diagonal in T and cannot change
its first nonzero off-diagonal entry. In the unnormalized orbit basis
that entry is 6/5; multiplying by sqrt(30) gives (4). The script implements
the correct orbit metric throughout. `microscopic_dark_series` in
`microscopic.py` independently performs the same elimination with
orthonormal floating-point matrices. Literal finite-epsilon spectra are
a third check. This certificate is finite algebra, not a proof of a
many-body continuum amplitude at arbitrary L.
