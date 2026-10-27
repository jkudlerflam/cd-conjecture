# D2 numerical search, conventions and assumptions

This is a finite search of an eight-dimensional continuous family. It cannot
establish a universal bound. No external files or literature were consulted.
All computations use periodic boundary conditions and Pauli matrices, not
spin operators. The uniform bond, including its field term, is fixed as L changes.

## Parameters

Write the real coupling matrix as J = S + A, with S symmetric and A
antisymmetric. A proper rotation diagonalizes S. An orthogonal diagonalizing
matrix with determinant -1 can be made proper by changing the sign of one
eigenvector. In that basis, there are three diagonal exchange entries, three
antisymmetric entries, and three field components. Removing the overall
positive scale leaves eight continuous parameters. Degenerate S gives redundant
coordinates, but does not remove any Hamiltonian from this chart.

The scan normalizes sum_ab J_ab^2 + sum_a h_a^2 = 1. This is a convention,
not an operator-norm constraint. Random sampling in this chart is not a uniform
measure on equivalence classes. Both signs of exchange are sampled. The
antisymmetric convention is

```
A = [[0, Dz, -Dy], [-Dz, 0, Dx], [Dy, -Dx, 0]].
```

The implemented density is

\[
h_j^{(s)}=\sum_{ab}J_{ab}\sigma_j^a\sigma_{j+1}^b
+\sum_a h_a[s\sigma_j^a+(1-s)\sigma_{j+1}^a].
\]

The default is s=1/2. The density in the question is s=1. Their sums are
identical on the periodic chain. Their nonzero Fourier modes are not identical.
This distinction matters for the proposed matrix-element estimator.

## Momentum resolution

Let T shift bit j to j+1, and let r be a representative of a translation orbit
of period p_r. At momentum k=2 pi m/L, an orbit is allowed if m p_r = 0 mod L.
Its normalized basis vector is

\[
|r;k\rangle=p_r^{-1/2}\sum_{t=0}^{p_r-1}e^{-ikt}T^t|r\rangle.
\]

Applying T to this sum gives eigenvalue e^{ik}. If a computational-basis
transition in h_j takes r to s'=T^a r', its contribution to the momentum
matrix is

\[
\langle r';k|H|r;k\rangle\supset
h_j(s',r)\sqrt{p_r/p_{r'}}e^{ika}.
\]

This follows by expanding the source orbit sum, translating each output back
to r', and using the p_{r'}-term normalized destination sum. Repeated orbit
terms give the ratio of periods. The Fourier operator

\[
h_n=\sum_j e^{2\pi i n j/L}h_j
\]

satisfies T h_n T^{-1}=e^{-2 pi i n/L}h_n. Thus it maps m to m-n. The same
matrix-element construction uses the destination momentum in e^{ika} and
includes the additional phase e^{2 pi i n j/L} on each local transition.
In particular h_{-2} maps m to m+2.

The code checks the complete L=6 spectrum against a separate full-space
construction, and compares H|psi> and h_{-2}|psi> in full and orbit bases for
every momentum. Errors are about 1e-14. An early test caught a float16 square
root caused by storing orbit lengths as int8. The normalization now casts to
float64 before taking the square root. The recorded checks use the corrected
code.

Generic candidate refinements search every momentum for the ground state at
each L. The calibration and zero-field XXZ controls use their expected ground
momentum, 0 for L divisible by 4 and pi for L=2 mod 4, with Ising at 0. This
momentum convention is checked against all sectors at small L; its extension
to larger L is an explicit benchmark assumption. Fixed half filling is used
only for U(1)-symmetric, zero-field reference chains. A charge restriction is
not a valid option for a generic Hamiltonian.

## Entropy

Reshape the normalized ground-state vector into a matrix A for a block of
ell contiguous sites. If s_i are its singular values, the eigenvalues of the
reduced state are p_i=s_i^2 and

\[
S(\ell)=-\sum_i p_i\log p_i.
\]

For a periodic CFT vacuum the assumed scaling form is

\[
S(\ell)=\frac{c}{3}\log\!\left[\frac{L}{\pi}
\sin\!\left(\frac{\pi\ell}{L}\right)\right]+s_0+\cdots.
\]

The raw search statistic c_eff is the linear-fit slope using ell=2,...,L/2,
with the abscissa divided by 3. We also inspect ell>=3 and ell>=floor(L/4),
and the scaling of S(L/2) between sizes. The regression standard error measures
only fit residuals. It does not bound finite-size bias or establish a CFT.
Cut-window variation, size drift, and inconsistent spectral estimates are
reported separately. Entropy within a degenerate ground space depends on the
chosen state. The small-size solver resolves translation in its returned
ground multiplet, and records the translation defect and returned degeneracy.
It need not return an entire extensively degenerate multiplet.

## Energy and velocity

Assuming a single conformal velocity,

\[
E_0(L)/L=e_\infty-\pi c v/(6L^2)+b/L^4+\cdots.
\]

A fit of E_0/L against 1, L^{-2}, L^{-4} gives coefficient a=-pi c v/6.
Hence c=-6a/(pi v). The code records the lowest energies at relative momenta
plus and minus 2 pi/L, with candidate velocities

\[
v_{\pm1}(L)=\frac{L}{2\pi}[E_{\min}(k_0\pm2\pi/L)-E_0].
\]

Identifying the lowest such state as dimension one is an assumption and fails
already for Ising. A stress state of dimension two instead gives

\[
v_{T,\pm}(L)=\frac{L}{4\pi}(E_{T,\pm}-E_0).
\]

For unequal left and right velocities of one CFT, the Casimir expression uses
(v_L+v_R)/2. For several CFTs, the energy coefficient measures sum_a c_a v_a,
not sum_a c_a. Dividing by one observed velocity does not give c_T in general.
The report therefore calls energy estimates diagnostic unless the branch and
scaling assumptions are supported. Incommensurate finite-size changes can also
produce oscillatory Casimir terms, which the polynomial fit does not describe.

## Fourier-density matrix elements

If h_{-2}|0> has the continuum stress component
(2 pi v/L)L_{-2}|0>, Virasoro normalization gives
<0|L_2 L_{-2}|0>=c/2 and E_T-E_0=4 pi v/L. Substitution yields

\[
8\frac{|\langle T|h_{-2}|0\rangle|^2}{(E_T-E_0)^2}
=8\frac{(2\pi v/L)^2(c/2)}{(4\pi v/L)^2}=c.
\]

The needed assumptions are (i) the chosen state is a stress state, (ii) the
density has the stated stress projection without a dominating improvement
term, and (iii) a finite-size extrapolation is valid. Largest overlap alone
does not prove (i) or (ii). Both chiral momenta are computed.

The two chiral estimates are not added to define c_T. For a single CFT each
estimates the same central charge. For several sectors one would need the
sum over identified sector stress contributions within one chirality.

Degenerate eigenvectors are arbitrary. Within a numerically degenerate energy
cluster we orthonormalize the vectors and sum their squared overlaps. The
normalized projection of h_n|0> into that cluster is the largest-overlap
combination. Starting the target eigensolver with h_n|0> avoids missing the
coupled combination when the eigenvalue has multiplicity greater than the
number of degeneracy partners returned by an arbitrary Krylov start. Every
returned cluster, its weight, and the total captured weight are saved. The
default searches the lowest 10 levels in calibration and 12 in refinements;
states outside this window are not tested. Captured weight is a diagnostic,
not proof of the stress identification.

Several independent CFTs could have distinct stress energies. Selecting one
largest-overlap state would then miss other contributions to c_T. The saved
cluster data allow inspection of these branches. Unconditionally summing all
matrix-element ratios would include non-stress excitations and is not justified.

## Located failure of the literal largest-overlap prescription

For the critical Ising density -X_j X_{j+1}-Z_j, the largest-overlap state at
relative momentum 4 pi/L is not the stress state in these calculations.
At L=20 its inferred c_MV is 0.87046, whereas the lower stress-like state gives
0.49383. Symmetrizing the field selects the latter branch and gives 0.48174,
which extrapolates toward 1/2. This is a change in the local density, not in H.

Algebraically,

\[
h_n^{(1)}-h_n^{(1/2)}
=\frac{1-e^{-2\pi i n/L}}2\sum_j e^{2\pi i n j/L}
\,\mathbf h\cdot\boldsymbol\sigma_j.
\]

Although its prefactor is O(1/L) at fixed n, a local scaling field of dimension
x has a Fourier matrix element scaling as L^{1-x}. The improvement then
scales as L^{-x}. At x=1 it competes with the O(1/L) stress amplitude.
Thus its relative contribution need not vanish. The observed Ising failure
does not contradict its known c=1/2. It invalidates the unconditional use of
the literal largest-overlap prescription for an arbitrary density.

## A finite-density control for the energy fit

For J=diag(1,1,0) and field (0,0,h_z), Jordan-Wigner fermions have hopping
amplitude 2 and one-particle energy 4 cos k - 2 h_z, with the constant h_z L.
For N down spins, the spin-periodic boundary condition gives periodic fermions
for odd N and antiperiodic fermions for even N. Thus

\[
E_0(L)=\min_{0\leq N\leq L}\left[h_zL+
\sum_{\text{lowest }N}\bigl(4\cos k-2h_z\bigr)\right],
\quad k=\frac{2\pi}{L}(j+\eta_N),
\]

where eta_N=0 for odd N and 1/2 for even N. This follows by filling the N
lowest single-particle levels in each parity sector. The code checks these
energies against full-space ED at L=8,10,12.

For |h_z|<2 there are two simple Fermi points and the linearized dispersion
has speed v=4 sqrt(1-h_z^2/4). The supplied quasi-free mode-counting result
identifies this as c=1. Nevertheless, fitting its exact energies at
L=8,10,...,22 to the unmodified Casimir polynomial gives c_E=-0.38670 at
h_z=0.3 and c_E=1.71226 at h_z=1.1, even when the exact velocity is used.
The minimizing integer N changes irregularly with L. The resulting zero-mode
energy is also of order 1/L and is not included in that polynomial ansatz.
These values are estimator failures, not a negative central charge or a
counterexample.

The largest individual Fourier-density overlap also oscillates in this
finite-density control. At nonzero field, band curvature splits states
that become degenerate in the CFT, and the finite ring need not realize the
identity charge sector. Thus neither exact-degeneracy clustering nor choosing
one largest-overlap state guarantees the vacuum stress normalization.
Identifying the conformal family and accounting for its zero mode are additional
requirements. We do not use disagreement of these raw estimators alone to
prove a generic candidate nonconformal.

## A nonconformal control with a logarithmic entropy

For J=diag(1,1,-1), an alternating on-site pi rotation about z maps the
even periodic chain to -sum_j sigma_j dot sigma_{j+1}. This on-site unitary
preserves every block entropy. Since sigma_j dot sigma_{j+1}=2P_{j,j+1}-1,
every fully symmetric state has bond energy -1. Its total energy is -L and
there is no Casimir coefficient. Every nonzero h_n annihilates such a state.

The symmetric state with M=L/2 down spins has Schmidt probabilities

\[
p_a=\frac{\binom{\ell}{a}\binom{L-\ell}{M-a}}{\binom{L}{M}}.
\]

These follow by grouping the equal-amplitude spin configurations by the
number a in the block. Direct evaluation yields a raw chord-fit c_eff around
1.36 for L=12,...,22. This exceeds one without conformal dynamics.

On the one-spin-flip basis, H-E_0 has diagonal 4 and hopping -2 to each
neighbor, so its eigenvalue is 4-4 cos q. At q=2 pi/L this scales as
8 pi^2/L^2. The dynamical exponent is two and the inferred velocity tends
to zero. The half-chain entropy has leading coefficient (1/2) log L from the
width of the hypergeometric distribution, so that estimator would suggest
3/2. The exact energy, dispersion, and Fourier matrix elements show why
neither entropy slope is a conformal central charge.

## Search and error policy

The initial score is c_eff(L2) minus 0.6 times the entropy drift and
0.12 min(|z_gap-1|,5), where z_gap=-log(g2/g1)/log(L2/L1). A translation-defect
penalty is also included. The first nonzero gap above a fixed numerical
threshold 1e-7 is used for this screening score. Exact or exponentially small
multiplet splittings can therefore be omitted. The raw lowest gap and levels
are retained. This is a heuristic ranking, not a criticality test. In
particular, near the ferromagnetic point a magnetization-changing gap can
scale as 1/L while finite-momentum excitations remain in a quadratic
crossover regime. The detailed velocity tests are needed.

Reported extrapolation spreads are sensitivity estimates from fit windows and
powers, not rigorous confidence intervals. Eigenvector residuals quantify
numerical accuracy separately. A point is not claimed to have c>1 unless its
entropy, spectrum, and justified stress identification support a common
large-L description. Failure of that consistency test does not prove it is
gapped, nor exclude a longer crossover beyond L=22.
