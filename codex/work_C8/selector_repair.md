# A uniform selector-defect gap from a local unitary repair

Status: proved for the binary-arrow selector with a strengthened
penalty. This changes only a forbidden-sector coefficient of the
existing explicit C8 bond. The entire legal Hamiltonian, including
all finite-coupling effects, is unchanged.

## 1. Binary arrows and the repair unitary

The selector alphabet is the set of all two-bit pairs (a,b):

    #=(1,0), L=(0,0), X=(0,1), R=(1,1).

A bond is legal when its outgoing bit b_i equals the next incoming
bit a_(i+1). Its violation projector is

    F_i=1_(b_i != a_(i+1)).

All remaining onsite factors are spectators for this definition.
Let X_in,j flip a_j without changing b_j or any other factor. In
letter notation it swaps # with L and X with R.

[proved] For any definite defect pattern f_i in {0,1}, the unitary

    U_f=product_(i:f_i=1) X_in,(i+1)

maps that sector onto the completely legal sector. It changes F_i
to 1-F_i at each selected bond and leaves every other F_j invariant.
This follows because each incoming arrow occurs in just one F_i.
There is no global parity obstruction. Even an isolated single
forbidden bond can be repaired. The map need not preserve marker
number, which is appropriate when comparing the full ground state.
It does preserve any charge acting only on the nonselector factors.

## 2. An exact energy comparison with no volume error

Write a translation-invariant nearest-neighbor ring Hamiltonian as

    H=K sum_i F_i+sum_i g_(i,i+1),  ||g||<=B.

Assume H preserves every F_i, as the C8 bond does. Let E_legal be
its minimum over all zero-defect states, including marker-free ones.

[proved] In a sector with k=sum f_i defects,

    H|f >= [E_legal+(K-4B)k] I.                           (1)

For a normalized state psi in that sector, U_f psi is a legal trial
state. At most 2k bonds meet the k sites changed by U_f, and on each
such bond ||U_f dagger g_i U_f-g_i||<=2B. Therefore

    E_legal <= <U_f psi,H U_f psi>
             <= <psi,H psi>-Kk+4Bk.

This proves (1), including for entangled states and arbitrary
spectator dimension. No gap in the legal Hamiltonian and no estimate
of its extensive ground energy are needed. K>4B gives a positive
defect gap uniform in ring length.

The conserved-defect assertion for C8 can be read directly in arrow
bits. Its nontrivial selector changes are LX <-> XR and #L <-> R#.
Each flips both inner arrows b_i and a_(i+1), keeping their mismatch
and both outer arrows unchanged. All other bond operations leave the
selector letters fixed. Hence every overlapping F_j is conserved on
the full physical Hilbert space, not just on promised legal words.

## 3. An explicit strengthened C8 coefficient

Use the flat wrapper, unit U and J0, and the original core at unit
scales. For 0<=epsilon<=epsilon_max<=1 the bond other than K F obeys

    B <= B0+epsilon_max B_flat,
    B0=|mu_selector|+t+49,
    B_flat=1+2m+4(k+2)+409600 pi^4 b^4 m^3+320 pi^2 b^2 m.  (2)

All quantities are those in periodic_clean_bond.py. The marker density
has norm at most |mu_selector|, head hopping at most t, the complete
flat wrapper bond at most 45 by flat_vacancies.md, and the two internal
helix terms and two cut pins cost at most four in a triangle bound.
Equation (2) uses the full-space perturbation bound, not projected
visibilities.

[proved] Replace the old selector coefficient by

    K_cert=max(K_old,8[B0+epsilon_max B_flat]).            (3)

Then (1) excludes every forbidden selector sector by at least
4[B0+epsilon_max B_flat] per defect throughout that coupling interval.
This replacement is a single repeated diagonal two-site term. It
does not change the on-site dimension, symmetries, any legal-sector
matrix element, the clock, wrapper, atomic gap, or projected target.
In particular none of the legal estimates has to be rederived.

The executable implementation is CertifiedPeriodicBond in
certified_periodic_bond.py. It subclasses the existing sparse oracle
and changes only Ksel. Its formula is the specification at arbitrarily
large finite q; numerical construction remains intended for modest q.
The coefficient is fixed before L varies. No claim is made that the
old smaller coefficient already satisfies this new proof.

## 4. Marker-free legal rings

[proved] A legal selector ring with no # is either all L or all R.
There cannot be an X, since an X changes arrow 0 to 1 and only a #
can change it back. Every perturbation term requires a marker, so
the perturbation vanishes in these two selector sectors. There is no
head hop. The periodic flat-wrapper Hamiltonian and helix terms are
nonnegative by the inherited complete-sector bounds. Thus their
absolute energy is at least zero.

The legal positive-marker trial of marker_density_bounds.md has
absolute energy at most e_base L+epsilon B_L, with the fixed
e_base=t e_c<0. It is negative for all sufficiently large L.
Consequently marker-free sectors are also excluded. Combined with
(1)-(3), the full microscopic ground state of the strengthened bond
lies in the positive-marker legal sectors at large L.

The density and factorial-tail estimates now concern the actual
full ground sectors, subject to their stated small-coupling bounds.
At any fixed low-energy window, illegal and marker-free sectors are
absent once the window is below their respective gaps. This does not
identify the infrared theory inside the legal sectors. The all-order
low-frequency elimination and Gaussian infrared identification remain
unproved.

## 5. Checks

selector_repair_checks.py exhausts binary-arrow repair on short rings,
tests the energy inequality in quantum selector-plus-spectator chains,
and checks the strengthened full-space sparse columns. Its random
quantum control preserves each defect pattern but has a general
complex interaction within every permitted block. It is not a
classical-energy check. Numerical sizes, norm comparisons, and the
explicit q=4 coefficient are recorded in selector_repair_checks.json.
