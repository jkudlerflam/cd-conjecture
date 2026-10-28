# A momentum-restoration threshold for the half-filled free charge chain

Status: [proved] free-fermion inequality and its thermodynamic
application. This is a useful input to vacancy matching, not a proof
of the interacting matching theorem.

## 1. Arbitrarily many excitations

Let epsilon(k)=-2t cos(k), with chemical potential zero and Fermi
points +/-pi/2. For k on the Brillouin circle write

    k=s pi/2+d,     s in {+1,-1},     |d|<=pi/2,

choosing a closest Fermi point. Since sin(x)>=2x/pi on [0,pi/2],

    |epsilon(k)|=2t sin|d| >= (4t/pi)|d|.               (1)

Consider any number-conserving collection of particles outside and
holes inside the half-filled sea. Its excitation energy is the sum
of |epsilon(k)| over the changed occupations. Its total momentum has
the form

    P=pi m+sum_(particles) d-sum_(holes) d,   m in Z.   (2)

Indeed the Fermi-point contributions are (pi/2) times the difference
of two sums of equally many +/-1 signs, hence an integer multiple of
pi. Combining (1), (2), and the triangle inequality on the circle gives

    E_exc >= (4t/pi) dist(P,pi Z).                     (3)

This is independent of the number of particle-hole pairs. It permits
zero threshold at P=0 or pi modulo 2pi, which are the usual forward
and 2k_F channels. It forbids a zero threshold at other fixed physical
momenta. The bound is deliberately weaker than the exact lower edge.

## 2. Internal twists do not evade physical momentum conservation

On a ring of L sites, the allowed charge momenta in an internal
translation sector K are

    k_j=(2pi I_j+K+theta_parity)/L.

Changing K by Q shifts a continuously followed sea's total physical
momentum by N Q/L. Changing the sea's current branch adds
2pi N m/L, while its energy changes by O(t/L) for fixed branch m.
At N/L=1/2, these statements become Q/2 modulo pi and O(t/L).

An interaction preserving physical lattice momentum cannot connect
the two twisted sea ground states alone when Q/2 is nonzero modulo
pi. Its charge part must supply the missing physical momentum.
Applying (3) to that excitation gives

    liminf_(L->infinity) E_charge,needed
      >= (4t/pi) dist(Q/2,pi Z).                      (4)

The limits use N/L->1/2 and bounded total excitation energies. For
N=L/2 on compatible even rings, the exact finite-grid sea may differ
from the chemical-potential-zero sea by at most endpoint occupations;
their energies and momenta approach zero offsets. For the actual
N rounded to multiples of b, N-L/2=O(b), the same endpoint correction
vanishes as L grows at fixed b.

Allowing the current branch m to grow with L also does not evade (4).
A bounded-energy half-filled current shift has |m|=O(sqrt(L)); its
extra momentum caused by N/L-1/2=O(1/L) is O(1/sqrt(L)). It vanishes.
Equivalently, the elementary finite-grid energy bound controls the
sum of distances of changed occupations from the limiting Fermi
points. The O(1/L) filling offsets cannot generate a fixed nonzero
restoration momentum at bounded vanishing energy.

This argument is special to the rational half-filled charge sea.
At an irrational limiting density, increasingly large current
windings can approximate additional momenta. One must not transfer
(3) with pi Z to an irrational density without further analysis.

## 3. The candidate's uniform positive threshold

For a payload momentum q in cell units and cell harmonic ell,

    Q=(q+2pi ell)/b mod 2pi.

If q is separated from 2pi Z by delta_int, then

    dist(Q/2,pi Z)
      = (1/(2b)) dist(q+2pi ell,2pi b Z)
      >= delta_int/(2b).

Therefore the charge channel needed by any such nonforward low-order
payload process has threshold

    E_threshold >= 2t delta_int/(pi b).               (5)

At b=24 and delta_int>3.1558e-7 this is greater than
8.371e-9 t. Pure nontrivial period-b order harmonics have the larger
bound 4t/b. The q=0 uniform component and the ordinary soft charge
channels remain in the low-energy action, where their forward and
umklapp operator dimensions must be examined separately.

Equations (3)-(5) repair the incorrect inference that a change of
internal twist itself costs a finite charge energy. It costs a finite
energy only when the full process must restore physical momentum and
the twist is nonresonant. Small nonzero J requires a local matching
argument around this free threshold; the free inequality alone does
not establish the dressed one to all orders.
