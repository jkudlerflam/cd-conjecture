# A finite-energy charge threshold for nonforward material transfer

**Status: [proved free-charge spectral bound; independently checked by
exact momentum optimization].** The theorem resolves a specific concern
about fast-charge matching: different material momentum sectors impose
charge boundary twists whose unconstrained sea energies differ only by
`O(1/L)`, but those seas have different **physical momenta**. A uniform
microscopic interaction cannot connect them without compensation.

The lower bound covers arbitrarily many particle-hole excitations. It
does not assume that the material-coordinate adjacency operator is local
in physical coordinates, or approximate that operator by a density mode.

## 1. Exact momentum bookkeeping

Use the definitions and signs in [README.md](README.md). Set `L=2M`, so
the occupied density is exactly one half. In the internal cyclic sector
`R=exp(i theta_r)`, `theta_r=2 pi r/M`, charge orbitals have momenta

    k_n = [2 pi n - theta_r - pi(M-1)]/L,
    epsilon(k) = -2t cos k,  n=0,...,L-1.

For any set `I` of `M` occupied orbitals, the physical translation
phase of the complete state is

    exp[(2 pi i/L)(r+M(M-1)/2-sum_(n in I) n)].

A charge ground sea occupies `M` consecutive integers starting at `ell`,
where `ell` minimizes `|theta_r-2 pi ell|`. Its physical momentum label is
`r-M ell` modulo `L`. The two values of `ell` at a half-flux degeneracy
are retained as separate momentum branches.

Consider an initial charge ground sea in sector `theta` and a target
sector `theta'`. Relative to a target ground sea, a target excitation
carries charge momentum `P`. Equality of the complete physical momenta
forces, modulo `2 pi`,

    P = (theta'-theta)/2 - pi(ell'-ell).

In particular,

    dist(P, pi Z) = (1/2) dist(theta'-theta, 2 pi Z).

The `O(1/L)` difference between the unconstrained sea energies therefore
does not by itself provide a low-energy channel at conserved physical
momentum.

## 2. Bound for every charge excitation

At half filling, a ground sea fills every negative-energy orbital and
chooses the required occupancy among any zero-energy orbitals. Every
number-preserving excitation can be described by particle orbitals `p`
and hole orbitals `h`, with the same number in both sets. Its energy
above that sea is exactly

    Delta E = 2t [sum_p |cos k_p| + sum_h |cos k_h|].

Let `d(k)` be the distance on the momentum circle to the two Fermi edges
`{pi/2,-pi/2}`. Then `0<=d<=pi/2` and

    |cos k| = sin d(k) >= (2/pi) d(k).

Write each momentum as a nearest Fermi edge plus a deviation of magnitude
`d(k)`. Because particle and hole numbers agree, the signed sum of all
chosen Fermi-edge momenta is an integer multiple of `pi`. The triangle
inequality on the circle therefore gives

    dist(P, pi Z) <= sum_p d(k_p) + sum_h d(k_h).

Combining the inequalities yields the exact finite-ring bound

    Delta E >= (4t/pi) dist(P, pi Z).

This holds for **every occupation set**, including states containing an
extensive number of particle-hole excitations. Thus higher charge
harmonics or the nonlocal material-coordinate string cannot evade it.

## 3. Comparing different twists

The exact charge ground energy at half filling is

    E_sea(theta) = -[2t/sin(pi/L)] cos(delta(theta)/L),
    delta(theta) = theta reduced to [-pi,pi].

Its range over twists is

    max E_sea - min E_sea = 2t tan(pi/(2L)).

Consequently every target state in sector `theta'` with the same physical
momentum as an initial ground sea in sector `theta` satisfies

    E_target - E_initial
      >= (2t/pi) dist(theta'-theta,2 pi Z)
         - 2t tan(pi/(2L)).                         (1)

For any fixed nonzero twist transfer, the thermodynamic lower bound is
strictly positive. The exceptional transfers are precisely zero modulo
`2 pi`, including the usual soft current/umklapp branches already present
in the forward sector.

For a material cell with `b` occupied sites, cell crystal momentum `q`
has one-site material branches

    theta'-theta = (q+2 pi ell)/b,  ell in Z.

Minimizing (1) over every branch yields

    E_target - E_initial
      >= [2t/(pi b)] dist(q,2 pi Z)
         - 2t tan(pi/(2L)).                         (2)

Thus a fixed finite set of nonforward material vertices with mismatch
margin `delta` has a charge threshold at least
`2t delta/(pi b)` in the limit. The cell-size factor must be retained.
This is the relevant threshold for virtual material momentum transfer,
rather than the unconstrained boundary-twist splitting.

## 4. Exact optimization checks

`charge_threshold.py` uses dynamic programming to minimize charge energy
at fixed particle number and exact physical momentum, over **all** orbital
occupation sets. For `L=16` it independently enumerates every occupation
to verify the dynamic program. The larger calculations run through
`L=256`; the largest search represents `choose(256,128)` states without
truncating excitation number.

At `L=256`, `t=1`, transfers `Delta theta/pi = 1/4,1/2,3/4,1` have
fixed-momentum excitation costs approximately

    0.76750, 1.42085, 1.85844, 2.01227,

while their proved lower bounds are

    0.48773, 0.98773, 1.48773, 1.98773.

All twenty cases satisfy (1). The observed limiting values agree with
`2t sin(dist(Delta theta,2 pi Z)/2)`, the sharper free-fermion lower edge;
that equality is numerical evidence here, not needed for the proved
linear bound. See [charge_threshold_results.json](charge_threshold_results.json).

## 5. Scope for Wilsonian matching

This bound supports eliminating **nonforward charge excitations at a
fixed mismatch**. The retained forward charge boson remains gapless and
must stay in the effective theory. Vertices whose material momentum
approaches zero are part of that low-energy sector, not covered by a
uniform positive threshold.

The theorem is for free charge at exactly half filling and an initial
charge sea. Low-energy oscillator momenta and energies produce the
expected vanishing corrections to a fixed mismatch. The bounded-density-rounding extension is proved below. Interacting charge
and full local effective-action matching still require an argument; neither
follows from a finite-volume global Schrieffer-Wolff gap.

At irrational density, arbitrarily high current harmonics have momenta
`2 pi rho w` dense modulo `2 pi`. The clean global half-filling bound
should not be asserted there without restricting the operator class or
controlling harmonic dimension. Rational half filling is doing substantive
work in this construction.

## 6. Bounded density rounding on arbitrary physical rings

The same positive limiting threshold holds when `M=L/2+a`, with `a`
bounded independently of `L`; `L` may be odd. This covers rounding to a
nearby allowed particle number divisible by a fixed material cell size.
Set `A=ceil(|a|)+1`.

In a target twist sector let `F_0` fill every strictly negative-energy
orbital, and let `S` be an M-particle ground sea. Let `D_S` be the sum of
distances to `{+pi/2,-pi/2}` over the orbitals in their symmetric
difference. The sets differ by at most `A` orbitals. The nearest available
orbital on either side of an edge is within `2 pi/L`; the m-th one is
within `2 pi m/L`. Therefore

    D_S <= pi A(A+1)/L <= 2 pi A^2/L.

For any M-particle target occupation `I`, define `D_I` similarly against
`F_0`, and let `W_S` be the corresponding sum of `|cos k|` for `S`.
Both I and S have the same net particle number relative to F_0. Their
signed nearest-edge contributions consequently differ by an integer
multiple of pi. Thus for their momentum difference P,

    dist(P,pi Z) <= D_I + D_S,
    E(I)-E(S) >= (4t/pi) D_I - 2t W_S
               >= (4t/pi) dist(P,pi Z) - (4/pi+2)t D_S.

Here `W_S<=D_S` was used. This is an exact finite-ring inequality, whether
F_0 contains floor(L/2) or ceil(L/2) particles.

Reduce the initial and target twists separately to `[-pi,pi]`. Their
sea momentum difference is rho times the difference of those reduced
twists, modulo `2 pi`. Since the reduced difference has magnitude at most
`2 pi`, its distance to `pi Z` is at least

    (1/2) dist(theta'-theta,2 pi Z) - 2 pi |a|/L.

At general M the sea-energy range is
`2t sin(pi M/L) tan(pi/(2L))`, bounded by the same half-filled range.
Combining the estimates gives

    E_target-E_initial
      >= (2t/pi) dist(theta'-theta,2 pi Z)
         - t[8|a|+(8+4pi)A^2]/L
         - 2t tan(pi/(2L)).                         (3)

The remainder vanishes uniformly for bounded a. Minimizing over the b
material branches again divides the leading mismatch by b. Although
conservative and sometimes negative at small L, (3) provides a positive
thermodynamic threshold for every fixed nonforward cell momentum on all
bounded-rounding ring sequences.

`rounded_threshold.py` checks both the displayed uniform bound and the
sharper evaluated D_S bound by exact all-occupation optimization. Fourteen
cases include odd physical sizes and occupied numbers divisible by 24.
All pass; see [rounded_threshold_results.json](rounded_threshold_results.json).

Provenance: the coordinating agent proposed the sine-distance inequality;
an independent GPT-6-family worker completed the finite-ring twist and
momentum proof and exact all-occupation checks, 2026-10-28. No human review
is claimed.
