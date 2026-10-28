# Why the tiny loop produces a visible momentum mismatch

[checked] There is no normalization inconsistency. The decisive point is
that **all 46 Fermi crossings have extremely small physical slopes**, even
though the hopping matrix has a much larger operator norm. Comparing the
loop with the largest hopping entry therefore gives the wrong estimate of
the root displacement.

The source formula is

    T/g^2 = -4 alpha^4 Re(U^dagger R V),
    mu/g^2 = -4i alpha^4 Im(U^dagger R U + V^dagger R V).

After dividing by E_unit=2 alpha^4 beta^2 g^2, the saved integer matrices
represent T=-F/(2S^2) and mu=iM/(2S^2). Direct evaluation of these ratios
from the saved integer numerators gives:

| Quantity | Value |
| --- | ---: |
| beta | 9.0865873e-73 |
| t_d/beta = 2^-136 | 1.1479437e-41 |
| t_d | 1.0431e-113 |
| ||T/E_unit|| | 4.4187664e-3 |
| ||mu/E_unit|| | 7.3370954e-64 |
| ||mu|| / ||T|| | 1.6604397e-61 |

Thus the onsite loop contribution is actually *much smaller* than the
normalized hopping norm. Its norm divided by t_d/beta is only
6.3915e-23, reflecting the remaining propagation factors in the source
contractions. The explanation is not an unexpectedly large physical
onsite term.

The tree filter Gram matrix has about 60 orders of spectral conditioning:
the certificate places its largest eigenvalue near 10^-2.35 and its
smallest between 10^-63.62 and 10^-62.25. In the stable congruence used for
root certification, the same matrices instead have

    ||T_whitened||  = 0.52970087,
    ||mu_whitened|| = 0.00084713536.

This is a congruence, not a unitary change of basis. It preserves roots
and inertia, but physical eigenvalue slopes require the accompanying
positive metric. If h_tilde=C h C^dagger, B=C C^dagger and u is a normalized
null vector of h_tilde at a root, then

    physical slope = (u^dagger h_tilde' u)/(u^dagger B u),
    physical onsite expectation = (u^dagger mu_tilde u)/(u^dagger B u),
    delta k = -(u^dagger mu_tilde u)/(u^dagger h_tilde' u).

With mu set to zero while retaining the same real hopping and common
field, every center is zero to 1.2e-16. The physical root velocities of
that reference range only from 4.43e-63 to 3.76e-62 in E_unit units.
There are no Fermi crossings whose velocity is comparable to ||T||.
The band eigenvectors change across the Brillouin zone; a large bandwidth
does not imply a large slope at a crossing of this tiny common field.

For example, at negative root 7 the normalized physical slope is
-2.2795711e-62 and the loop onsite expectation is -1.3299898e-66.
Their ratio predicts a displacement of -5.8343863e-5. The numerical
finite-loop pair center is -5.8341350e-5. The large metric factor at this
root is 4.3941391e61; the corresponding whitened slope and onsite
expectation are -1.0016752 and -5.8441603e-5.

The actual stored root centers were checked, including the pair-transfer
witness. Numerical centers of the finite-loop g^2 reference are:

    center 7 = (k_7+k_38)/2 = -5.8341349967e-5,
    center 8 = (k_8+k_37)/2 = -5.8141530700e-5.

Their pair sums differ by 3.9963853e-7. More importantly, the exact rational
root brackets in `root_certificate_exact_n23.json`, which also enclose
the finite-g quadratic projection, imply the following intervals
(displayed here in rounded decimal form):

    center 7 in [-5.8351350e-5, -5.8331350e-5],
    center 8 in [-5.8151531e-5, -5.8131530e-5].

The intervals remain separated. Across all 23 centers the minimum interval
separation is about 1.79819e-7; the independent all-pair chord test gives
the more conservative certified delta_4>3.1558764e-7. The numerical centers
range from -5.8341350e-5 to +1.2016796e-4.

`check_loop_scale.py` reads the saved integer matrices and exact root
brackets, reconstructs these scales, and compares the onsite-free and
finite-loop roots. Its output is `loop_scale_results.json`. The velocity
and linear-response numbers are numerical diagnostics; the interval
certificate supplies the strict nonresonance statement. No full
certificate regeneration was required for this explanation.
