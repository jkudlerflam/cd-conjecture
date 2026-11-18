# Microscopic amplitudes: what is established

## 1. Full leading effective interaction

In the original spin frame Delta=1, imaginary hopping w=i/8, repulsion
V=4, bulk mu=1/2. In the controlled formal cusp expansion one instead takes
Delta/V=epsilon much less than one, with w and mu of order Delta²/V.
Let P forbid nearest-neighbour occupations. Up to a constant the leading
bulk effective Hamiltonian is

P sum_j [w b^dagger_{j+1}b_j+h.c.
 -Delta²/(2V)(b^dagger_{j+2}b_j+h.c.)
 -Delta²/V n_j n_{j+2}-Delta²/(3V)n_j n_{j+3}
 -(mu-3Delta²/V)n_j] P.

The blockade, both attractions and chemical shift must be retained in
any matching calculation. At the V4 values, |t2|=1/8, the interactions are
-1/4 and -1/12, and the shifted chemical potential is -1/4. They are not
weak perturbations to the two-valley band. The full leading model was
previously checked; see [C173](../../codex/work_C173/PROOF.md) and the
[prior V4 algebra](../v4_interpretation_audit_20261117/algebra.json).
The field-theory relative-backscattering amplitude g_s depends on this
whole interaction and dressing. Its value at V4 is not derived here;
no weak Born estimate or old V40 scattering length substitutes for it.

## 2. Leading number-changing source

The mixed second-order paths create an adjacent pair and move one endpoint.
They give separation-two amplitude -Delta(w+w*)/V. The two denominators
are equal for every compatible blockaded background; the amplitude vanishes
for the imaginary hopping of this bond. A uniform even-separation source
is also ruled out exactly by Sigma, beyond this perturbative cancellation.

At third order two adjacent pairs are created and the middle pair removed:
create (j,j+1), create (j+2,j+3), annihilate (j+1,j+2), or exchange the two
creation steps. Costs are V then 3V. Thus the complete leading pair source is

H_pair,3 = g3 P sum_j (b^dagger_j b^dagger_{j+3}+h.c.) P,
g3=2 Delta³/(3V²).

checks.py enumerates all three-pair paths on a ten-site periodic chain,
including every initial blockaded configuration. There are no additional
correlated Delta N=2 paths at this order; the matrix matches the expression
above exactly. The imaginary-hop mixed term also cancels for all these
backgrounds. This specializes [C179 sections 1--2](../../codex/work_C179/PROOF.md)
without re-running its much larger spectral extraction.

The formal V4 value is g3=1/24, or g3/|t2|=1/3. It is a **leading expansion
coefficient**, not an exact V4 matrix element. At epsilon=1/4 the omitted
virtual corrections have no established small error bound, and bulk mu/V
is 1/8. Physical boundaries have different processes and denominators.

## 3. Generated quartet: permission is not an amplitude calculation

Expanding the separation-three source in bosonic valleys gives same-valley
pair components at 2k_A and 2k_B and an inter-valley component at pi.
The latter has a form factor proportional to exp(3ik_A)+exp(3ik_B),
which is generically nonzero; all these individual pair sources oscillate
at generic incommensurate momenta. Products of two inter-valley pairs, and
products of the two same-valley pairs, carry 2pi and can generate the
quartet. Products with density harmonics generate the mixed ledger entries.

A matching calculation at second order in g3 is a connected cumulant in
**the full leading interacting Hamiltonian above**. Schematically the pi
channel contains g3² sum_r integral d tau (-1)^r times an interacting
short-distance/OPE kernel. The other pair channels and harmonics also
contribute. The alternating kernel sum can cancel; symmetry does not fix
its value or guarantee it is nonzero. Neither dropping the attractions nor
retaining only a free-band denominator derives the actual coefficient.

If that kernel is finite and nonzero with its energy scale of order |t2|,
power counting permits

Gamma4 = C4 g3²/|t2| + ...
       = C4 (16/9) epsilon² |t2| + ...
       = C4 (8/9) Delta^4/V^3 + ... .

C4 is an **uncomputed matching coefficient**, possibly suppressed or zero,
not an asserted order-one number. C179 section 5.2 wrote 16/9 also in the
last physical-units expression; that equality misses the factor 1/2 in
|t2|=Delta²/(2V). Its earlier g3²/|t2| expression has the correct power
counting. The formal V4 scale would be C4/72, not a measured coupling.
No locking length, gap or sign is predicted from it.

This low-energy generated amplitude differs from a direct Delta N=4 matrix
element of a canonical effective Hamiltonian on the *whole* blockaded
space. C179's path count places possible direct terms at sixth weighted
order or later, and its small-ring coefficients showed further cancellation.
Those results do not compute the interacting continuum C4. Conversely an
allowed continuum operator is not proof that a particular direct microscopic
coefficient is nonzero at the earliest power-counted order.

## 4. Prior claims retained, corrected, or withheld

| Source | Retained | Qualification for this audit |
|---|---|---|
| C179 | odd-separation pair source, generic quartet/backscattering ledger | quadratic charge/relative decoupling is not enforced; quartet matching coefficient uncomputed |
| C180 | exact symmetry stratum, generic integer selection, allowed drift, zero-drift dimensions | drift covariance fails a free control; fermionic resonance labels cannot use bosonic peak centres |
| C189 | bosonic operator frame and pair/quartet relation | density Ward identity applies to its conserved effective charge, not raw V4 density; old fitted K not imported |
| C191 | chiral Gaussian covariance and number-breaking threat | its other-density/V estimates and heuristic locking lengths do not determine V4 |
| Previous V4 audit | two symmetry-related peaks and parity/accuracy-sensitive entropy | no new exponent constraint or central-charge determination follows |

The possibility of emergent U(1) despite microscopic pairing is real in
other constrained chains, as demonstrated in
[Verresen, Vishwanath and Pollmann, opened](https://arxiv.org/abs/1903.09179).
It is not a guarantee here: the dimensions and amplitudes of this bond's
allowed perturbations have to be established independently.
