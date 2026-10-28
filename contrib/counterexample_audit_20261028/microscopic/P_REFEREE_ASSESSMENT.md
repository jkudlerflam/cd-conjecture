# Fresh assessment of the referee's input P

**Assessment:** the referee correctly identifies the remaining continuum
inference as conventional physics rather than a consequence of BDL alone.
After rechecking that inference, I do not find an additional unverified
bond-specific assumption that the frozen chain has c=23. Its microscopic
reduction and the hypotheses of the usual weak-coupling Gaussian RG
argument can be exhibited separately. This is not a mathematical
lattice-to-CFT theorem, and it does not settle vacancy matching M.

This is a fresh response to `contrib/referee_d24/README.md`, written by the
microscopic-audit GPT-6-family worker on 2026-10-28. The prior audit verdict
is not used as a premise. No human check or measured central charge is
claimed.

## What P actually needs

1. **An actual local reduction, with a massive complement.** In the selected
   marker sector the physical bond blocks into finite cells, with exact
   rank-2^23 code projectors and gap Delta_0>0. The port is bounded and
   nearest-cell. The data field belongs at order lambda^2. The local BDL
   construction gives an exact unitary identity with a full-space,
   exponentially local remainder O(lambda^3), at fixed order two. The
   block-diagonal reference has a threshold Delta_0/2 above its ground
   energy in the auxiliary complement. These are microscopic statements,
   not assumptions about the infrared central charge. The applicable
   source is [BDL, Sections 4.1–4.4](https://arxiv.org/pdf/1105.0675).

2. **A genuine free-fermion reference with a finite low-energy window.**
   The leading g^2 source contraction is quadratic. C33's rational
   certificate gives 46 simple roots of its quadratic projection,
   together with strict bilinear and quartic mismatch margins. Simplicity
   implies a nonzero velocity and a single null eigenvector at each root.
   Compactness of the Brillouin zone then supplies a positive energy
   cutoff outside sufficiently small, disjoint root neighborhoods. No
   numerical lower bound on that cutoff is required for an existence
   argument, although its size enters the unknown constants.

3. **A small genuine interaction, not just a small quadratic error.**
   The complete two-cell resolvent has an operator-norm remainder
   R_g bounded by alpha^4 g^4/4 after removing its g^0 and g^2 terms.
   The nonquadratic residual of the Hilbert–Schmidt quadratic projection
   is also small, with a fixed projection constant. Explicitly, on
   m=2n fermion modes, the off-diagonal quadratic coefficients of R obey
   |h_ij|<=||R||, the diagonal coefficients obey |h_ii|<=2||R||, and
   the constant coefficient has magnitude at most (m+1)||R||. Therefore

       ||R - quadratic_projection(R)||
           <= [1+(m+1)^2] ||R||.

   At n=23 this factor is 2210. It cannot invalidate the limit: after
   normalizing by the smallest g^2 payload scale, this bound tends to zero
   as g^2. Reducing g preserves the certified root brackets when the
   common field is rescaled with g^2.

4. **The correct operator algebra and exhaustive dangerous-vertex test.**
   Exact internal U(1) becomes total logical fermion number. A neutral
   local spin operator is an even finite polynomial in Jordan–Wigner
   fermions: exterior parity strings cancel, and any remaining string
   stays within its support. Thus fractional fermion-number vertices or
   standalone disorder strings are not missing local Hamiltonian terms.
   In the free chiral theory a primary with integer charge vector z has
   dimension sum(z_i^2)/2. Neutral primaries of dimensions one and two
   are precisely the bilinear and quartic classes checked by C33.
   Pair sums with replacement are included. Surviving zero-momentum
   quartics have identical creation/annihilation multisets and reduce to
   products of Abelian chiral densities. All other neutral primaries
   have dimension at least three; descendants add derivatives.

5. **Control throughout the finite matching interval.** Momentum mismatch
   does not mean that a detuned operator has no effect at short scales.
   It means its potentially dangerous flow stops after a finite matching
   interval. That interval, the band projectors, curvature bounds and
   all conversion factors are fixed and finite at n=23. Choosing g first
   makes the dimensionless interaction arbitrarily small throughout it.
   Subsequent root/filling shifts must be smaller than the certified
   mismatch margins; this is another finite tolerance, not an omitted
   symmetry protection. With g fixed, taking lambda smaller makes
   C lambda^3/(lambda^2 w_g)=C lambda/w_g arbitrarily small as well.

The fifth point would have been a real gap if the argument only compared
the coupling with the largest hopping, or if its matching interval grew
without bound as the interaction was reduced. Neither occurs: the
normalized g^2 band geometry is fixed, and the order of choices is g,
then lambda. The all-length selector inequalities persist as lambda is
reduced. There is no contrary lower bound on either coupling.

## What remains conventional RG input

After the finite matching interval, the dimension-two interactions are
Abelian current products. Retaining them gives a quadratic bosonic
Hamiltonian. Sufficiently small corrections relative to the smallest
velocity keep its stiffness positive. Non-Gaussian primaries begin at
dimension three, and a relative change of the scaling-dimension form
smaller than one third keeps them above two (the campaign uses one
quarter). Higher-charge vertices are protected by the same relative
bound, not by a finite enumeration of all orders. Irrelevant terms can
renormalize the Gaussian coefficients but do not supply a missing
non-Abelian marginally relevant channel.

The conventional inference is that this open Gaussian phase persists
under sufficiently small symmetry-preserving exponentially local
perturbations, including coupling to the massive auxiliary sector.
BDL does not prove that continuum conclusion. It does, however, supply
the actual local-smallness statement to which it is applied. No exact
all-order diagonalization of the extensive code band is needed.

I would therefore describe P as **a derived microscopic reduction
followed by the stated conventional weak-coupling RG inference**, not
as an independent assumption that the frozen bond happens to have
c=23. If the desired standard requires a mathematical CFT-limit theorem,
P remains unproved. At the repository's stated theoretical-physics
standard, I found no further bond-specific missing input in P.

## Quantitative hierarchy and the referee's checks

The hierarchy criticism is substantial for feasibility and confidence.
The smallest second-order velocity has log10 approximately -431.692
before the factor lambda^2. Applying the conservative factor 2210 above
to the reported full resolvent remainder gives log10 approximately
-442.563, or about 1.35e-11 relative to that velocity at g=2^-365.
This numerical comparison is illustrative; it does not compute the
complete matching constants. One may reduce g further if needed.

For the BDL remainder the condition is

    lambda << w_g/C_2,
    log10(w_g) approximately -431.692 at g=2^-365.

Thus a scale such as lambda near 10^-435 is compatible with the direction
of the hierarchy, but **cannot be certified sufficient without C_2 and
the RG tolerance**. An even smaller fixed positive lambda is allowed.
The existence argument establishes a nonempty asymptotic family; it
does not provide one numerically certified final parameter tuple.

A human audit of the highly conditioned C33 certificate and a calibrated
small-n central-charge measurement would materially increase confidence.
They are reasonable acceptance requirements for a referee. They are not
additional physical hypotheses in the weak-coupling argument. Conversely,
a successful n=1 or n=2 central-charge measurement would not itself
validate the n=23 conditioning, mismatch margins or scale hierarchy.

The independent literal-cell g^2 calculation in this audit already tests
the bond-to-source-formula identification at n=1,2, and the n=23 integer
certificate was regenerated and its acceptance logic inspected. Those
are meaningful checks but are not a replacement for either human
scrutiny or a central-charge measurement.

**Conclusion for the coordinator:** I found no new false lemma or hidden
bond-specific obstruction in P. Keep its continuum step explicitly
labelled conventional RG, keep the final couplings existential rather
than numerically certified, and do not use this assessment to relabel
the separate gapless vacancy-matching input M as a theorem.
