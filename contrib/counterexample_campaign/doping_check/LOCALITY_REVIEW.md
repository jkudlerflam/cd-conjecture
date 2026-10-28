# Independent review of charge compactification and marker harmonics

**Status: [derived field-theory selection rule; exact finite-ring control
model checked].** The operator dictionary below assumes the standard
long-wavelength number-field representation of the stable vacancy fluid.
It does not prove that the full d=24 Hamiltonian flows to that theory.

## 1. Cyclic-color control model

Take `b` occupied colors and one vacancy. Penalize an adjacent occupied
pair unless its colors advance by one modulo `b`, and allow vacancy swaps
of amplitude `-t`. In the perfect ordered-color sector, `M` must be a
multiple of `b`; there are `b` cyclic origins. Every occupied adjacency
already has the prescribed colors, so its color penalty vanishes at all
positions. The perfect sector is invariant under hopping.

Its exact spectrum is the union of free spinless-fermion spectra with
internal twists `theta=2 pi r/b`, `0<=r<b`, including the fermion parity
shift. A particle crossing the physical cut rotates the material colors,
so the `b` origins are tied to charge flux. They are not `b` disconnected
copies of the same physical-ring spectrum.

`color_gluing.py` verifies every eigenvalue in six finite perfect-sector
rings with `b=2,3`, including prime physical lengths 5 and 7. Dimensions
range from 20 to 105; the largest discrepancy is `5.33e-15`. Every tested
fixed-M sector has one ground state, in the zero-twist charge sea. The
other origin combinations appear at charge zero-mode energies. See
[color_gluing_results.json](color_gluing_results.json).

This exact control tests gluing. It does not assert a uniform spectral
gap to imperfect material-color words. For finite positive color penalty,
all such words must be treated separately when asserting the full model's
phase.

## 2. Fractional harmonics are real and must be included

Let `j(x)` count occupied particles up to physical coordinate `x`.
Normalize the charge boson so that

    j(x) = rho x - phi_c(x)/pi + oscillatory terms,
    Delta(exp(2 i a phi_c)) = a^2 K_c.

In the cyclic-color control, a local color Fourier component contains
`exp(2 pi i j(x)/b)`. Its physical wavevector is `2 pi rho/b`, and its
charge scaling dimension is `K_c/b^2`. These fractional harmonics are
allowed local observables. Restricting the charge vertex lattice to
integer powers of `exp(2 i phi_c)` would miss them.

For a cell model with `b` occupied sites per material cell, a material
operator of momentum `q` (modulo `2 pi` in cell coordinates) can carry any
of the branches

    Q = q + 2 pi ell,  ell in Z.

Embedding at the particle positions gives physical phase and charge
vertex

    exp(i rho Q x/b) exp[-i Q phi_c/(pi b)].

Its charge contribution to scaling dimension at the decoupled point is

    Delta_c = (Q/(2 pi b))^2 K_c.

The cyclic-color example is `q=0, ell=1` and reproduces the fractional
marker harmonic exactly at the level of the number-field dictionary.

## 3. Why half filling repairs the resonance test

A physical lattice Hamiltonian can retain a uniform vertex only if

    rho Q/b = 2 pi h,  h in Z.

At `rho=1/2`, this forces

    Q = 4 pi b h,  q = 0 mod 2 pi,
    Delta_c = 4 h^2 K_c.

Thus all inherited material vertices with nonzero crystal momentum stay
oscillatory, even after every fractional marker harmonic is included.
For `h=0`, the remaining selection is exactly the internal zero-momentum
selection. For `h != 0`, the decoupled charge dimension is at least four
at `K_c=1`, and remains above two for `K_c>1/2` before accounting for
charge/internal Gaussian mixing. The primitive color harmonic is very
relevant as an observable but has nonzero physical momentum; it is not
an allowed uniform Hamiltonian perturbation at this density.

This supports the proposed half-density dictionary. It also explains why
one must not infer stability from an integer-only charge compactification:
the correct vertex set is larger, even though the final resonance bound
at this particular rational density survives.

## 4. Remaining global obligations

Exact occupied-particle conservation excludes number-changing dual-charge
vertices from the Hamiltonian. The original material model also conserves
markers. In the selected phase the relation `M=b N_markers` is an energetic
locking condition, not an operator identity throughout the microscopic
Hilbert space. One must establish a positive cost for the corresponding
mismatch defects. Otherwise an additional compressible species could
change the low-energy field content.

At physical lengths not compatible with half filling or with `M=0 mod b`,
choose the nearest allowed occupied number; the density change is `O(b/L)`.
A full finite-ring argument must still check the defect cost and the
resulting zero-mode/gluing sectors. The exact small-ring control and
fractional-harmonic calculation do not prove those properties for the
23-state material model.

Finally, the internal velocities at n=23 are extremely small. Bounding
bare charge/internal interactions by `J/t` alone is insufficient: their
size must be small in the Gaussian form normalized by the smallest
retained velocity. Choosing `J/t` sufficiently small can supply that
hierarchy at fixed n, but the local finite-J transfer requires its own
argument. No gapless charge mode has been integrated out by the finite
projection check.

Author: independent GPT-6-family research agent, 2026-10-28. The half-filling
proposal was supplied by the coordinating agent and reviewed here through
the fractional-color control; no human review is claimed.
