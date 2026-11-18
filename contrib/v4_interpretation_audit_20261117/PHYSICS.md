# Exact symmetry and conditional effective description

## Microscopic symmetry

A uniform rotation U=exp[-i(pi/8) sum Z] takes

H0 = sum [(XX-YY)/2 + (XY-YX)/16 + ZZ] - (7/4) sum Z

to the displayed 9XY/16+7YX/16 model. Set n_j=(1-Z_j)/2; b_j annihilates a
down-spin occupation. Pair terms violate N=sum n conservation. P=(-1)^N
is exact.

In the original frame set R=exp[-i sum_j (pi/2+pi*j)n_j], with zero-based j,
and A=R K, with K complex conjugation in the Z basis. Then A is antiunitary,
A^2=1 and A H A^-1=H. Adjacent pair terms have R phase one; imaginary
hopping changes sign under both R and K. On even periodic rings,
A T A^-1=P T. Hence k maps to -k in the even sector and pi-k in the odd
sector, modulo 2pi.

The observed L=16 even +/-pi/4 partners and odd pi/4,3pi/4 partners are
therefore consistent with exact symmetry partners, not extra mode counts.
The uniform rotation carries A to the displayed frame. algebra.py directly
checks these identities and the nonzero [H,N]. The XX+YY bond is A-odd;
its zero expectation in a symmetry-preserving finite state does not exclude
spontaneous order. Its connected channel was measured in the preceding audit,
without a converged order parameter or infrared exponent.

## Complete leading projected model

Original-frame Jordan-Wigner parameters are Delta=1, w=i/8, V=4,
mu_bulk=1/2, and mu_endpoint=mu-V/2=-3/2. See
[C190 section 2](../../codex/work_C190/PROOF.md).

Expanding about V sum n_j n_j+1 and projecting away adjacent occupations,
the leading bulk Hamiltonian, after dropping a constant, is

\[
P\sum_j\left[\left(\frac{i}{8}c^\dagger_{j+1}c_j+\mathrm{h.c.}\right)
-\frac18(c^\dagger_{j+2}c_j+\mathrm{h.c.})
-\frac14n_jn_{j+2}-\frac1{12}n_jn_{j+3}+\frac14n_j\right]P.
\]

P enforces the nearest-neighbour blockade. The hard-core-boson form has
the same bulk projected hopping; periodic Jordan-Wigner boundary signs
must retain parity. The chemical shift is mu_tilde=mu-3Delta^2/V=-1/4.
In units |t2|=1/8, hopping strengths are i and -1, attractions -2 and -2/3,
and chemical potential -2. **Both attractions and the blockade are leading**,
not optional corrections. The direct eight-site calculation verifies the
second-order pair contribution, including constant -L/4, shift 3N/4, both
attractions and hopping. See [C173 section 3](../../codex/work_C173/PROOF.md).

This is a formal leading expansion, not an error-controlled description of
V4. Delta/V=1/4, and denominator corrections and higher-order number-breaking
terms are not bounded here. Near-ferromagnetic leg interactions can produce
slow collective velocities without extreme microscopic ratios; bare band
curvature is not a measured velocity. Open endpoints have different virtual
processes, so the bulk expression cannot replace their exact Hamiltonian.

Raw density includes virtual pairs and differs from the conserved effective
charge. Vacuum dressing estimates cannot be subtracted at density about
0.29 without a controlled calculation. K and velocities measured at V=40
or another chemical potential cannot be transplanted to this point.

## Conditional two-valley fields and locking

For c_j proportional to exp(ikj), the leading free dispersion is
 epsilon(k)/|t2|=2sin(k)-2cos(2k). Its minima are -s and -pi+s,
s=asin(1/4), with folded separation 0.83914pi and sum -pi. The measured
bosonic separation is about 0.882pi; interactions and density matter.
Extracting a renormalized hopping from that separation alone is unjustified.

If b(x) contains exp(ik_A x) Psi_A + exp(ik_B x) Psi_B, density contains
interference at k_B-k_A and intravalley harmonics at 2pi*rho_A,B. This
accounts for a density peak at the valley separation without calling it
the fundamental 2pi*n harmonic of a one-component liquid. It does not
establish that both components are independently gapless.

For a hypothetical Gaussian fixed point define
phi_c,s=(phi_A +/- phi_B)/sqrt(2), and similarly theta_c,s. Relative
backscattering cos(sqrt(8)*phi_s), dimension 2K_s, can gap that sector if
relevant. Psi_A Psi_B carries momentum pi, whereas its square carries
reciprocal-lattice momentum and is allowed by translation and parity.
The quartet vertex exp(2i*sqrt(2)*theta_c), dimension 2/K_c, can lock the
total charge if relevant and nonzero, despite the apparent U(1) of the
leading projected model. These are conditional dimensions, not measured
K values or claims that either flow occurs at V4. An actual proposed fixed
point also requires checking any other commensurate resonances.

See [C189](../../codex/work_C189/PROOF.md) for operator-frame and emergent-
charge cautions, and [C191 section 8](../../codex/work_C191/PROOF.md) for the
allowed quartet. Their old locking-length estimates at different parameters
are not imported. One surviving Gaussian sector gives c=1; two stable ones
give c=2. A long crossover may resemble the latter while the ultimate limit
is different. Neither that flow nor any additional sector is resolved here.

## Simplicity accounting

Physical d=2, without markers or blocking. The displayed basis has four
real coefficients and four literal Pauli monomial types; after removing a
common positive scale there are three independent ratios. Largest/smallest
absolute displayed coefficient ratio is 4. The original basis has six
monomial types and ratio 28; notation does not change the physical model.
Pairing-to-imaginary-hopping ratio is 8 and ZZ-to-DM coefficient ratio 16
in these conventions. Few displayed terms do not remove strong effective
interactions or an emergent slow sector.
