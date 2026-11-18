# Margins, physical errors and Hamiltonian complexity

## Established object

The all-n object is a finite-dimensional, leading lambda^2 g^2 quadratic
logical symbol plus a physically matched order-lambda^2 g^2 L1 field.
Its carrier has 2n+2 vertices but its **physical on-site dimension is nine**.
No marker, encoded replacement dimension, extra rail label or vacancy is
used. The statement concerns every finite n with its own finite parameters.
It does not take n to infinity inside one fixed Hamiltonian.

The original nine-state contribution certifies exactly
n=9,10,11,12,14,16,20. Those seven certificates remain evidence for their
original couplings. This contribution establishes an analytic family at
modified couplings; it does not retroactively certify an all-n version of
the old parameter schedule. Its new proof also resolves the moment
nonvanishing obstruction in the eleven-state analysis, but only after that
analysis's distinct physical hypotheses are retained.

## Finite-cluster error schedule

Fix n and all the band parameters through the nonzero loop in
DERIVATION.md. Let delta be the smallest positive angular root/pair margin.
It is bounded below by the minimum of half the two tree margins and
kappa eta/2. These are three different scales; the finite express delay
creates the third and cannot be omitted as an approximation error.

A strictly positive operator tolerance tau can be obtained directly from
P(z)=Tz^2+(M+F)z+T^T. Surround its 2n simple roots with disjoint small disks,
of radius less than delta/100, and bound |det P| below on their boundaries.
Also enclose the unit-circle complement and bound ||P|| above. Multilinearity
of the determinant and Rouche's theorem give a positive norm tolerance.
Hermitian reciprocal symmetry keeps the isolated perturbed roots on the
circle, and small enough disks preserve all pair comparisons. Every
constant here is finite and computable at the chosen n. No uniform-in-n
lower bound on tau is claimed.

The exact grammar removes all finite-K errors, including any zero-port
payload splitting. At h_L=0 the biased clock code is exactly degenerate
and its clock gap is at least 1/2. The charged carrier has a gap greater
than the conservative value used in the original nine-state conversion
tail estimate. Its two-cell source norm is sqrt(2 zeta), and the pair
conversion norm is at most 2|g|. The same elementary even Neumann series
therefore gives the safe bound

    ||remainder at order lambda^2|| <= lambda^2 zeta g^4/4,
    |g|<=1/4.

For example the conservative zero/one-particle coefficient extraction
factor 64n^2 and the two nearest-cell directions allow the sufficient
normalized quadratic condition

    16 n^2 D^(2n+2) g^2 < tau/3.                    (1)

The remaining nonquadratic operators have a comparable local bound; they
are not discarded in a many-body conclusion. Formula (1) only controls
the extracted quadratic symbol on a fixed cluster.

Next set h_L=lambda^2 g^2 zeta D^(-2n-2) epsilon/q^2. A fixed-cluster
analytic perturbation expansion, keeping this actual field in the
Hamiltonian, has a finite order-lambda^4 bound C_lambda lambda^4. This
includes the field-induced distortion of the vacuum clock and its source
visibility. Port parity removes odd powers. Choosing lambda after g to
obey

    C_lambda lambda^2 D^(2n+2)/(g^2 zeta) < tau/3     (2)

is sufficient for that finite-cluster quadratic correction. The remaining
budget may be assigned to numerical enclosures or other declared local
quadratic corrections. These are existence bounds, not measured budgets
for a many-cell calculation, and no old numerical port constant is imported
without rechecking its hypotheses for the biased clock.

## Number of terms and coefficient spread

The extra bias can be simplified using completeness before counting terms:

    -sum_(y not L) P_y tensor P_W
      =-I tensor P_W+(P_L0+P_L1) tensor P_W.

In the translational sum, its first term just shifts the W on-site
coefficient. There is no omitted physical interaction in this rewriting.
After this reduction and the choice K=0, a basis consisting of on-site
diagonal matrix units and off-site matrix units gives:

| Contribution | Literal count |
| --- | ---: |
| Vacuum head hops | 2 Hermitian transition pairs |
| Active data hops | 4 pairs |
| Express hops | 4 pairs |
| Left switches in all neighbour contexts | 14 pairs |
| Right closing links in all neighbour contexts | 14 pairs |
| Conversion | 1 pair |
| Ports | 2 pairs |
| Reduced bias | 2 diagonal two-site projectors |
| Nonzero on-site diagonal entries | 6 projectors |

Thus there are 41 undirected transitions, or 82 directed off-diagonal
matrix units, and eight diagonal monomials: **90 literal monomials** in
this stated mixed one-/two-site matrix-unit basis. Writing a sum over
neighbour contexts as one symbol does not turn fourteen transitions into
one operator monomial. If instead each complex transition is decomposed
into standard real-symmetric and imaginary-antisymmetric Hermitian matrix
units, the generic n>=2 count is 45 off-diagonal plus eight diagonal,
or **53 real Hermitian basis elements**. The two counts refer to different
explicit bases, not different physical Hamiltonians.

After fixing the common positive scale, Delta=3, the three fixed phases,
and the bias/switch/express unit strengths, the adjustable positive
magnitudes are

    q,w,D,nu_head,g,lambda,t_d,h_L: eight real controls.

The selector fixes nu_head from q,n, leaving seven freely selectable
continuous inputs to the existence schedule. The expressions q=(wu)^2
and h_L=W epsilon/q^2 merely parameterize those choices; they do not remove
literal operators. There is no extra freely tunable phase for each band.

For the real Hermitian basis above and generic n>=2, all nonzero absolute
coefficients belong to

    {1,q,w, sin(pi/(2n)),cos(pi/(2n)),
     t_d,g,lambda/sqrt(2),
     nu_head-1,nu_head+3,nu_head+3+D,h_L}.

The two reduced bias terms have coefficient +1. On-site W has
nu_head-1=q^2-(n+2)E_n>0 for this schedule; A,B have nu_head+3; XA,XB
have nu_head+3+D. No identity constant, telescoping bond term or overall
scale remains in these counts. The normalized spread is the maximum of
this set divided by its minimum. For n=1 the express sine/cosine list is
adjusted by omitting the zero cosine component.

The proof does **not** improve this spread. It uses small q,w,eta,g,lambda,
a potentially large D, and a very small h_L. The resulting small velocities
are also suppressed by zeta and D^(-2n-2), but these emergent suppressions
do not replace the microscopic coefficient ratios just listed. The biased
clock gives a robust auxiliary gap; its length-selector density gap and
logical nonresonance margins can still be very small. Compared with the
original nine-state family this is a proof-oriented enlargement of the
control choices, not a claim of a simpler or experimentally robust model.

## Remaining physical work

The theorem removes the arbitrary-n *band nonresonance* obstacle for this
modified nine-letter family. It leaves:

- nonzero-port selection of the intended cell pattern in the microscopic
  ground sector, with a common bulk subtraction;
- volume-independent elimination and a controlled interacting effective
  Hamiltonian;
- stability against all allowed many-body scattering processes, and
  identification of a unitary compact infrared CFT and its Virasoro c.

A vacancy would be a tenth physical state. No mobile dressing estimate is
proved here, and the L1-only readout is not the conserved-charge field
studied in the separate mobile assignment. The new zero-port selector
cannot be used as a finite-J/t mobile theorem.

The most decisive next analytic calculation is a volume-uniform local
elimination and interacting stability budget at these nonresonant bands.
For making this particular theorem reviewable, the next computational
deliverable is an exact interval implementation of the nested parameter
schedule at one n beyond the previous seven certificates, with a report of
the actual normalized coefficient spread; a root scan alone would not add
to the all-n proof.
