# Velocity, momentum and the remaining infrared theorem

## 1. Compare with the logical scale, not only the clock gap

Fix the finite n band parameters from the all-n construction. Write

    s_0=zeta D^(-2n-2), W=lambda^2 g^2 s_0,
    v_min,cell=W vhat_min>0,
    v_max,cell=W vhat_max,
    delta=min(delta_2,delta_4)>0.                    (1)

The vhat values are physical eigenvalue slopes of the normalized Bloch
matrix, with Euclidean-normalized physical flavour vectors. A whitening
congruence does not preserve these slopes without the corresponding
metric correction. Delta_4 includes every different unordered pair sum,
including the opposite-pair comparisons resolved by finite express delay.

In physical-site units, one cell has b=n+1 sites, so

    v_min,site=b W vhat_min, delta_site=delta/b.

The fixed frozen sector has period b; its momentum selection rules are
modulo 2pi in cell units, or 2pi/b in site units. A continuum interaction
coefficient also gains b when switching to site coordinates. Consequently
the dimensionless interaction/velocity ratios below are unchanged. The
original microscopic Hamiltonian still has one-site translation symmetry;
the translated frozen sectors are separate vacua.

Subtract each connected operator's vacuum scalar before discussing
interactions. This costs at most a factor two in the full operator bounds.
Using kappa=kappa_CAR from UNIFORM.md, define the conservative nontrivial
matching error

    I <= 2 e^(2kappa) zeta lambda^2 g^4/[16(1-g^2/4)]
         +2 C_4 lambda^4+2 C_5 |lambda|^5,           (2)

where C_4 is the explicit connected bound (11) and C_5 is the unprojected
remainder constant (9). The last term remains an auxiliary coupling;
it is included in the budget, not asserted to be a logical interaction
already eliminated exactly. Scalar subtraction here is an energy reference
within the selected sector. DEFECTS.md uses the exact dressed bulk instead
and does not throw away sector-dependent scalar dressing.

The first term of I/W is proportional to D^(2n+2)g^2, not merely g^2
compared with an order-one auxiliary gap. The fourth-order term contributes
lambda^2 C_4/(g^2 s_0), and the full remainder contributes
|lambda|^3 C_5/(g^2 s_0). These are the additional hierarchies.

## 2. What a multichannel weak-coupling check requires

Use chiral fields normalized so a free fermion is exp(i phi_a),
x(ell)=||ell||^2/2 at the free point, and

    H_G=(1/4pi) integral dx (partial_x phi)^T V_G (partial_x phi),
    J=diag(+1 for right movers, -1 for left movers).

Microscopic charge conservation is exact for this nine-state model and
gives sum_a ell_a=0 for a local vertex. It is not an imposed extra symmetry.
Currents and their products correspond to ell=0 and must be retained.
The normal-ordered interactions in (2) determine the actual forward
matrix E in V_G=V_0+E; no Gaussian metric is chosen independently of them.
A loose normalization bound is ||E||<=4 I within the finite-order chiral
projection (the coefficient one-norm includes all branch labels).

For a positive V_G the equal-time Gaussian covariance is

    G=V_G^(-1/2) |V_G^(1/2) J V_G^(1/2)| V_G^(-1/2),
    x(ell)=ell^T G ell/2.                            (3)

This follows by diagonalizing the real symmetric chiral matrix in absolute
value. It includes mixing and unequal velocities; it reduces to I for a
diagonal free velocity matrix. In particular forward scattering changes
both the velocities and the scaling dimensions rather than disappearing
from the calculation.

Put kappa_v=vhat_max/vhat_min and y=4 I/(W vhat_min). If y<=1/4,
V_G is positive. Elementary positive-square-root perturbation estimates
give the conservative bound

    ||G-I||<=8 kappa_v^2 y.                          (4)

For clarity, use ||sqrt(A)-sqrt(B)||<=||A-B||/
(sqrt(lambda_min(A))+sqrt(lambda_min(B))). Apply it first to V_G,V_0,
then to the squares of V_G^(1/2) J V_G^(1/2) and V_0 J. Their smallest
singular values are at least lambda_min(V_G) and v_min, respectively.
For ||E||/v_min<=1/4 the resulting bound is smaller than
6 kappa_v^(3/2)||E||/v_min, which implies (4).

If y<=1/(48 kappa_v^2), every integer vertex with ||ell||^2>=6 has
dimension at least 5/2 in this finite-order Gaussian theory. The only
nontrivial neutral vectors with smaller squared norm are the bilinear
ones (norm squared 2) and four-fermion ones (norm squared 4). Their
nonzero oscillatory momenta are protected by delta_2 and delta_4, except
the ell=0 forward terms already included in V_G. This accounts for the
entire integer lattice near the free point, including higher commensurate
vertices, rather than checking only one selected Cooper operator.

A conservative **finite-order admissibility window** is therefore

    y < theta,
    theta=min{1/(48 kappa_v^2), delta/100,
                    1/[100(1+log(1/min(delta,1)))]}.                  (5)

The second condition leaves room for momentum shifts; the third leaves
room for marginal running until the oscillatory cutoff 1/delta. Equation
(5) is a quantitative perturbative matching/positivity test, **not a proved
all-scale RG theorem**. Bounding the actual self-energy shifts and all
running vertices requires an infrared construction for the full model.
Nonforward quartics cannot simply be discarded at the matching scale,
and their cutoff condition must use dressed momenta in that construction.

## 3. A finite nonzero schedule does fit every prescribed positive budget

Let tau>0 be the earlier finite-band quadratic tolerance, and fix a
positive admissibility threshold theta (for example (5)). Put

    sigma=min(tau/4, vhat_min theta/4).

For g<=1/4, the first term in (2), divided by W, is bounded by
(8/63)e^(2kappa)D^(2n+2)g^2. Thus the simple sufficient choice

    0<g<=min(1/4, exp(-kappa)D^(-n-1)sqrt(sigma))    (6)

puts it below sigma/3. Fix chi=g^2 s_0 epsilon/q^2, then compute C_4,C_5
and r_0 from the formulas in UNIFORM.md. Choose

    0<lambda<=min{
       r_0/2,
       lambda_bd,
       [sigma g^2 s_0/(6 C_4)]^(1/2),
       [sigma g^2 s_0/(6 C_5)]^(1/3),
       (1-q^2)^2 q^(2n+4)/(10sqrt(2)) } .            (7)

All quantities are positive and finite at each fixed n; coefficients
remain fixed as L grows. The second and third terms each contribute at
most sigma/3 to I/W. The last condition gives exact nonzero-port frozen
sector selection. Set h_L=lambda^2 chi without changing any earlier band
parameter or nu_head.

This proves **algebraic feasibility of a volume-independent local smallness
schedule measured against velocities and momentum margins**. It proves
the finite-order forward Gaussian test, conditional on using its correctly
projected couplings. It does not prove that a finite-order residual can be
ignored to arbitrarily low energies. The further premise needed for a CFT
claim is a volume-independent exact low-energy integration/renormalization
theorem that turns these bounds into controlled dressed momenta, covariance
and irrelevant remainders. Such a theorem is not supplied here.

## 4. Numerical scale audit without another root scan

`checks.py` reuses the saved 70-digit n=2 roots from the previous contribution
at D=10000, q=1/16, w=1/2, epsilon=1e-6 and eta=1e-6. It computes slopes
in the physical flavour metric. At this validation point:

    vhat_min = 7.9790e-7,
    s_0 = 1.50218e-29,
    delta_2 = 0.0053241660,
    delta_4 = 2.9770755e-7.

Taking the deliberately conservative fixed couplings g=2^-100 and
lambda=2^-320, the logarithms base ten of the three bounds in (2),
divided by W, are approximately

    conversion -27.81, fourth order -93.48, full residual -42.02.

They are much smaller than vhat_min times the diagnostic threshold
delta_4/100. The normal-form radius and exact frozen selector conditions
also hold. This scale audit is numerical and does not produce a new
interval certificate for tau or an RG certificate. The point is used
only to expose the physical size of the sufficient bounds.

In site units the smallest velocity is approximately 10^-287.31 and
h_L approximately 10^-285.28 in the microscopic unit-hop energy scale.
The largest normalized microscopic coefficient divided by h_L is about
10^289.28. These are not moderate ratios. Much of this extremity comes
from deliberately loose full-Fock and exponential-locality bounds, but
the calculation does not demonstrate that they can be relaxed safely.

The microscopic operator count and eight adjustable real magnitudes are
unchanged from the prior contribution (53 real Hermitian basis elements
after its diagonal reduction and removal of the common scale). No new
control or physical state was introduced. The progress is a controlled
local estimate and a ground-sector theorem, not a simplicity improvement.

## 5. Strongest conclusion and next calculation

Established here: connected nonquadratic matching through fourth port
order, a complete finite-g conversion remainder, an all-n attractive
readout interaction, a uniform local normal-form remainder, and exact
nonzero-port frozen ground-sector selection on compatible rings.

Still open: an all-scale low-energy matching/RG theorem for the selected
interacting chain, including the residual auxiliary coupling. This is now
the decisive missing calculation. It must control the marginal forward
sector and the dressed oscillatory cutoffs, and yield a unitary compact
Virasoro theory. Neither c=n nor microscopic unbounded c is proved here;
no central charge was fitted or measured.
