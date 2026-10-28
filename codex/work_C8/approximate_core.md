# Vacancy lifting of a positive approximate code

Status: proved finite-cell spectral and observable bounds. This relaxes
exact frustration freeness in the finite-cell input. It does not supply
a uniform simulator or a thermodynamic infrared theorem.

## 1. Hypotheses and spectrum

Let h0,a,b be positive, and let

    H_core=sum_(j=1)^(k-1) h0_j+a_1+b_k.

Assume its lowest D eigenvalues lie in [0,epsilon], with isometry V
onto that spectral subspace, and every other eigenvalue is at least
gamma>epsilon. The terms need not annihilate V. Use the ordered-tag
vacancy construction of tagged_vacancy.md at N=k+p, with charge
ground vector chi, energy e_ch, and gap delta>0. Then

    H_ext-e_ch=A+B,
    A=(H_ch-e_ch) tensor I>=0,
    B=sum_j D_j tensor h0_j+I tensor(a_1+b_k)>=0.          (1)

Put

    p_min=min_j <chi|D_j|chi>>0,
    M=(k-1)||h0||+||a||+||b||,
    g0=delta p_min gamma /
           [(sqrt(delta)+sqrt(M))^2+p_min gamma].        (2)

If epsilon<g0, the lowest D eigenvalues of H_ext-e_ch lie in
[0,epsilon], and the (D+1)-st eigenvalue is at least g0. Thus the
low code band is isolated by at least g0-epsilon.

Proof. Let W0=chi tensor V and P=W0 W0^dag. Charge compression of
B is sum_j p_j h0_j+a_1+b_k, which lies between p_min H_core and
H_core as a quadratic form. Therefore P(A+B)P<=epsilon P.
The min-max principle gives lambda_D<=epsilon; positivity gives
lambda_1>=0.

For a normalized psi perpendicular to range(W0), decompose
psi=x+eta with x=chi tensor phi and eta perpendicular to the charge
ground state. Then phi is perpendicular to V. If E=<psi|(A+B)|psi>,

    ||eta||^2<=E/delta,
    <x|B|x>>=p_min gamma ||x||^2,
    sqrt(<x|B|x>)<=sqrt(E)+sqrt(M)||eta||.

These are the same elementary inequalities as the exact-core gap
proof, but only H_core>=gamma on V^perp was used. Eliminating eta
gives E>=g0. Hence Q(A+B)Q>=g0 Q for Q=I-P. Applying min-max on
this codimension-D space gives lambda_(D+1)>=g0. QED.

## 2. Code angle and boundary errors

Let P_low be the actual lowest-D spectral projector. Then

    ||P_low-P||<=min{1,sqrt(M epsilon)/(g0-epsilon)}=:s.  (3)

To prove it, note AP=0 and 0<=B<=M, so

    ||Q(A+B)P||^2<=||BP||^2
                     =||PB^2P||<=M||PBP||<=M epsilon.

For an orthonormal low-band eigenframe U with eigenvalue matrix E_U,
the Q component satisfies the Sylvester equation

    Q(A+B)Q (QU)-(QU)E_U=-Q(A+B)P(PU).

The first operator has spectrum >=g0 and the second <=epsilon.
Its solution can be written as an integral of
exp[-t Q(A+B)Q] times the right-hand side times exp(t E_U),
giving ||QU||<=sqrt(M epsilon)/(g0-epsilon). For equal-rank
projectors ||QU|| is their largest principal-angle sine, proving (3).

Align the logical bases by the polar decomposition of U^dag W0,
and call the resulting frame W. Principal angles give

    ||W-W0||<=sqrt(2)s.

For every physical operator F,

    ||W^dag F W-W0^dag F W0||<=2sqrt(2)s ||F||.           (4)

This applies to both endpoint maps. If two neighboring cells change
length and the exact reference codes have projected transport beta I
under a local operator C, the approximate codes obey

    ||transport_actual-beta I||
                   <=2sqrt(2)(s_1+s_2)||C||.             (5)

There are four isometry differences, two before and two after the
cut move; the product-isometry triangle inequality proves (5).
The right side may be large, so a useful approximation requires
epsilon small compared with g0^2/M, not merely below the core gap.

## 3. Adding the filter and complete finite-cell sectors

Tensor the occupied-coordinate helix filter onto the approximate core.
Its desired length-k ground vector is unique. The input low band
becomes V tensor g_k, with the same epsilon and with complement
threshold at least min(gamma,J_filter gamma_filter).

The filter's rejection estimates at wrong occupied numbers use only
positivity of the payload terms, so they still hold. The arrow-wrapper
endpoint and illegal-sector bounds likewise use positivity, not exact
annihilation. Thus on each physical length N=k+1,k+2 there are exactly
D eigenvalues in

    [e_p-4U, e_p-4U+epsilon],

and every other eigenvalue is at least

    e_p-4U+min{g0,U-v(N-1)}.                             (6)

Here the g0 for each length uses the combined payload/filter core.
Equations (3)-(5) control the two-length code and its port/transport
errors. The helix equalizer can be applied to the reference maps;
actual errors are then multiplied by the norms of its compensators.

This makes the finite-cell interface compatible with positive
approximate history simulations, PROVIDED such a simulation actually
preserves both endpoint observables to the required accuracy.
The lemma is not itself a construction of such a simulation.

## 4. Why this is not a volume-uniform low-band theorem

On an array of M_cells the retained band can have width M_cells epsilon,
while one auxiliary excitation costs a fixed amount. Its full spectrum
can therefore overlap the excluded bands as M_cells grows. The exact
base ground-space theorem of selector_wrapper.md cannot be quoted
unchanged when epsilon>0.

One can use the exact low-band projectors of each finite cell to define
a reference Hamiltonian that flattens those bands. The difference is
a local-in-cell perturbation of norm at most epsilon. This is a useful
setup for a local effective-theory argument, but the flattening term
need not be a microscopic nearest-neighbor bond. It is an analysis
device, not a proposed physical modification. Infrared control at
fixed error remains an independent obligation.

## 5. Checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 approximate_core_checks.py

The core is the k=3 graded qubit-swap parent plus tau times the total
occupation. It has the vacuum at energy zero and its symmetric
one-particle code state at energy tau, with the complement beginning
at 0.5+tau. The local decomposition uses positive bonds and positive
endpoint terms. It is not frustration free on the complete D=2 code.

At p=1,2, U=1, v=0.6 and tau=1e-6,1e-7,1e-8, all six complete
ordered-tag matrices, dimensions 32 and 80, obey (2)-(4). Their
proved complement thresholds are about 0.02403 and 0.004927, while
their actual third eigenvalues are 0.23246 and 0.033488. The largest
code angle is 2.5e-7 and the largest endpoint error is 1.2e-8. The
bounds are deliberately conservative. At tau=1e-6 the larger angle
bound is 0.2871, well above the observed angle. The charge hopping
used here tests the finite-cell lemma, not the separate wrapper
number-selection inequality v<U/(k+1).

Output is approximate_core_checks.json. No full periodic low-band
or conformal conclusion is drawn from these checks.
