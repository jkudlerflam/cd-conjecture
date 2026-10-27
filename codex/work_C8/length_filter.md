# A two-state length filter and occupation selection at two cell lengths

Status: proved. The selection theorem is for the invariant sector with
one L and one R endpoint tag in the correct order. It does not select
that sector against arbitrary tag patterns, and is not an all-ring
or conformal construction.

## 1. An initialized helix accepts one desired length

Fix k>=2 and set phi=2pi/(k+2). Write

    |+_theta>=(|0>+exp(i theta)|1>)/sqrt(2),
    |s_phi>=(|01>-exp(-i phi)|10>)/sqrt(2),
    h_phi=|s_phi><s_phi|,
    a=I-|+_0><+_0|,
    b=I-|+_{(k-1)phi}><+_{(k-1)phi}|.

On an r-site open chain use

    F_r=sum_{j=1}^{r-1}(h_phi)_{j,j+1}+a_1+b_r.           (1)

Every term is positive. The ground kernel is nonzero if and only if

    r = k mod (k+2),                                     (2)

and at those lengths it is one-dimensional, spanned by

    |g_r>=tensor_{j=0}^{r-1}|+_{j phi}>.                  (3)

In particular the only accepted length in 2<=r<=k+2 is r=k.

Proof. The product of on-site rotations diag(1,exp(i j phi)) maps
the usual spin-half ferromagnetic singlet projectors to h_phi on
each bond. The kernel of the untwisted sum is the symmetric subspace,
because the adjacent transpositions generate all permutations.
The first boundary condition pins its first spin to |+_0>. Symmetry
then pins every spin, so the simultaneous kernel of the bulk and
left boundary term is exactly (3). Its right boundary expectation is

    beta_r=<g_r|b_r|g_r>=sin^2((r-k)phi/2).               (4)

It vanishes exactly when (2) holds. Positivity excludes other zero
states. QED.

The phase and boundary projector are fixed for the selected k and do
not depend on the physical ring length. Using an irrational phi/(2pi)
instead would accept exactly r=k among all positive lengths. The
rational choice above has a useful explicit rejection bound on the
finite range needed here.

## 2. Explicit rejection energy

For r>=2, the bulk plus left boundary Hamiltonian has gap at least

    kappa_r=1/[r(2r-1)]                                  (5)

above its unique ground state (3). Here and below the bounds are
deliberately conservative.

Proof of (5). Untwist and rotate |+_0> to |0>. Let

    H_left=n_1+sum_{i=1}^{r-1} h_i,
    n_j=|1><1|_j,   h_i=(I-S_{i,i+1})/2.

For any normalized psi of energy E, move site 1 to site j by the
j-1 adjacent swaps. If U is that product, n_j=U n_1 U^dag and

    ||n_j psi|| <= ||n_1 psi||+||(U^dag-I)psi||
       <= sqrt(<n_1>)+2 sum_{i<j} sqrt(<h_i>)
       <= sqrt((4j-3) E).

The product-minus-identity telescopes with the remaining unitaries
on the left of each factor, so its norm is bounded by the stated
sum on the original vector. Also (I-S_i)=2h_i and h_i is a projector.
The last inequality is Cauchy-Schwarz. Since n_j is a projector,

    <sum_j n_j> <= sum_j(4j-3) E = r(2r-1) E.

On the orthogonal complement of |0>^r, sum_j n_j>=I. This proves
(5). Unitary conjugation gives the twisted statement. QED.

Adding the mismatched right boundary yields

    lambda_min(F_r) >=
      kappa_r beta_r / [(1+sqrt(kappa_r))^2+beta_r].        (6)

To prove this, let psi=c g_r+eta, eta perpendicular to g_r, and
E=<psi|F_r|psi>. Then ||eta||^2<=E/kappa_r. Because b_r is a
projector,

    sqrt(beta_r)|c| = ||b_r c g_r||
       <= ||b_r psi||+||b_r eta||
       <= sqrt(E)+sqrt(E/kappa_r).

Together with |c|^2>=1-E/kappa_r, rearrangement gives (6).

For all 2<=r<=k+2 with r!=k, define

    kappa_* = 1/[(k+2)(2k+3)],
    beta_* = sin^2(pi/(k+2)),
    epsilon_* = kappa_* beta_* /
                  [(1+sqrt(kappa_*))^2+beta_*] > 0.        (7)

Then F_r>=epsilon_* I. The minimum nonzero sine-squared on the
relevant integer phases is at least beta_*, and kappa_r>=kappa_*.
The right side of (6) increases with either kappa_r or beta_r, which
proves the uniform claim. In particular epsilon_* is of order k^-4.
This is a finite interval of rejected lengths, not a uniform
spectral gap over every r.

## 3. Selecting k occupied sites at physical lengths k+1 and k+2

Take any positive uniform payload bond h_pay on q states, with
positive one-site boundary terms a_pay,b_pay. Suppose its k-site
core is frustration free with complete code V. Tensor one filter
qubit onto every occupied payload site. In occupied coordinates use

    h_core = h_pay tensor I + J I tensor h_phi,
    a_core = a_pay tensor I + J I tensor a,
    b_core = b_pay tensor I + J I tensor b,                (8)

where the tensor ordering on a two-site bond is understood by grouping
payload and filter factors at each site. The occupied-site dimension
is 2q. At length k the complete code is V tensor |g_k>.

Apply the endpoint-tag construction of tagged_vacancy.md to these
operators. Its local dimension is now 6q+1. Fix physical length
N=k+1 or k+2, but allow EVERY occupied number r=2,...,N in the
union of ordered-tag sectors L B^(r-2) R. The same repeated bond is
used at both lengths. Choose

    J >= 2U/epsilon_*,      0<w<U/(k+1).                 (9)

Then the global ground space within that entire union is in the
r=k sector. It is the complete code

    chi_{N-k} tensor V tensor |g_k>,                    (10)

and every other r sector lies at least

    U-w(N-1)>0                                          (11)

above it. There is no vacancy chemical potential in this theorem.

Proof. First set w=0. A charge configuration contributes U per
occupied-vacant interface. If it has at least two interfaces, its
energy is at least 2U because all payload and filter terms are
positive. If it has exactly one interface, all occupied sites form
one contiguous interval touching a physical endpoint. Its filter
Hamiltonian is then the full initialized F_r. For r!=k, (7)-(9)
give filter energy at least 2U, in addition to its interface energy.
The only configuration with no interface in the present r>=2
sector is full occupation, r=N. This also has r!=k and its initialized
filter energy is at least 2U. Therefore every wrong-r sector has
energy at least 2U at w=0.

The hopping operator on each physical bond exchanges a vacancy and
an occupied site and has norm one before multiplication by w.
There are N-1 bonds, so switching on hopping lowers any sector's
ground energy by at most w(N-1). Every wrong sector thus has

    E_wrong >= 2U-w(N-1).                                (12)

At r=k the initialized payload and filter are frustration free.
The endpoint-tag theorem proves the exact code (10) and energy
E_target=e_{N-k}. An endpoint droplet of the N-k vacancies is a
normalized charge trial state with expectation U and zero hopping
expectation. Thus E_target<=U. Combining with (12) proves (11).
Condition (9) makes it positive for both physical lengths. QED.

The first and last payload maps retain their factors r_p, because
the physical observable can act as F tensor I_filter and the filter
state has norm one. The cut-shift transport coefficient remains
beta I on the payload codes. Neither operation changes the number
of occupied sites in a code cell. Thus (10) supplies the occupied-
number selection that was absent from the original vacancy theorem,
provided the ordered one-L/one-R sector is available.

## 4. What this resolves and what it does not

The earlier failure of a vacancy chemical potential did not establish
that occupied-number selection is impossible. The filter changes the
fully occupied sector by a positive frustration cost, and the
interface energy prevents wrong-length particles from avoiding that
cost cheaply by splitting into separate clusters. This is the
mechanism of (11).

The construction also explains why the boundary-free finite-witness
theorem does not bound the initialized length selector. The filter
uses boundary terms whose relative orientation depends on k. It
can accept an arbitrarily late first length at fixed filter
dimension two. The repeated tagged bond still has many other tag
sectors, including zero-cost product sectors. The finite-witness
theorem applies to the complete kernel and is not contradicted.

Three restrictions remain. First, the theorem assumes the one-L/one-R
ordered tag sector. Competing tag counts, invalid ordering, and
marker sectors in a full ring need separate control. Second, it
uses physical cell lengths k+1 and k+2, whose selection and mobile
coupling must be combined consistently with the payload. Third,
condition w<U/(k+1) places this charge Hamiltonian in a strongly
attractive regime if it is used on an unrestricted long ring.
The calculation therefore cannot be repurposed as a proof of a
homogeneous gapless charge fluid. Its use is inside finite cells.

There is still no uniform NN high-channel core, and no proof of
finite-coupling conformal behavior for the combined model.

## 5. Checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 length_filter_checks.py

The script checks k=3,4, every filter length 2<=r<=k+2, and all
occupied-number sectors 2<=r<=N at N=k+1,k+2. It uses U=1,w=.1
and J=2/epsilon_*, approximately 433.953 and 797.985. The filter
has only the specified zero state at r=k, its mismatched energies
obey (6)-(7), and the pinned-left gaps obey (5).

Every number-sector inequality (12) holds. The measured wrong-sector
separations for the four cells are about .99135, .89036, .99010,
and .86476, above the respective guarantees .7,.6,.6,.5. The largest
sector matrix has dimension 240. Ground-code residuals are below
4.4e-13 despite the large conservative filter coupling. The output
is length_filter_checks.json. These are exact finite-cell checks,
not central-charge measurements.
