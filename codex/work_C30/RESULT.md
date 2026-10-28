# C30. Symmetric switch-cell selection and mobile padding

[Proved] A nine-letter nearest-neighbor selector forms exactly
`P_L L^m M R^m P_R` with one switch and two ports. One chemical potential
selects any integer m>=1. It has an explicit positive energy margin
against every other complete cell, including unequal arms, and a
volume-independent bound against arbitrary illegal sectors. With the
carrier spaces of routeA_compile, its frozen local dimension is 120.

[Proved] The two-arm selector must not be placed at its adjacent-m
crossing and called a symmetric two-length fluid. Four arm pairs then
coexist, including both unequal-arm pairs. Also, adding an R site adds
a logical fermion, so that crossing does not preserve C8's payload.

[Proved, using the cited C8 finite-cell theorems] There is a positive
initialized symmetric core with local dimension 84 that can instead be
put inside C8's ordinary padding wrapper. It has k=2m+3 occupied sites
and m logical modes at both physical lengths k+1 and k+2. The complete
mobile selector has dimension 2028, or 8112 with the two endpoint
compensator helices. These are finite microscopic constructions and
alphabet counts. They are not a proof of c=m+1 for the carrier model.

[Assumed] C8's wrapper, occupation filter, and marker-transport theorems
are reused from `../work_C8/{length_filter,arrow_wrapper,selector_wrapper,
helix_equalizer}.md`. Their hypotheses are verified below for the new
core. C13/C20 supply conditional infrared criteria, not a finite-coupling
theorem for this cell. The source convention and carrier interface are
those of `../../contrib/routeA_compile/README.md`, Sections 1 and 4.

## 1. The repeated selector bond

[Proved construction] The following nine structural letters are arrows
on five auxiliary vertices. Arrows are labels, not additional sites.

| Letter | Arrow | Carrier dimension |
| :-- | :-- | --: |
| P_L | 0 to 1 | 4 |
| L_0 | 1 to 1 | 4 |
| L_X | 1 to 2 | 4 |
| L_1 | 2 to 2 | 4 |
| M | 2 to 3 | 4 |
| R_0 | 3 to 3 | 32 |
| R_X | 3 to 4 | 32 |
| R_1 | 4 to 4 | 32 |
| P_R | 4 to 0 | 4 |

[Proved construction] An ordered pair is legal when its two arrows
concatenate. Let Q be the projector onto all other pairs. With t>0,
put on every bond the same operator

    h_sel = K Q + (mu/2)([P_L] tensor I + I tensor [P_L])
            -t (|L_0 L_X><L_X L_1| + h.c.)
            -t (|R_0 R_X><R_X R_1| + h.c.),
    K=32t.                                                   (1)

The rewrites are identities on the carrier spaces at both physical
sites. The three L letters all have the same carrier space, as do the
three R letters. Thus (1) is a well-defined Hermitian two-site operator
on dimension 120. Before adding carriers the dimension is nine.

[Proved] Each rewrite preserves its external arrow vertices and is
legal on its internal bond. Consequently every individual Q_i commutes
with the Hamiltonian, including for overlapping bonds. The port and
switch positions are conserved as well. This is exact sector separation,
not an infinite-penalty limit.

[Proved] Every complete legal cell is

    P_L L_0^j L_X L_1^(a-1-j) M
        R_0^h R_X R_1^(c-1-h) P_R,
    a,c>=1, 0<=j<a, 0<=h<c.                                (2)

The only legal periodic words without P_L consist of one repeated loop
letter, L_0, L_1, R_0, or R_1. Indeed, every directed cycle that uses a
non-loop arrow must traverse all five vertices in order. This also
proves that missing ports, extra switches, and missing heads require
an illegal bond unless the entire ring is one of those four loop words.

[Proved] In a fixed cell the two heads are independent path walkers.
Their full spectrum and ground vectors are

    E_(p,q)(a,c)=mu-2t cos(pi p/(a+1))-2t cos(pi q/(c+1)),
    p=1,...,a, q=1,...,c,
    v_a(j)=sqrt(2/(a+1)) sin(pi(j+1)/(a+1)).                 (3)

This follows by substituting the sine vector in the tridiagonal hopping
matrix, including its zero boundary values. The sine vectors are an
orthonormal basis. The cell length is a+c+3, including both ports and M.
In particular the physical arm letters never fluctuate between L and R.
Only the neutral head position fluctuates inside an arm.

## 2. One chemical potential selects both arms

[Proved] Define

    f(r)=-2t cos(pi/(r+1)),     d_r=f(r+1)-f(r), r>=1.

The continuous extension of f is decreasing and strictly convex for
r>=1. With x=r+1>=2,

    f''(r)=4t pi sin(pi/x)/x^3 + 2t pi^2 cos(pi/x)/x^4 > 0.

Thus d_r<0 and d_(r+1)>d_r. For m>=2 choose

    e_m=(d_(m-1)+d_m)/2,
    mu_m=(2m+3)e_m-2f(m),
    eta_m=(d_m-d_(m-1))/2 > 0.                             (4)

For m=1 use e_1=-2t, mu_1=-10t, eta_1=t.

[Proved] For every r>=1, telescoping the increasing discrete slopes gives

    f(r)-f(m)-e_m(r-m) >= eta_m |r-m|.                      (5)

For r>m, each increment d_j-e_m is at least d_m-e_m. For r<m,
each reversed increment e_m-d_j is at least e_m-d_(m-1).
These are eta_m. For m=1 only the first case is needed and
d_1-e_1=t. Adding the two inequalities gives the explicit bound

    E_0(a,c)-e_m(a+c+3)
       >= eta_m (|a-m|+|c-m|)
       >= eta_m |a-c|.                                    (6)

Equality at zero occurs only at a=c=m. Every wrong complete cell costs
at least eta_m relative to the common energy-density reference.
For m>=2 the single-arm nearest competitors attain eta_m, so this
coefficient cannot be increased. All arbitrarily long cells are covered.

[Proved] This is genuinely selection by one chemical potential with a
fixed selector bond otherwise. Define

    C_m=(2m+3)d_m-2f(m).

Then C_1=-5t and

    C_(m+1)-C_m=(2m+5)(d_(m+1)-d_m)>0.

For m>=2, every mu in (C_(m-1),C_m) selects m; for m=1,
every mu<C_1 selects 1. To see this, put
`e=(mu+2f(m))/(2m+3)` in (5). The interval is precisely
`d_(m-1)<e<d_m`. The concrete choices (4) satisfy
`-10t<=mu_m<4t` and e_m<0. The upper bound follows from
`-2f(m)<4t` and the negative slope term; the lower bound for m>=2
follows from C_(m-1)>=-5t.

[Numerical, arithmetic from (4), t=1] The following values are also
checked against all arm pairs 1<=a,c<=48.

| m | mu_m | e_m | eta_m |
| --: | --: | --: | --: |
| 1 | -10.0000000000 | -2.0000000000 | 1.0000000000 |
| 2 | -2.9497474683 | -0.7071067812 | 0.2928932188 |
| 3 | 0.0472741754 | -0.3090169944 | 0.1051965680 |
| 4 | 1.4879631289 | -0.1589186226 | 0.0449018038 |
| 5 | 2.2687272593 | -0.0919518735 | 0.0220649453 |
| 8 | 3.2424077913 | -0.0271769838 | 0.0044491928 |

[Proved consistency checks] Taylor expansion gives
`eta_m=3 pi^2 t/(m+1)^4+O(t/m^5)`, `mu_m -> 4t` from below,
and `e_m -> 0` from below. All energies scale linearly with t.
No coupling or alphabet size depends on the total ring length.
The shrinking margin at large m is allowed, but must be included in
any later perturbation hierarchy.

## 3. All illegal sectors and the complete frozen code

[Proved] Let z>=1 illegal bonds occur on a ring of length L. Cut at
these conserved bonds. The hopping vanishes across every cut.
Each resulting legal open string consists of complete cells plus at
most two incomplete fragments. Each fragment contains at most two
heads and at most one P_L. Its kinetic energy is at least -4t,
and its chemical contribution is at least -mu_-, where
`mu_-=max(0,-mu_m)`. This statement also covers fragments with no heads
or no ports. A head path, even if truncated, has hopping energy >=-2t.

[Proved] Complete cells obey (6). Since e_m<0, a fragment of s sites
has energy at least `e_m s-4t-mu_-`. Therefore

    H_sel >= e_m L + (K-8t-2mu_-) z
           >= e_m L + 4t z.                             (7)

This uses only mu_-<=10t and K=32t. It covers arbitrarily many adjacent
violations and arbitrary malformed segments. The four legal no-port
rings have energy zero, strictly above e_m L.

[Proved] Write b=2m+3. On every ring with b dividing L, the selector
has exactly b ground states, the translations of the length-b tiling
with the unique two-head ground vector in each cell. Its gap is bounded
below, uniformly in those L, by

    Delta_sel >= min {eta_m, omega_m, 4t, |e_m| b},
    omega_m=2t[cos(pi/(m+1))-cos(2pi/(m+1))], m>=2.        (8)

For m=1 set omega_1=+infinity because each walker has only one state.
Every competing legal tiling has a wrong cell, every head excitation
costs at least omega_m, and (7) handles violations. These alternatives
exhaust the Hilbert space. This is a frozen construction with broken
translation, not the desired mobile phase on all rings.

[Proved construction] Tensor the compile carrier operators with the
identity on the arm head labels. They depend only on the underlying
letter P_L,L,M,R,P_R. Thus the selector and carrier Hamiltonians commute
inside every structural sector. The interface has n=c logical data
orbitals in a cell (a,c); the selected cell has n=m and code dimension
2^m. There is no first-order port compression, as in the compile code.

[Proved sufficient positivity condition] A convenient explicit carrier
bound is

    delta_car=e_min-3 w_max-2 tau_max>0,
    H_car >= delta_car N_car.                            (9)

Here N_car counts non-code dimers, e_min is the smallest carrier cost,
w_max bounds the magnitude of each carrier hop, and tau_max bounds
each conversion. The express/data graph has maximum degree three.
A carrier hop connects two configurations differing by one carrier
move. Add its absolute amplitude to both corresponding diagonal
entries. The resulting two-by-two matrix is positive, including for
complex hopping phases. The total diagonal charge per configuration
is at most 3 w_max per carrier. Conversions preserve N_car and require
at most 2 tau_max per carrier by the same two-by-two bound. The unused
onsite energy gives (9). Fermionic signs change phases, not this bound.

[Proved] This argument also supplies positive local terms, if needed.
Assign to each letter and lane a fixed incoming and outgoing diagonal
hop budget. For E(M) they are w_max and 2w_max; for ordinary internal
lanes they are at most w_max each. On each physical bond the assigned
budgets dominate all its hopping blocks. Missing graph edges leave
nonnegative unused budget. Assign 2 tau_max on each converting lane
to the onsite conversion blocks. The remaining onsite diagonal is at
least delta_car per carrier. All these operators annihilate the code.
Endpoint half-onsite terms can be restored by positive boundary terms
when the bond is used on an open chain.

[Proved] The kernel in (9) is all dimers in their code state with the
data orbitals free. This proves the complete carrier code and a gap
at least delta_car for arbitrary a,c, not just in the one-carrier
sector. At e=1,w_max=0.3, the sufficient condition is tau_max<0.05.
It includes the small-conversion switch regimes quoted in compile.
It is sufficient, not a claimed sharp bound for its larger-tau ED.

[Proved] Nonnegative carrier terms do not weaken (6) or (7), and vanish
on the selected code. Thus the complete frozen ground space has dimension
`b 2^(mL/b)` and a gap to its complement at least
`min(Delta_sel,delta_car)`. The microscopic dimension is

    d_frozen = 4+4+4 + 3*4 + 3*32 = 120.                (10)

All head labels are neutral. The compile fermion-number symmetry is
preserved exactly. Intercell coupling is set to zero in this exact
selection theorem; a finite intercell coupling is a separate low-energy
problem.

## 4. Why the adjacent-m crossing is insufficient

[Proved] At mu=C_m the common energy density is e=d_m, and

    mu+f(a)+f(c)=e(a+c+3),  a,c in {m,m+1}.              (11)

Thus `(m,m)`, `(m,m+1)`, `(m+1,m)`, and `(m+1,m+1)` all lie on the
same supporting line. The middle two have physical length 2m+4.
This is an exact obstruction to the independent-arm crossing scheme,
not a no-go for all local symmetry selectors. No energy rescaling of
these two independent identical walkers removes it while preserving
both symmetric endpoints of the crossing.

[Proved] Even imposing equality by a further constraint would leave
a different issue. The two symmetric codes would have dimensions
2^m and 2^(m+1). C8's scalar transport identity acts between two copies
of one fixed logical code. It cannot identify these two different
cell Hilbert spaces. A mechanism that adds a pinned inactive data
site or transports a logical mode would be an additional construction.

[Proved] Symmetry of the arms and symmetry of their hopping amplitudes
are separate conditions. The compile test described as a 5 percent
arm asymmetry multiplies one arm's hopping by 1.05. It does not change
its number of sites. That test gives no tolerance for unequal arm
lengths. Exact symmetric length selection exists here, so a geometric
asymmetry fallback is unnecessary. A small hopping asymmetry can still
be added on the carrier track without changing the selector proof.

## 5. What a literal length-step-two fluid would be

[Assumed model, not a microscopic construction] Suppose an independent
mechanism supplied fixed-payload cells with lengths b and b+2, b odd,
and scalar exchanges `(b,b+2)<->(b+2,b)` of amplitude -J. This section
solves that projected graph. The nine-letter selector at (11) does not
supply it, and a bare move of a shared port by two sites would spoil
the symmetry of both neighboring cells unless their centers also move.

[Proved for the assumed graph] With M cells and r long cells,

    L=bM+2r,    0<r<M,    M=L mod 2.                    (12)

Represent a state by a reference marker x in Z_L and a binary gap word
eta of weight r, with the identification

    (x,eta) ~ (x+b+2eta_0,S eta).                         (13)

The dimension is `(L/M) binomial(M,r)`. The identification is free since
no partial sum of positive gaps can equal L. Open-chain exchanges connect
all gap words at fixed x. A reference hop changes x by two. When L is
odd this connects all x. When L is even it connects each parity class,
and (13) identifies those two classes because b is odd. Hence the
physical graph is connected. There is no required twofold disconnected
parity sector in this odd-b case.

[Proved] Fourier transforming x with P=2pi ell/L and gauging the gap
wavefunction by `exp[-2iP sum_j j eta_j/M]` gives ordinary hard-core
hopping with twist 2P. Jordan-Wigner gives the complete spectrum

    q_a=(2pi/M)[I_a+(r-1)/2],
    E=kappa M-2J sum_(a=1)^r cos(q_a+2P/M),
    sum_a I_a+r(r-1)/2+ell=0 mod M,                     (14)

with distinct I_a in {0,...,M-1}. To check the last constraint, the
gauged wavefunction under S acquires `exp(-iPL/M)` from (13), whereas
its Slater determinant acquires `exp(i sum q_a)`. Equating them gives
(14). These changes of basis are unitary on the quotient, so the list
is complete, not just a set of eigenvectors. Perron-Frobenius on the
connected physical graph puts the ground state at P=0 and yields

    E0=kappa M-2J sin(pi r/M)/sin(pi/M).                  (15)

[Proved] With nu=r/M and rho=M/L=1/(b+2nu), the energy density and
its curvature are

    f(rho)=kappa rho-(2J/pi)rho sin(pi nu),
    f''(rho)=(pi J/2)rho^(-3) sin(pi nu)>0,
    v=2J rho^(-1) sin(pi nu),
    K_marker=4rho^2.                                   (16)

The last formula uses `f''=pi v/K_marker` with marker charge M.
At fixed L only changes of M by two are allowed. If the charge is
instead normalized as M/2, its stiffness parameter is rho^2. The
current spectrum has center displacement `2pi D/L`, where
`D=bj+2t_int` and D is any integer because b is odd. In the marker-charge
normalization its winding is D/2. These parity and winding conditions
must accompany K_marker; copying the step-one gluing would be wrong.

[Assumed continuum identification] Expanding (15) gives the Casimir
term `-pi v/(6L)`, and the particle-hole towers in (14) are those of
one compact boson. A unit of extra length is a pair of sites, but
there is one length mode, not two. Since gcd(b,b+2)=1, all sufficiently
large L admit such tilings. These facts describe the hypothetical
projected exchange graph, not its missing symmetric microscopic parent.

## 6. A positive symmetric core for C8

[Proved construction] To obtain a mobile fixed-payload model, keep
k=2m+3 occupied sites and use vacancies for the physical length
fluctuation. Replace the two head tracks inside this occupied core
by one neutral qubit on each L and R site. No qubit is needed on
P_L,M,P_R. Its alphabet dimension is

    q=4+4+4+2*4+2*32=84.                                (17)

Use the compile five-letter pattern grammar with forbidden-pair
penalties, and positive endpoint penalties requiring P_L at the first
occupied site and P_R at the last. A legal cell has at least one L and
one R because neither P_L M nor M P_R is allowed.

[Proved construction] Set p=2k+1, phi=2pi/p and define

    |+_theta>=(|0>+exp(i theta)|1>)/sqrt(2),
    h_phi=|s_phi><s_phi|,
    |s_phi>=(|01>-exp(-i phi)|10>)/sqrt(2),
    a=I-|+_0><+_0|,
    b_m=I-|+_((m-1)phi)><+_((m-1)phi)|.                 (18)

On L,L and R,R bonds put J_h h_phi. On P_L,L and M,R bonds pin
the right qubit by J_h a. On L,M and R,P_R bonds pin the left qubit
by J_h b_m. All are repeated two-site positive operators conditioned
on the letter pair. The carrier Hamiltonian is unchanged, with the
positive decomposition of (9). The helices are independent of it.

[Proved] An r-site arm has helix parent

    F_r=a_1+sum_(j=1)^(r-1) h_phi(j,j+1)+(b_m)_r.

Untwisting by the site rotations diag(1,exp(i j phi)) maps the bulk
kernel to the symmetric spin subspace. The left pin then forces all
spins to |+_0> in that frame, so the only possible zero vector is
`g_r=tensor_(j=0)^(r-1)|+_(j phi)>`. Its right-pin expectation is

    beta_r=sin^2((r-m)phi/2).                             (19)

Thus `ker F_r` has dimension one when r=m mod p and is zero otherwise.
Among 1<=r<=k there is exactly one accepted length, r=m. Positivity
is essential here: it identifies the kernel with the intersection
of the term kernels.

[Proved explicit rejection bound] The bulk plus left pin has gap
`kappa_r>=1/[r(2r-1)]`. For completeness, untwist and rotate the pin
to n_1=|1><1|. Moving that spin to site j by adjacent swaps gives

    <n_j> <= (4j-3)<H_left>,
    <sum_j n_j> <= r(2r-1)<H_left>.

The first inequality follows by telescoping the product of swaps,
using `I-S_i=2h_i`, and applying Cauchy-Schwarz to the left-pin
amplitude and the j-1 swap amplitudes. The number operator is at least
one on the orthogonal complement of the all-zero vector, proving the
gap. This also holds for r=1, with kappa_1=1.

[Proved] Write a normalized vector as c g_r+eta with eta perpendicular
to g_r and energy E under F_r. Then

    ||eta||^2<=E/kappa_r,
    sqrt(beta_r)|c|<=sqrt(E)+sqrt(E/kappa_r).

Eliminating |c| proves

    F_r >= kappa_r beta_r/[(1+sqrt(kappa_r))^2+beta_r] I
                                                        (20)

at rejected lengths. Uniformly for 1<=r<=k, r!=m, use

    kappa_*=1/[k(2k-1)],   beta_*=sin^2(pi/p),
    epsilon_*=kappa_* beta_*/[(1+sqrt(kappa_*))^2+beta_*]>0.
                                                        (21)

At r=m the gap is at least kappa_m because the matched right pin is
positive and annihilates the already unique ground vector.

[Proved] The complete k-site initialized core has only one legal
zero-energy structural word, P_L L^m M R^m P_R. A legal word with
fixed endpoint letters is a concatenation of complete cells. Every
zero-energy arm must have length m, so every such cell has length k.
Only one fits. Any other legal word has a rejected arm and costs at
least J_h epsilon_*. Illegal pairs and wrong endpoint letters cost
at least their positive penalty scales K_pat and K_end. Carrier
excitations cost at least delta_car. Therefore the complete core code
has dimension 2^m, with the explicit lower gap bound

    gamma_core >= min {K_pat,K_end,J_h epsilon_*,
                       J_h kappa_m,delta_car}>0.        (22)

Each code vector is a fixed carrier/data code tensored with the two
product helices. Neutral helices do not change the unwrapped port
source resolvent, because every carrier and port term acts as their
identity and keeps the structural word fixed.

[Proved scope] Unlike the frozen walker selector, (18) encodes m in
helix parameters. It is a positive initialized core for C8, not a
claim that the whole mobile family varies only one affine chemical
potential. The one-chemical-potential theorem is (1)-(6). C8 already
uses k-dependent occupation filters and endpoint compensators.

## 7. The actual mobile length fluid and its margins

[Proved application of the C8 theorems] Apply C8's occupied-site filter
and three-tag vacancy wrapper to the positive q=84 core. The occupied
sequence stays fixed at k=2m+3 while the physical cell has one or two
vacancies. Write

    B=k+1=2m+4,      physical lengths B and B+1.

The complete cell codes are `chi_p tensor V_m tensor g_filter`,
p=1,2, where V_m is the initialized code above. Both carry the same
m logical modes. A local H_0-to-H_2 vacancy conversion in the marker
hop transports the identity on that code. Symmetry is exact in the
ordered occupied coordinates; individual instantaneous vacancy
positions need not be symmetric in physical coordinates.

[Proved, explicit inherited margins] Here are parameters that make the
composition concrete. Choose U>0 and vacancy hopping
`0<w<U/(k+1)`. C8's occupation-filter qubit uses angle
`2pi/(k+2)` and strength `J_f>=2U/epsilon_occ`, where

    kappa_occ=1/[(k+2)(2k+3)],
    beta_occ=sin^2(pi/(k+2)),
    epsilon_occ=kappa_occ beta_occ/
                    [(1+sqrt(kappa_occ))^2+beta_occ].

Set the wrapper endpoint reward to 2U and its illegality penalty to
6U. Its complete selected ground energies are
`epsilon_B=-4U+e_1`, `epsilon_(B+1)=-4U+e_2`, with 0<=e_p<=U,
and its energy at every other open length is bounded below by -4U.
Here e_p is the smallest eigenvalue of the finite vacancy-interface
Hamiltonian at p vacancies. It is unambiguously defined by C8's local
bond, not an adjustable fitting parameter.

[Proved, explicit inherited margins] Define dimensionless quantities

    a_l=2(1-cos(pi/l)),
    delta_B=(B+1)a_B-B a_(B+1),
    e_c=a_(B+1)-a_B<0,   mu_c=2-delta_B>0,
    eta_B=min{(a_(B-1)-2a_B+a_(B+1))/2,
              (a_B-2a_(B+1)+a_(B+2))/3}>0.

Take the outer C7 selector hopping and chemical potential to be

    t_out=2U[1/eta_B+B/mu_c+1/|e_c|],
    mu_out=t_out mu_c+4U-(B+1)e_1+B e_2,
    e_out=t_out e_c+e_2-e_1,
    K_out=4t_out+9U.                                    (23)

Then e_out<0, the two selected lengths have energy e_out times their
length, every other complete length l has excess at least

    (t_out eta_B-U)(1+|l-B|)>0,                         (24)

and z outer-selector violations cost at least Uz above e_out L.
These are the bounds of C8 selector_wrapper, with its q-state input
now supplied by (17)-(22). They exclude arbitrary competing lengths
without an assumption on the core at those lengths.

[Proved, positive finite gaps] For each selected vacancy number p,
C8 gives a within-wrapper gap at least

    min{ delta_p p_min gamma /
          [(sqrt(delta_p)+sqrt(M0))^2+p_min gamma],
         U-w(N-1) },                                   (25)

where N=k+p, delta_p is the positive vacancy ground-state gap,
p_min is the minimum expectation of a neighboring-occupied-pair
adjacency projector, gamma is the initialized core-plus-filter gap, and
M0 is the sum of the norms of its bond and boundary terms. The vacancy
graph is finite and connected for w>0, so delta_p and p_min are strictly
positive. Equation (22) and the filter bound make gamma explicit.
Together with (24), U, the outer head gaps, and the marker-free energy
cost, (25) gives a volume-independent gap above the complete two-length
manifold. All sufficiently large L have a B/(B+1) tiling.

[Proved] The projected marker hopping is

    J_len=u A_B^2 beta_wrap beta_helix>0,
    A_B=sum_(j=1)^(B-1) sqrt(2/B)sin(pi j/B)
                           sqrt(2/(B+1))sin(pi j/(B+1)).

The vacancy and helix overlaps are strictly positive at the C8
parameters; without the optional equalizer helices set beta_helix=1.
The projected length graph is therefore the step-one C8 graph,
`L=B M+r`, rather than (12). It has one length mode under the usual
free-fermion continuum identification. Padding resolves the change
of logical dimension and the length-step-two issue at the same time.

[Proved alphabet count] The count, before any optimization, is

    initialized symmetric core                  q = 84
    occupied length-filter qubit                2q = 168
    three occupied tags plus three vacancies    6q+3 = 507
    outer four-letter C7 selector               24q+12 = 2028
    two physical endpoint-compensator qubits    96q+48 = 8112.   (26)

The last line is the mobile d0 used here for source comparison.
The 2028-dimensional version already supplies the unperturbed mobile
selector and scalar payload transport. Neither count includes an
unbuilt computation alphabet. All added tracks are neutral, and
vacancy transport preserves the charge carried by the occupied labels.
Exact fermion number survives the construction.

## 8. Source spectra, C25, and the infrared qualification

[Proved] C8's equality of compressed endpoint operators is not equality
of complete second-order source resolvents. A vacancy can remove an
occupied bond while a carrier propagates. This changes the effective
carrier Hamiltonian seen at p=1 and p=2. Thus exact occupied-word
symmetry does not imply exact left/right or length-independent source
spectra at finite vacancy hopping. C25 exhibits this problem in its
C18 test; the same distinction is necessary here.

[Proved finite-dimensional ordered-limit statement] The source mismatch
is tunable for this positive carrier core without assuming commuting
carrier dynamics. Fix w/U>0 first and scale all core terms by s>0.
Let P project onto the fast vacancy and auxiliary-helix ground spaces,
retaining the carrier Hilbert space in the fixed structural sector.
Write the source Hamiltonian as A+sB with
`A P=0`, `A>=g(1-P)`, `B>=0`, and `||B||=b0`.
Here g includes the small vacancy tunneling gap. On the carrier-excited
source sector let `C=PBP>=gamma P`, gamma>0. The positivity and
completeness proof above, together with positive adjacency weights,
ensure gamma>0 for every fixed w>0.

[Proved] The Schur complement on P is

    sC-s^2 PBQ (QAQ+sQBQ)^(-1) QBP,   Q=1-P.

If `s b0^2/g<=gamma/2`, its difference from sC has norm at most
`s^2 b0^2/g`, and the inverse identity gives

    || s P(A+sB)^(-1)P-C^(-1) ||
                  <= 2s b0^2/(g gamma^2).              (27)

For general bounded source vectors, the inverse block formulas show
that the scaled cross blocks and scaled Q,Q block also converge to
zero. For example their bounds are respectively
`2s b0/(g gamma)` and
`s/g+2s^2 b0^2/(g^2 gamma)`. These prove convergence of the entire
second-order bond, including compensator-created auxiliary excitations,
not only a ground-to-ground matrix element. The exact zero-carrier
code is removed before taking the source inverse.

[Proved] As w/U tends to zero at fixed finite k, the vacancy ground
vector tends to the reflection-even sum of the two endpoint droplets.
The charge diagonal assigns U to an endpoint droplet and at least 2U
to every other configuration. The connected hopping ground vector is
unique, positive, and reflection even. Projection onto the other charge
configurations has norm O(w/U), by its gap U and a finite matrix norm
bound. Both endpoint droplets leave every occupied core bond adjacent.
Since the adjacency indicators are diagonal, every averaged internal
bond coefficient is `1-O((w/U)^2)`. Hence

    C=H_source,bare+O((w/U)^2)

at fixed k. Boundary terms attached to occupied endpoint tags remain
present. C8's endpoint compensators have expectation 1/r_p, canceling
the projected source amplitude r_p for each p. The inverse identity
then bounds the difference of the corresponding resolvents by the
Hamiltonian difference divided by their positive source gaps.

[Proved] Consequently the compensated bonds obey the ordered limit

    s K2_(p,q)=K2_bare+O((w/U)^2)+O(s/g),  p,q=1,2,       (28)

with finite constants depending on k, source gaps, and physical source
norms. First choose w/U small and nonzero, then choose s relative to
the resulting g. Equation (27), not a uniform limit in both ratios,
justifies the second step. This extends the finite-matrix mechanism
behind C25 to the noncommuting carrier sources. No numerical constants
for the full wrapped switch resolvent have been computed here.

[Proved] The mixed coefficient of order u lambda vanishes for this
carrier implementation. Let N_car now count non-code dimers across
all cells. The unperturbed core and wrapper preserve it. The marker
hop transports occupied labels and moves only a vacancy across a cut,
so it also preserves N_car. A compile junction acting on the code
creates two carriers, either H,D or T,T, and has no zero-carrier
component. Therefore neither `P V_move R V_link P` nor its adjoint can
return to the zero-carrier code. This selection rule is specific to
this number-of-carriers-preserving construction and does not establish
vanishing for a general future port gadget.

[Proved limitation] Scalar, length-dependent second-order terms are
part of the length Hamiltonian. They cannot be discarded as a global
energy shift. Their mismatch is tunable in (28), but their size must
be compared with J_len. Likewise any traceless mismatch must be small
relative to the weakest retained payload scale, not just the largest
carrier hopping or the largest logical eigenvalue.

[Assumed infrared inputs] Applying C13/C20 requires local elimination
of the massive sectors at fixed nonzero couplings, a positive full
stiffness matrix, actual payload and length momenta avoiding every
relevant or marginal non-Gaussian resonance, and the corresponding
translation and winding-sector gluing. Exact total fermion number
forbids net pairing but does not forbid interflavor pair transfer.
The number of retained gapless payload bands must first be established.
The small logical eigenvalues in compile make the required hierarchy
strongly m dependent.

[Numerical input from compile, not repeated here] Its symmetric bare
switch gives the Gaussian band count c=n for n=1,...,4, and its
hopping-asymmetric phase-bearing n=2 example displaces the seas.
Neither result proves the same band count and nonresonance margins
for all m. C30 supplies selection and a finite mobile-core construction,
not that missing all-m spectral or infrared analysis.

## 9. Exact-diagonalization checks and reproducibility

[Numerical] `checks.py` uses only numpy and scipy. It imports no other
project code and writes only `checks_results.json` in this directory.
Run with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py

[Numerical, sizes and estimators] The completed checks are:

* All 162 surrounding-letter contexts for the two hopping rewrites
  preserve all three overlapping legality indicators.
* Complete two-walker spectra for every 1<=a,c<=8, dimensions ac<=64,
  agree with (3) to 6.22e-15 in maximum eigenvalue error. Selection
  inequalities were also evaluated for m=1,...,12 and a,c=1,...,48.
* The entire nine-letter ring spaces at L=3,4,5, dimensions 729,6561,
  59049, were decomposed into exact hopping components and diagonalized.
  Every illegal component obeyed (7). At m=1 and L=5 the ground energy
  is -10 and the ground multiplicity is five. At L=3,4 the four loop
  words have energy zero and are the ground states, as expected before
  a complete length-five cell fits.
* Arm-helix Hamiltonians at m=1,...,4 and every
  1<=r<=min(2m+3,9) were diagonalized, 30 matrices of dimensions up to
  512. Only r=m had a zero eigenvalue. All rejection bounds (20) and
  matched-length gap bounds passed with tolerance 1e-10.
* All legal fixed-endpoint pattern sectors at k=5,7,9,11 were checked
  by combining the arm spectra. There are 1,3,5,11 patterns, including
  multiple-cell patterns when they fit. Only the intended symmetric
  pattern has zero energy in each case.
* For the separately assumed step-two exchange graph, full spectra at
  (b,M,r)=(5,4,2),(5,5,2),(5,6,3), dimensions 36,58,120, agree with (14)
  to 7.55e-15. Each physical graph is connected. This checks the graph
  calculation only, not a microscopic symmetric pair-transfer rule.

[Proved scope of verification] The finite proofs cover arbitrary arm
lengths and illegal sectors. The ED checks small instances rather than
replacing those proofs. No many-cell carrier ED, wrapped-switch source
ED, central-charge fit, or thermodynamic simulation was performed.
No network or cluster job was used. No file outside this directory was
modified. The project-level STATE file is outside the permitted write
scope; `STATE.md` here contains the handoff.

[Proved conclusion] Symmetric frozen selection is possible with one
chemical potential at d0=120. A mobile fixed-payload realization uses
symmetric occupied cores and ordinary one-site padding, at d0=8112
including endpoint compensators. A literal length-step-two fluid has
one boson with modified charge and winding gluing if its exchange
parent can be supplied; that parent is not the independent-walker
crossing. [Assumed remaining step] The route-A conformal conclusion
still needs the all-m carrier band and resonance analysis and controlled
finite-coupling infrared transfer.
