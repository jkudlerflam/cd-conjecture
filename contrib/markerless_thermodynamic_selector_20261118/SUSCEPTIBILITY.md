# Exact half-occupancy contact response: what is and is not bounded

2026-11-18. Author/AI assistance: OpenAI Codex, GPT-6 Astra / Max assignment.
The proofs below await independent review. This concerns the ORIGINAL
corrected d=10 mobile model, not the modified biased all-n clock.

## 1. Scope and normalization

Keep the ten physical states, Hermiticity, microscopic nearest-neighbour
range, one-site translation symmetry, and coefficients independent of L.
Initially lambda=f=0. In the legal all-vacuum-head block, the exact operator
is the previous contribution's

    H_w=H_ch+J sum_(c,s) A_(a_c+s) h_(b_c,s),
    H_ch=-t sum_j(S_(j-1) S_j^*+h.c.),
    A_j=|g_j=0><g_j=0|,  sum_j g_j=V.                  (1)

There are N material gaps. Here V=N and physical L=2N; hence rho=1/2.
The chemical potential is a common scalar at fixed N and is removed.
The susceptibility is defined with no factor of two:

    G_ij=<x_i,(H_ch-E_ch)^(-1)_perp x_j>,
    x_j=(A_j-alpha_N)Omega,
    gamma_N(q)=sum_j exp(iqj)G_j0.                       (2)

The missing bound from the preceding round is sup_(N,q) t gamma_N(q)<infty.
**This full bound is neither proved nor refuted here.** We prove a uniform
bound for q=0, calculate its thermodynamic value exactly, give an exact
finite-N determinant expression for every q, and isolate why a tempting
free-fermion two-particle-hole proof does not extend to q!=0. RETARDATION.md
then treats a nonzero-J dynamical kernel without inventing a bath gap.

## 2. An exact spectral formula at every material momentum

The gap problem is the phase model. Its Bethe boundary condition is
e^(i L p_a)=(-1)^(V-1)e^(iP), P=sum_a p_a mod 2pi. This agrees with
Bogoliubov, Izergin and Kitanine, [equations (5)-(7)](https://arxiv.org/pdf/solv-int/9612002),
after removing their number-dependent constant and rescaling the hopping
to t. The following contact-projector determinant is derived directly.

For even N=V (in particular every compatible b=10 ring), put L=2N and
let I range over the L half-integers

    -L/2+1/2, ..., L/2-1/2.

For l=0,...,N-1 and each N-element subset I with sum I=l modulo L, set

    q=2pi l/N,
    p_a=(2pi I_a+q)/L,
    E(p)=-2t sum_a cos p_a.                              (3)

Represent a gap configuration by the ordered boson coordinates
0<=x_1<=...<=x_N<=N-1 and define y_a=x_a+a-1. These are distinct
integers in {0,...,L-2}. An unnormalized eigenvector is

    Psi_p(x)=det_(a,b) exp(i p_a y_b).                    (4)

Interior moves are ordinary fermionic hops in the y coordinates. At a
cyclic move the reordering sign and the common translation give precisely
the Bethe boundary condition above; thus (4) has energy (3). Cauchy-Binet
gives its norm:

    ||Psi_p||^2
      =det[sum_(y=0)^(L-2) exp(i(p_b-p_a)y)]
      =det[L 1-u u^*]=N L^(N-1).                       (5)

Within one q sector distinct sets are orthogonal: different sets with
the same total momentum cannot differ in just one occupied momentum,
and a rank-one deletion of the complete L-site overlap cannot connect
sets differing in two or more momenta. Different q sectors are
orthogonal by translation. Shifting every I by one modulo L changes
sum I by N modulo L and pairs residues l and l+N. Therefore the number
of states retained in (3) is binomial(2N,N)/2=binomial(2N-1,N), exactly
the gap-space dimension. This proves completeness in the case used here.

Let p^0 be the symmetric filled ground sea. The contact A_0 requires
x_1>=1, which removes y=0 as well as the already missing y=L-1.
Consequently the normalized contact form factor is exactly

    F(p,p^0)=det M/[N L^(N-1)],
    M_ab=sum_(y=1)^(L-2) exp(i(p_b-p_a^0)y).              (6)

In terms of these finite determinants,

    gamma_N(q)=N sum_(p in sector q, p!=p^0)
                   |F(p,p^0)|^2/[E(p)-E(p^0)].          (7)

The factor N is from the normalized Fourier source N^(-1/2)sum_j
e^(-iqj)A_j. Equation (7) is an exact positive spectral sum, not a
free density Lindhard function or a truncation in particle-hole number.
Its equality to the independent gap-coordinate resolvent is checked on
all states for N=4,6, to 2.8e-14 or better.

For q=0, both seas use the SAME momentum grid. M is the complete
fermion overlap with two rank-one deletions, so only at most two
particle-hole changes survive. For q!=0 the whole grid shifts by q/L;
the complete overlap is a Cauchy matrix, and this rank-two shortcut fails.
For example, at N=6,q=pi, sets differing in more than two Bethe integers
contribute 0.0004500937213/t to gamma. This classification of integer
sets is a grid convention, not a new physical quasiparticle count.
It suffices to disprove a proposed same-grid two-particle-hole truncation.
It does not disprove the desired susceptibility bound.

## 3. A volume-uniform theorem for the exact q=0 channel

The uniform contact sum is the physical adjacent-occupied-pair operator

    W=sum_j A_j=sum_(x=0)^(L-1) n_x n_(x+1).              (8)

At q=0 the root gluing introduces no extra twist. Thus its response is
that of N=L/2 free spinless fermions. With k on the antiperiodic grid,
the only excited matrix elements are two-particle-hole states with
p_1+p_2=h_1+h_2 modulo 2pi. For unordered hole and particle pairs,

    <p_1 p_2;h_1 h_2|W|FS>
      = (2/L)[cos(p_1-h_1)-cos(p_1-h_2)]                 (9)

up to a basis phase, and

    gamma_N(0)=(1/N)sum |matrix element|^2/
               [epsilon_p1+epsilon_p2-epsilon_h1-epsilon_h2]. (10)

**[Proved] For every even N>=2, 0<=gamma_N(0)<=3/(4t).**

Proof: write L=4M. On each side of either Fermi point, distances are
(2r+1)pi/L, r=0,...,M-1. Each positive particle or hole energy is at
least (8t/L)(r+1/2), by sin x>=2x/pi for x in [0,pi/2]. The matrix
element in (9) has square at most 16/L^2. Replace the unordered sum by
one quarter of the ordered sum over h_1,h_2,p_1, with p_2 fixed by
momentum; discard its occupancy restriction and its positive energy.
There are at most eight choices of Fermi-side labels. Hence

    gamma_N(0)<=8/(tL^2) sum_(r,s,u=0)^(M-1)
                                      1/(r+s+u+3/2).

Grouping by k=max(r,s,u) bounds the sum by

    sum_(k=0)^(M-1) (3k^2+3k+1)/(k+3/2)
       <=sum_k 3(k+1/2)=(3/2)M^2.

This proves the stated 3/(4t), with no volume-independent charge gap.
It is conservative; the exact limiting value is much smaller.

## 4. Thermodynamic value and the commensurate low-energy branch

**[Derived from the exact XXZ energy; checked against (10)]**

    lim_(N->infinity) t gamma_N(0)
                     =(24-pi^2)/(18pi^3)
                     =0.02531816334516650....            (11)

Here is the normalization and differentiation. The scalar reference
H_ch+uW maps to the half-filled spinless t-V chain. If e_X(Delta) is
the energy per site of sum(SxSx+SySy+Delta SzSz), its energy per MATERIAL
site is

    e_gap(u)=4t e_X(u/(2t))+u/2.

The input energy formula is [Kato et al., equation (11)](https://arxiv.org/pdf/cond-mat/0304475):

    e_X=cos gamma/4-sin gamma I(gamma),
    I(gamma)=integral_0^infinity [1-coth(pi x)tanh(gamma x)]dx,
    Delta=cos gamma.

At Delta=0, I=1/pi and

    I''(pi/2)=8/pi^3 integral_0^infinity
                      y^2[2 sech^2 y-sech^4 y]dy
               =8/(9pi)+8/(3pi^3),
    e_X''(0)=(pi^2-24)/(9pi^3).

The integrals use integral y^2 sech^2 y=pi^2/12 and
integral y^2 sech^4 y=(pi^2-6)/18. Since gamma(0)=-e_gap''(0)/2,
(11) follows. The finite sums (10) are dominated Riemann sums: the
three-dimensional 1/(r+s+u) corner singularity is integrable, as the
bound in §3 explicitly shows. Thus this value is not an extrapolation
of a few small rings.

The first commensurate current branch is also explicit. Moving the
two leftmost occupied momenta to their translates by pi gives, for
one of the two current directions,

    contribution to gamma_N(0)
      =8 sin^4(pi/L)/[t L^3(sin(pi/L)+sin(3pi/L))]
      ~pi^3/(32t N^6).                                 (12)

The branch has energy 4t[sin(pi/L)+sin(3pi/L)], which tends to zero.
Its susceptibility weight nevertheless vanishes as N^-6, due to the
antisymmetric two-fermion matrix element. It is therefore not a
logarithmic Peierls divergence. This is an exact test of the first
commensurate branch, not a bound on every q near zero.

## 5. Long wavelength, canonical zero mode, and actual patch sources

For vacancy density per gap r=V/N, rho=1/(1+r), the free equation of
state and contact expectation give

    e(r)=-(2t/pi)(1+r)sin[pi/(1+r)],
    alpha(r)=rho-sin^2(pi rho)/(pi^2 rho),
    e''(1)=pi t/4, alpha'(1)=-(1/4+1/pi^2).

The hydrodynamic density-relaxation coefficient is therefore

    C=[alpha'(1)]^2/[2e''(1)]
     =0.07857579749228010.../t.                          (13)

**[Hydrodynamic inference, not a uniform lattice theorem]** If the
remaining irreducible static response is regular at q=0, the small
nonzero-q limit is gamma(0)+C=0.1038939608374466.../t, whereas the
exact canonical q=0 value is (11). The smooth finite-ring kernel is
C(delta_ij-1/N). This zero-mode distinction was already identified
in `route_a/DEFECT_TRANSFER.md`; (11) supplies the new microscopic
canonical contribution. Neither value is by itself an upper bound
at arbitrary momentum.

The bounded gap-resolvent calculations give:

| N=V | gap-space size | t gamma_N(0) | max_q t gamma_N(q) |
| ---: | ---: | ---: | ---: |
| 4 | 35 | 0.02453690 | 0.09627416 |
| 6 | 462 | 0.02497254 | 0.09977670 |
| 8 | 6435 | 0.02512405 | 0.10137896 |
| 10 | 92378 | 0.02519402 | 0.10220945 |

These support the predicted finite long-wavelength limit; they do not
certify a supremum over N or establish a microscopic interacting phase.

For the original n=9 sources, let p be the period-10 profile of the
unchanged target and d=sigma_B-sigma_A the UNWRAPPED local replacement,
including its two distinct exterior contacts. Put e_10=epsilon_10/10
and Delta=sum_j d_j. The exact q=0 part of the bulk-subtracted second
coefficient, including the cross term, is

    D_2,zero=-gamma_N(0)[2 e_10 Delta+Delta^2/N].         (14)

The clock transition sources sum to zero, so they contribute nothing
to (14). At fixed patch and N->infinity this is favorable, and equals
3.0171370245e-6/t for the balanced replacement and 1.5085685123e-5/t for either
head-number patch; full digits are in checks.json. It is a finite
defect term because the unchanged background has mean e_10. Dropping
the cross term would incorrectly remove its thermodynamic limit.

In the acoustic approximation, restrict to |q|<Lambda<2pi/10. The
nonzero Fourier harmonics of p then do not contribute to the cross
term. Since d(q)=Delta+O(q) and the clock form factors are O(q^2),
the local-source forward contribution is finite as Lambda goes to zero:

    -(C/pi) Delta^2 Lambda+O(Lambda^3),
    (C Lambda^3/(3pi))[sum_removed ||D_b^perp v_b||^2
                          -sum_inserted ||D_b^perp v_b||^2]+O(Lambda^5).

Thus the actual patch response has no acoustic infrared obstruction
in this matched low-energy description. This does not bound its
ultraviolet contribution or prove the hydrodynamic remainder.

## 6. A weaker projected target than the unweighted supremum

The full supremum in the preceding contribution's (22) is sufficient,
but stronger than necessary. Let

    T_N=(1/N)sum_q gamma_N(q)=G_00,
    M_N=(1/N)sum_q 4sin^2(q/2)gamma_N(q),
    P_N=sum_(r=0)^9 |p_hat_r|gamma_N(2pi r/10).

Here p_hat_r=(1/10)sum_(j=0)^9 p_j exp(-2pi i r j/10), and N is
a multiple of ten so that the exterior target fits on the ring.
For each cell let Y_k be the partial sums of its centered source
vectors xi_s, and B_b=sum_k ||Y_k||. Summation by parts gives
F_b(q)<=4sin^2(q/2)B_b^2. Positivity and translation symmetry then
prove the EXACT sufficient second-order estimate

    |D_2|<=2||d||_1 P_N+||d||_1^2 T_N
                         +M_N sum_changed B_b^2.        (15)

Here changed includes removed and inserted clocks. No unchanged clock
is charged to a defect. Also M_N<=4T_N. Consequently uniform local
susceptibility T_N and the ten background-harmonic responses suffice,
even if an unweighted supremum were to fail in a narrow momentum region.
Equation (7) is an explicit exact determinant sum for these remaining
quantities. The q=0 term in P_N is now uniformly controlled by §3.
The nonzero harmonic values and the required uniform bound on T_N are
NOT certified by the bounded numerical table.

The actual original-clock sources give the following numerical
coefficients in |D_2|<=a_patch P_N+b_patch T_N after M_N<=4T_N:

| Replacement | ||d||_1 | a_patch | b_patch |
| --- | ---: | ---: | ---: |
| (10,10) -> (9,11) | 0.1494668170 | 0.2989336339 | 9.5740430650 |
| 9 length-10 -> 10 length-9 | 2.3120643489 | 4.6241286978 | 48.9101320840 |
| 11 length-10 -> 10 length-11 | 2.4157865728 | 4.8315731456 | 57.7885615909 |

These displayed coefficients are rounded, not outward-rounded interval
certificates. Their exact definitions are a_patch=2||d||_1 and
b_patch=||d||_1^2+4sum_changed B_b^2. The two exterior contacts are kept
distinct in computing d. The script records all ten |p_hat_r| and the
partial-sum norms; it never wraps a patch into an isolated periodic ring.

The original task of establishing/refuting the full bound (22) remains
open. Its exact outstanding mathematical content can now be stated as
a bound on (7), or the weaker moment conditions in (15), rather than
an unjustified free-fermion density bubble.

## 7. A genuine nonzero-source bound in a scalar reference

For the homogeneous scalar reference H_ch+uW, the same thermodynamic
XXZ formula gives an explicit finite-u bound:

    0<=-e_gap''(u)/2<=4/t,    |u|<=t.                   (16)

To see this put gamma in [pi/3,2pi/3]. Using coth(pi x)<=1+1/(pi x),

    I<=3log2/pi,
    |I'|<=(9log2+3)/pi^2,
    I''<=9/(2pi)+18log2/pi^3.

In e_X''(Delta)=e_(gamma gamma)/sin^2 gamma
                    -e_gamma cos gamma/sin^3 gamma,
these give |e_X''|<8; the exact upper expression is checked in the
script. Concavity supplies the sign. This is a thermodynamic
homogeneous-source result; it does not assert a finite-volume many-clock
curvature bound. Setting u=J e_10 treats the target's scalar mean
reference at nonzero J without any gap denominator. The periodic source
and operator-valued clock terms cannot simply be added to its energy,
and (16) does not make that substitution valid.
