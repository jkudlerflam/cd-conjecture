# Connected coefficients, readout distortion and fermionic normalization

2026-11-18. **[Proved identities and bounds; numerical literal checks,
not independently reviewed.]** The physical model is exactly
[the modified nine-state bond](../markerless_all_n_20261118/model.py).
No physical term or state is added. All statements below hold at fixed
n and fixed couplings, independently of the number of selected cells.

## 1. Expansion convention

Use a selected frozen sector with n data letters and one head per cell.
Its full cell Hilbert space has dimension 5(n+1)2^n; its code has dimension
2^n. These are cell bookkeeping dimensions, not the physical dimension,
which remains d=9. Subtract the constant selected zero-port cell energy.

Take g into the single-cell Hamiltonian exactly and write

    H(lambda)=H_0(g)+lambda V+lambda^2 F,
    H_0(g)=sum_x h_x(g),  h_x P_x=0,  h_x>=delta Q_x,
    delta=1/2,
    V=sum_x V_(x,x+1),  ||V_(x,x+1)||=b=sqrt(2),
    F=chi sum_x N_(L1),x,
    chi=g^2 zeta D^(-2n-2) epsilon/q^2,
    zeta=|v_0 v_n|^2.                                (1)

Here v is the actual biased W-clock ground vector. All the earlier
band parameters, including the loop and epsilon, have already been fixed
in the proved order. Assume q<=1/4, w<=1/2, Delta=3, |g|,|t_d|<=1/4,
D>0. The charged pair produced by a port has gap at least G=3; at g=0
one may use G_0=4. Each cell conversion has norm at most one. Nu_head
is a constant on a selected fixed-head sector and has been subtracted,
not changed. The exact readout coefficient is h_L=lambda^2 chi.

On a finite connected cluster put P=product_x P_x and R=Q H_0^(-1) Q.
P contains every logical Fock state. A port source has norm

    ||V_e P||=beta=sqrt(2 zeta),                      (2)

because both head endpoints must be occupied by the vacuum clock.
The source norm is distinct from the full port norm b.

The following effective operators use the Hermitian direct-rotation
frame: the overlap between the exact low subspace and P is chosen
positive. Section 4 of UNIFORM.md relates it to a local normal form.
This frame fixes the otherwise ambiguous fourth-order code rotation.

## 2. The first conversion-generated four-fermion operator

At order lambda^2 and with the field contribution separated,

    K_2(g)=-P V R(g) V P.

This is a sum of independent nearest-cell bond terms: two different
ports cannot return a vacuum source to P at this order. With C the sum
of the two unit conversions and R_0 the charged-pair inverse at g=0,

    K_2(g)=K_20+g^2 K_22+g^4 K_24+...,
    K_24=-P V R_0 C R_0 C R_0 C R_0 C R_0 V P.       (3)

Odd coefficients vanish by charge conservation. K_20 is scalar, and
K_22 is the already proved quadratic hopping/on-cell operator. Four
conversion insertions give at most four logical fermion operators after
the established parity gauge, so K_24 is scalar plus quadratic plus
quartic. Its quartic part must not be replaced by its one-particle block.

For modes in physical cell order, our normal-order convention is

    H=E_vac+sum_ij h_ij c_i^dagger c_j
       +sum_(i<j,k<l) U_(ij,kl) c_i^dagger c_j^dagger c_l c_k+... .  (4)

In the two-particle basis |ij>=c_i^dagger c_j^dagger|0>, the matrix U
is exactly the two-particle block after subtracting the vacuum scalar
and the full Fock lift of h. There is no extra 1/4 in this restricted-pair
sum. In an unrestricted antisymmetric-index sum the coefficient is 1/4.
`fock_extract` reconstructs every Fock sector, including three and four
particles in the n=2 example.

At the declared diagnostic point q=.21, w=.4, Delta=3, D=2, t_d=.08,
with the matched phases and g treated as the expansion variable:

| n | Literal even-token states, two cells | Norm of quartic part of K_24 |
| ---: | ---: | ---: |
| 1 | 272 | 6.5151204922e-9 |
| 2 | 2448 | 6.8295373915e-11 |

For n=1 the coefficient of n_1 n_2 is positive, 6.5151204922e-9.
For n=2 the tensor has both signs and nonzero off-diagonal entries; its
entire 6-by-6 pair matrix is saved in checks.json. The largest negative
density entry is about -6.8295373915e-11. These numbers multiply the
**physical** factor lambda^2 g^4. The n=2 remainder after complete
scalar/quadratic/quartic Fock reconstruction has norm 3.1e-23.
These are numerical nonvanishing checks, not interval certificates or
claims about a phase at these diagnostic parameters.

There is a bound on the complete operator, in every logical particle
sector. With beta^2=2zeta, ||C||<=2 and ||R_0||<=1/4,

    ||K_24|| <= zeta/32,
    ||K_2(g)-K_20-g^2 K_22||
       <= zeta g^4/[32(1-g^2/4)].                  (5)

After retaining g^4, the corresponding g^6 tail is at most
zeta g^6/[128(1-g^2/4)]. The finite-g inverses at g=.16,.08 obey this
bound in both examples. No one-particle extraction factor is used in (5).

## 3. The readout creates an attractive interaction for every n>=2

Write D_i=diag(1_(j>i)) on clock positions j=0,...,n. The exact field on a
fixed bit word b is h_L sum_i b_i D_i. Let

    R_v=(H_v-e_v)^+,  gamma_i=<v|D_i|v>,
    C_ij=<v|D_i R_v D_j|v>.

Ordinary nondegenerate clock perturbation gives, as an operator on the
whole code,

    h_L sum_i gamma_i n_i
       -h_L^2 [sum_i C_ii n_i+2 sum_(i<j) C_ij n_i n_j]+O(h_L^3).  (6)

**[Proved]** Every C_ij is strictly positive at q>0. A useful exact
formula avoids a spectral sum. Put pi_j=v_j^2, P_j=sum_(s<=j) pi_s and
c_j=q v_j v_(j+1)>0. The ground-state transform of H_v-e_v has Dirichlet
form sum_j c_j(f_j-f_(j+1))^2. For the centered threshold D_i-gamma_i,
the cumulative source current is

    J_i(j)=sum_(s<=j) pi_s (1_(s>i)-gamma_i)
          =-P_min(i,j) [1-P_max(i,j)].

Solving the one-dimensional current equation gives

    C_ik=sum_(j=0)^(n-1) J_i(j)J_k(j)/c_j >0.       (7)

All factors inside the products are strictly negative. This proves
attractive density interactions in (6) for every distinct pair of flavours,
not just the small examples. It does not prove an instability: the
coefficient is h_L^2 and must be compared with the much larger first-order
field, kinetic term and other interactions. Positive h_L still leaves the
all-zero microscopic clock minimum unchanged, which is used in DEFECTS.md.

For n=2,q=.21, C_01=0.0031694422820111. The density-pair coefficient is
-0.0063388845640222 h_L^2. Direct literal-cell projection and (7) agree
within 5.2e-16; the n=3 check is within 2.1e-15.

The source visibility also changes, rather than remaining equal to zeta.
For a bit word b, the normalized vector correction is

    v(b,h_L)=v-h_L sum_i b_i R_v D_i v+O(h_L^2),
    d zeta(b,h_L)/d h_L at 0
       =-2 zeta sum_i b_i [(R_v D_i v)_0/v_0
                          +(R_v D_i v)_n/v_n].      (8)

This exact formula is checked by re-diagonalizing the perturbed clock.
A small absolute vector correction need not give a comparably small
relative source correction when v_n is tiny. The mixed terms below retain
this effect without dividing a norm bound by a tiny endpoint amplitude.

## 4. All fourth-order port and field terms in one formula

Define B=P V R V P, A=P V R^2 V P, and F_0=P F P. The Hermitian expansion
through fourth weighted order is

    H_eff=lambda^2 (F_0-B)+lambda^4 K_4+O(lambda^6),

    K_4^port=-P V R V R V R V P + (AB+BA)/2,

    K_4^mixed=P V R F R V P+P F R V R V P+P V R V R F P
                      -(A F_0+F_0 A)/2,

    K_4^field=-P F R F P,
    K_4=K_4^port+K_4^mixed+K_4^field.               (9)

R excludes only the full code. In particular it retains intermediate
states with **no charged head but an excited W clock**. Dropping those
states would remove the connected three-cell term and distort (8).
The anticommutators are the normalization terms; omitting them gives a
non-Hermitian Bloch coefficient and spurious disconnected interactions.
Formula (9) includes the denominator shift from the readout, its endpoint
source distortion, and its own clock-mediated interaction (6).

An independent check uses the intermediate-normalization wave operator
Omega=P+sum_(r>=1) lambda^r Omega_r, P Omega_r=0:

    E_r=P[V Omega_(r-1)+F Omega_(r-2)],
    Omega_r=-R[V Omega_(r-1)+F Omega_(r-2)
                         -sum_(j=1)^(r-1) Omega_j E_(r-j)].

Negative subscripts are zero. Its Hermitian fourth coefficient is
E_4+[Omega_1^dagger Omega_1,E_2]/2, reproducing (9).
On the literal three-cell example the two algorithms agree within 1.4e-17.
An independent direct eigensolver and positive-overlap rotation on two
cells give errors 8.31e-10 and 1.30e-11 at lambda=.12,.06, respectively:
the ratio is the expected sixth-order factor, with residual/lambda^6
approximately 2.79e-4. The diagnostic field scale chi=.2 in this check
deliberately makes the field terms visible; it is not the small band field.

## 5. The first connected three-cell port operator

On an open line, each cell changes between W and charged whenever an
incident port acts. A returning word has even incidence at every vertex,
so every used edge of the line occurs an even number of times. At port
order four the only connected possibilities are one edge four times or
two adjacent edges twice each. Disconnected two-edge contributions cancel
in the normalized effective operator by additivity. Thus there is no
bulk four- or five-cell port term at this order.

For adjacent edges a,b, put A_e=P V_e R^2 V_e P and B_e=P V_e R V_e P.
The connected operator, already including normalization, is

    C_abc=-P[V_a R V_a R V_b R V_b
             +V_b R V_b R V_a R V_a]P
             +{A_a,B_b}/2+{A_b,B_a}/2.              (10)

In each first term the middle inverse is a W-clock inverse. At g=.13
and the n=1 diagnostic parameters, its full norm is 1.1808413907e-4,
dominated by a scalar. It also has a nearest-cell quadratic entry about
1.6280710828e-11 and density-pair entries about -4.2580e-14. The full
8-by-8 matrix is saved, including the smaller higher-body residual; no
entry is discarded on the assertion that it is scalar or quadratic.
Subtracting the two embedded edge-only fourth coefficients agrees with
direct edge switching on the same literal cluster to 3.8e-18.

Small periodic rings can have winding words: the two-cell ring already
has the familiar extra order-lambda^2 winding, and a four-cell ring can
have order-lambda^4 winding. Equations (9)-(10) specify the bulk open-cluster
expansion, or rings with more than four cells at this order. The general
local-norm theorem includes finite-ring wrapping terms; it does not mistake
those small rings for the bulk coefficients.

Useful bounds retaining the endpoint suppression are

    B_2=4zeta/(G^2 delta)+4zeta^2/G^3,
    ||K_4^port on one edge||<=B_2,
    ||C_abc||<=B_3=2B_2.                            (11)

They follow from the charged outer inverses, the clock middle inverse,
||V_e||=sqrt(2), and ||V_e P||=sqrt(2zeta). For the two-cell mixed readout
and the one-cell field square one may use

    B_mix=8n |chi| zeta/G^2+8n |chi| sqrt(zeta)/(delta G),
    B_F=n^2 chi^2/delta.                             (12)

These are bounds on full operators and include arbitrary finite-g
nonquadratic terms. Their weighted local sum is used in UNIFORM.md.
