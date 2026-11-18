# Strict negativity of every endpoint-moment coefficient

2026-11-18. **[Proved below; no independent human review.]**
This resolves the rational nonvanishing question in
[the eleven-state derivation, Section 6](../eleven_state_loop_splitting_20261118/DERIVATION.md)
for every finite n. The theorem is about those exact moments; its physical
application to the nine-state carrier is established separately in
[DERIVATION.md](DERIVATION.md).

## 1. Statement

For 0 <= i,j < n, define

    C_ij = (2i+2j+2)! / [(2i+1)!(2j+1)!],
    Q_ij = C_ij (i-j) / [2(i+j+1)],
    J_(n-2,n-1)=1, J_(n-1,n-2)=-1,

with all other entries of J zero. Let v_j be the unique vector supported
on j,...,n-1 with v_(j,j)=1 and C-orthogonal to every later coordinate.
Set d_j=v_j^T C v_j and

    L_j = sum_(l>j) (v_j^T Q v_l)(v_j^T J v_l)/d_l.

**Theorem.** For every n >= 2 and 0 <= j < l < n,

    sign(v_j^T Q v_l) = (-1)^(j+l),
    sign(v_j^T J v_l) = (-1)^(j+l+1).              (1)

Every summand in L_j is therefore strictly negative, and

    L_j < 0,  j=0,...,n-2.                         (2)

In particular there can be no cancellation in this sufficient
all-n nonvanishing test. As a normalization check,
L_(n-2)=-1/(16n-12), agreeing with the previously known last coefficient.

## 2. Reverse the polynomial variable, not merely the indices

Put b_i(t)=t^(2i+1)/(2i+1)! and
P_j(t)=sum_i v_(j,i)b_i(t). Then

    C_ij = integral_0^infinity b_i(t)b_j(t)e^(-t) dt,
    Q_ij = integral_0^infinity [b_i'(t)-b_i(t)/2]
                                      b_j(t)e^(-t) dt.

Consequently the P_j are mutually orthogonal. Define m=n-1-j and

    R_m(x)=t^(-(2n-1)) P_j(t),  x=t^(-2).

R_m is a degree-m polynomial with positive leading coefficient
1/(2j+1)!, and all the R_m are ordinary orthogonal polynomials for

    dmu_n(x)=(1/2) x^(-2n-1/2) exp(-1/sqrt(x)) dx,
    x>0.                                           (3)

Only this finite system through degree n-1 is needed. Moments through
degree 2n-1 exist; an infinite sequence with this fixed measure is
neither needed nor asserted. The usual elementary proofs of positivity
of the zeros and adjacent interlacing use only these finite moments.
For reference, the standard general results appear in
[NIST DLMF §18.2(vi)](https://dlmf.nist.gov/18.2#vi).
One can also obtain them from the positive three-term recurrence:
orthogonalize 1,x,...,x^(n-1); multiplication by x has positive adjacent
coefficients and positive quadratic form. Its successive Jacobi principal
matrices have simple positive eigenvalues with strict interlacing.

## 3. A signed resolvent identity

For any two different polynomials R_m,R_k in this finite system and s>0,

    I_mk(s)=integral R_m(x)R_k(x)/(x+s) dmu_n(x)
    has sign (-1)^(m+k).                            (4)

Here is a proof that does not assume a total-positivity claim. If m>k,
polynomial division and orthogonality give

    I_mk(s)=R_k(-s) integral R_m(x)/(x+s) dmu_n(x).

Dividing R_m(x)-R_m(-s) by x+s gives, again by orthogonality,

    integral R_m(x)^2/(x+s) dmu_n(x)
       =R_m(-s) integral R_m(x)/(x+s) dmu_n(x).

The integral on the left is positive. All zeros are positive and the
leading coefficients are positive, so R_r(-s) has sign (-1)^r. This
proves (4), including strictness, and symmetry handles m<k.

The elementary Stieltjes integrals

    x^(-1/2)=(1/pi) integral_0^infinity s^(-1/2)/(x+s) ds,
    x^(1/2) =(1/pi) integral_0^infinity s^(-1/2)x/(x+s) ds

now imply, for m != k,

    sign integral R_m R_k x^(-1/2) dmu_n =(-1)^(m+k),
    sign integral R_m R_k x^(1/2) dmu_n  =(-1)^(m+k+1).   (5)

For the second line, use x/(x+s)=1-s/(x+s) and orthogonality.
Interchanging the integrals is legitimate: both weighted absolute
polynomial products are integrable, as is seen directly from (3).
Equivalently, first truncate the s integral and then use the integrable
absolute product as a dominating function. The nonzero integrands in
(4) all have the same sign after the displayed factors are removed.

## 4. The derivative form has a strict checkerboard sign

Let V have rows v_j, let d=diag(d_j), and define

    A=V diag(2i+1) V^(-1),
    N_jl=integral P_j(t)P_l(t)e^(-t)/t dt.

The matrix A is upper triangular, A_jj=2j+1, because

    t P_j'(t)=sum_l A_jl P_l(t).

Integration by parts, whose endpoint terms vanish, gives

    A_jl d_l + A_lj d_j
      =integral (t-1)P_j(t)P_l(t)e^(-t) dt.

For j<l the second term is zero. The first line of (5), with t=x^(-1/2),
therefore proves

    sign A_jl=(-1)^(j+l), j<l,
    sign N_jl=(-1)^(j+l+1), j != l.                 (6)

The latter assertion uses the second line of (5), since 1/t=sqrt(x).
Both signs are strict. Integrating (P_j P_l)'e^(-t) also gives the exact
matrix identities

    AN+NA^T=d,
    V Q V^T=(AN-NA^T)/2.

For j<l this yields the particularly useful no-cancellation formula

    (V Q V^T)_jl = -sum_(i>=l) N_ji A_li.           (7)

Every term inside the sum has sign (-1)^(j+l+1): for i=l use
A_ll>0, and for i>l use both signs in (6). All are nonzero.
Equation (7) proves the first assertion of (1).

## 5. The endpoint skew form has the opposite sign

Let r=n-2, s=n-1, m=n-1-j and k=n-1-l. In the transformed polynomial,

    R_m(0)=v_(j,s)/(2n-1)!,
    R_m'(0)=v_(j,r)/(2n-3)!.

Thus, up to the positive factor (2n-3)!(2n-1)!,

    v_j^T J v_l
      =R_m'(0)R_k(0)-R_m(0)R_k'(0)
      =R_m(0)R_k(0) [R_m'(0)/R_m(0)-R_k'(0)/R_k(0)].

For m>k, strict interlacing of positive zeros implies

    sum_(zeros of R_m) 1/x > sum_(zeros of R_k) 1/x.

For adjacent degrees, each of the smaller m-1 zeros of R_m lies to the
left of its counterpart in R_(m-1), and R_m has one additional positive
zero; then iterate. The logarithmic-derivative bracket is therefore
strictly negative. The product R_m(0)R_k(0) has sign (-1)^(m+k).
This proves the second assertion of (1), and hence (2).

## 6. Verification and scope

`moments.py` performs direct exact reverse Gram-Schmidt, checks both
integration-by-parts identities and all the strict signs in (6)-(7), and
reproduces the old fractions for n=2,3,4. The bounded checks also use
n=5,6. These are independent arithmetic checks of the algebra, not the
reason for the all-n conclusion. No subagent or independent reviewer was
used. AI assistance: OpenAI Codex, requested GPT-6 Astra / Max.

The proof also resolves the same rational obstruction for the previous
eleven-state construction, if its separate physical reduction is retained.
It does not establish a many-body conformal limit for either model.
