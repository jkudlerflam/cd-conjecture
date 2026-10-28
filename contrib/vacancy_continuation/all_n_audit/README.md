# Independent audit of C39's all-n finite input

**[checked analytical argument; exact finite checks]** I found no error in
C39's generalized-eigenvalue, slope-separation, or pair-sum argument.
The order of its choices is sufficient for each fixed positive integer
n. The result does not require a lower bound uniform in n. An independent
exact verifier checks the determinant, simple roots and all unordered
pair comparisons through n=4. The n=4 reference certificate is new here.

This supports the finite band input for the separate mobile-completion
effort. It does not independently verify the inherited full selector,
prove the infrared transfer, or identify 46 with the alphabet of a
complete spontaneous-cell chain. C39 explicitly counts 46 initialized
core letters and 4,464 letters after its inherited wrapper.

## Claim and hypotheses

**[proved analytical lemma, checked application]** Fix finite n, a positive
definite real matrix A, and a real analytic symmetric matrix C(v) near
v=0 with C(0)>0. Let

    D(v)=diag(1,v^2,...,v^(2n-2)),
    B(v)=D(v) C(v) D(v),
    h(k)=-2 A cos(k)+2 eta B(v) sin(k)+mu diag(u^(4j)).

For sufficiently small positive eta, then u, then v, then mu, h has
exactly 2n simple Fermi roots, positive minimum speed, distinct sea
centers, and distinct unordered pair sums modulo 2pi, including pairs
with repeated labels. C39 uses A=C(w), fixed w=1/4, and C(v)=C(wv),
whose positive Gram limits satisfy the hypotheses. The mathematical
statement is about a finite quadratic reference, not its central charge
after microscopic interactions.

## Why the quantifiers close

**[proved; independent derivation of the core steps]** In the coefficient
of t^(n-k) in det(B(v)-t A), a B minor on row set I and column set J has
explicit power v^(2 sum(I)+2 sum(J)). The unique minimum is
I=J={0,...,k-1}. Its coefficient is the product of the leading principal
C(0) minor and complementary trailing A minor, so it is strictly
positive. Dividing by det A gives

    e_k(r)=c_k v^(2k(k-1))(1+o(1)),
    c_k=det C(0)[0:k] det A[k:n] / det A.

At scale t=v^(4j)y, only terms k=j and k=j+1 have minimum power.
Their nonzero root y=c_(j+1)/c_j is simple. This supplies all n positive
generalized eigenvalues

    r_j=kappa_j v^(4j)(1+o(1)),  kappa_j=c_(j+1)/c_j>0.

No cancellation, total positivity, or unproved genericity assumption is
needed. Ordinary positive definiteness supplies the required signs.

**[proved]** Eliminate the first j coordinates of B-r_j A. Their inverse
has norm O(v^(-4j+4)), while B_II^(-1)B_IT=O(v^2) and
r_j B_II^(-1)=O(v^4). After division by v^(4j), the tail equation is a
positive rank-one matrix against kappa_j A_TT. Its unique normalized
limit is A_TT^(-1)e_0 divided by its first coordinate. Thus the limiting
x_j has zero entries below index j and x_(j,j)=1.

The field derivative at the two roots theta_j plus/minus pi/2 is
minus/plus

    s_j=(x_j^T diag(u^(4i)) x_j)
        / [2 sqrt(1+eta^2 r_j^2) x_j^T A x_j].

At v=0 this is a positive constant times u^(4j)(1+O(u^4)). Choosing
u first makes every limiting slope ratio less than 1/64. With that u
fixed, choosing v sufficiently small makes every actual slope ratio
less than 1/32, while simultaneously separating all center ratios.
One must not set v=0 in the original singular eigenproblem or require
the slope limit to hold uniformly as u changes. C39 does neither.

**[proved]** For any two distinct unordered pairs, let h_(j,+/-) be their
signed label multiplicity difference and m_j=h_(j,+)+h_(j,-). Then
sum|h|<=4. Their zero-field pair difference is a lacunary center
combination plus an integer multiple of pi. If that integer is odd,
the center term, of magnitude below one, cannot cancel it modulo 2pi.
If it is even and m is nonzero, the first nonzero lacunary center
coefficient dominates the tail. If m=0, distinct label multisets imply
a nonzero integer slope combination with total absolute coefficient
at most four. That splits linearly in mu. This covers repeated-point
and same-flavor comparisons as well as four-distinct-point processes.

There are only finitely many comparisons at fixed n. Analytic root
branches have finite first and second derivative bounds on a common
neighborhood. One sufficiently small positive mu therefore retains all
nonzero constant gaps and splits every zero constant gap. The physical
matrix is D_w h D_w; its simple-root velocity is positive because D_w
is invertible. A congruence is used for root locations, not velocities.

**[checked parameter dependency]** The complete order is

    n,w -> eta -> u -> v -> mu -> D -> g -> lambda.

The first five choices fix a strict normalized reference neighborhood.
The large-D rational coefficient tends to that reference. After D is
finite, the normalized g^4 remainder tends to zero as g^2, though its
constant can grow severely with D and n. The code gap and readout norm
are finite when lambda is finally chosen. Shrinking lambda alone would
not reduce the g^4/g^2 ratio. C39 preserves the needed order and allows
all constants to depend on n. Its specific small-n parameter table is
not asserted to work at all n.

## Exact certificate checks and extension

**[exact finite checks]** `regenerate.py` reuses C39's rational Gram and
Sturm certificate generator, writing only in this folder. It reproduces
n=1,2,3 and extends its displayed parameter sequence to n=4:

    w=1/4, v=1/1024, u=1/64, eta=1/100, mu=10^(-(8n+10)).

`verify.py` then checks each saved certificate with different algorithms:

1. Expand the polynomial determinant directly as a permutation sum of
   products of quadratic polynomials, independently of C39's interpolation.
2. Check opposite exact signs at both endpoints of every saved interval.
   There are 2n disjoint such intervals for degree 2n, so all roots are
   simple and accounted for. The positive leading coefficient excludes
   the missing Cayley point k=pi.
3. Check all unordered pair comparisons using rational interval arithmetic
   on Cayley variables, independently of C39's arctangent implementation.

For x=tan(p/2), y=tan(q/2), z=tan(r/2), w=tan(s/2), let

    F=(x+y)(1-zw)-(z+w)(1-xy).

Then the sine of half the pair-sum difference has magnitude

    |sin((p+q-r-s)/2)|
      = |F| / sqrt((1+x^2)(1+y^2)(1+z^2)(1+w^2)).

Every interval for F excludes zero. A deliberately loose rational
denominator bound yields the following independent positive lower bounds
on angular separation. Displayed numbers are rounded down; JSON retains
the exact rational bounds.

| n | Simple roots | Pair comparisons | Independent delta4 lower bound |
|---|---:|---:|---:|
| 1 | 2 | 3 | 1.2798e-16 |
| 2 | 4 | 45 | 3.0909e-29 |
| 3 | 6 | 210 | 5.9076e-42 |
| 4 | 8 | 630 | 1.0437e-54 |

The n=4 generator's sharper angular bound is 4.1748e-54, and its
root and center separations both exceed 8.4678e-27. The independent
verifier intentionally sacrifices a constant factor to avoid square
roots and all trigonometric series.

**[scope]** The verifier inherits the rational Gram matrices from C39's
evaluator. It independently checks the determinant and the root-to-pair
certificate, not the microscopic derivation of that Gram formula. The
all-n result follows from the analytical argument above; n=4 does not
serve as evidence by extrapolation. None of these checks measure a
central charge or certify the finite-D diagnostic parameters used in
C39's many-body ED.

## Reproducing and provenance

From the repository root, using Python's standard library only:

    python3 contrib/vacancy_continuation/all_n_audit/regenerate.py
    python3 contrib/vacancy_continuation/all_n_audit/verify.py

The replay through n=4 took about 10 seconds here; independent verification
took about 0.12 seconds. Machine speed affects these timings. Saved data
are `certificate_n1.json` through `certificate_n4.json` and
`verification.json`. No C39 source or original certificate was changed.

Author and checker: an independent GPT-6-family Codex research worker,
2026-10-28, reviewing the C39 collaborator's argument. The separate
verification code and the analytical review were written for this audit.
No human review is claimed. The conclusion strengthens the finite input
available to C43; it supplies neither C43's selector nor its infrared
conclusion.
