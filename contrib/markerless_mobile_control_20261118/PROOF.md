# A finite mobile selector certificate and the remaining bulk comparison

2026-11-18. OpenAI Codex, GPT-6 Astra / Max assignment. Proofs and integer
certificates below have not had an independent human review. Floating-point
response coefficients are explicitly distinguished from the certificates.

## 1. The literal sector and ensemble

Use the corrected nearest-neighbour, Hermitian, one-site translation-invariant
ten-state bond in `../markerless_d9_d10/bond10.py`. Its material alphabet is
`L0,L1,R0,R1,W,A,B,XA,XB`; the tenth physical state is `V`. All material
energies, including both halves of the head chemical potential, are contact
weighted. The modified charge field is specified separately in FIELD.md.
No blocked or encoded dimension is substituted for physical d=10.

Fix the occupied number N and vacancy number V on a periodic physical ring
of length L=N+V. The occupied chemical potential is the common scalar -mu N
and is removed. A legal squeezed word has cells of physical material lengths
`w=(b_1,...,b_m)`, with sum b_c=N and one head per cell. All boundaries are
in ordered material coordinates. Fix the vacuum flag W on every cell and
initially set lambda=f=0. The data bits then give identical copies of the
same operator, so use all zeros. The extension to charged heads is in §6.

Let g_j be the number of vacancies after material particle j, let
A_j=1_(g_j=0), and let S_j raise g_j. The exact gap Hamiltonian is

    H_ch = -t sum_j (S_(j-1) S_j^* + h.c.),   sum_j g_j=V.       (1)

A cell starting at material site a and having b head positions r=0,...,b-1
has the following contact sources, denoted h_(b,s), s=-1,...,b-1:

    h_(b,-1) = (nu/2)|0><0|,
    h_(b,r)  = (nu/2)(|r><r|+|r+1><r+1|)
                            - |r><r+1| - |r+1><r|, 0<=r<b-1,
    h_(b,b-1)= (nu/2)|b-1><b-1|.                            (2)

Thus, with different cells acting on different clock factors,

    H_w = H_ch + J sum_c sum_(s=-1)^(b_c-1) A_(a_c+s) h_(b_c,s). (3)

Adjacent cells contribute separate halves to their common boundary contact.
In particular, the head energy is NOT nu times the head number in (3).
The sources in (2) agree entry by entry with the literal 81-by-81 bond on
words (2,3), (3,2,2), (10,10), and (9,11), with zero leakage. This check
keeps the endpoint contacts that a marker-based reduction could miss.

All off-diagonal entries of (3) are nonpositive. A rooted positive ground
vector can be summed over its material-root orbit to satisfy physical
cyclic gluing at momentum zero. Conversely, taking absolute values bounds
every twisted/glued sector below by the untwisted rooted minimum. Hence
the minimum over physical translations and momenta of this vacuum word is
exactly the rooted ground energy used here. No physical hopping link is
cut, and there is no O(t) boundary penalty. This argument is for the
stoquastic vacuum operator, not the complex finite-port Hamiltonian.

## 2. The frozen margin that must actually be compared

Put theta_b=pi/(b+1). The frozen cell ground vector and energy are

    v_(b,r) = sqrt(2/(b+1)) sin((r+1)theta_b),
    epsilon_b = nu - 2 cos(theta_b),
    e_b = epsilon_b/b.                                     (4)

For a target b, define a_b=2[1-cos(theta_b)] and
d_b=(b+1)a_b-b a_(b+1). The midpoint is

    nu_b = 2-(d_b+d_(b-1))/2
         = 2 cos(theta_b)+b[cos(theta_(b-1))-cos(theta_(b+1))]. (5)

Both adjacent wrong lengths have the SAME excess per cell:

    epsilon_(b-1)-(b-1)e_b
      = epsilon_(b+1)-(b+1)e_b
      = m_b := (d_(b-1)-d_b)/(2b).                         (6)

For n=9, b=10:

    nu_10 = 1.7702928472898476...,
    m_10  = 0.0020036046447729...,
    min(e_9-e_10,e_11-e_10) = m_10/11
                             = 0.0001821458767975....     (7)

The balanced replacement (10,10)->(9,11) costs 2m_10, not the energy
of deleting a head. The weak head-number replacements are

    (10 repeated 9 times)  -> (9 repeated 10 times), N=90,
    (10 repeated 11 times) -> (11 repeated 10 times), N=110. (8)

Each costs 10m_10. These replace the material inside a finite patch of
equal length and can be inserted into a larger all-10 word without moving
its exterior boundaries. At large target b, m_b=3 pi^2 b^(-4)+O(b^-5);
the nearest density gap is of order b^-5. Neither is the logical band's
normalized 2^-35 tolerance.

## 3. Exact finite-volume response with the vacancy fluid retained

The free gap ground state Omega is unique. Its contact expectation is

    alpha_(N,V) = N/(N+V)
      - sin^2(pi N/(N+V)) / [N(N+V) sin^2(pi/(N+V))].      (9)

Let Q_ch remove Omega and define the actual finite-volume susceptibility

    x_j=(A_j-alpha)Omega,
    G_ij=<x_i, Q_ch(H_ch-E_ch)^(-1)Q_ch x_j>.              (10)

G is real positive semidefinite and circulant. Its Fourier eigenvalues
are gamma(q)>=0. There is NO claim that the bath has a positive gap
uniformly in volume. Formula (10) is a finite-volume definition, not a
bound by a fictitious uniform inverse gap.

For each cell set

    a_(b,s)=<v_b,h_(b,s)v_b>,
    xi_(b,s)=(h_(b,s)-a_(b,s))v_b,
    sigma_j(w)=sum_(c,s:a_c+s=j) a_(b_c,s),
    F_b(q)=||sum_(s=-1)^(b-1) exp(iqs) xi_(b,s)||^2,
    R_b=(1/N)sum_q gamma(q) F_b(q).                      (11)

The exact cancellation sum_s xi_(b,s)=0 follows from (4). Exciting
different cells gives orthogonal clock states. Degenerate perturbation
theory first compresses onto Omega; its first-order ground vector is the
product of the v_b. It therefore gives, at FIXED finite N,V,w,

    E_w(J)=E_ch+J alpha sum_c epsilon_(b_c)
        -J^2 [sigma(w)^T G sigma(w)+sum_c R_(b_c)]
        +O_(N,V,w)(J^3).                               (12)

This is a scalar gap response plus explicit single-clock form factors.
It does not require diagonalizing a product of ten clocks to calculate
the coefficient for either word in (8). Its remainder is not proved
uniform in volume.

There is a particularly simple new closed form for the markerless clock.
Writing a_b(q)=sum_s a_(b,s)exp(iqs), one has

    F_b(q) = epsilon_b^2 cos^2(q/2)
        + 4 sin^2(theta_b)(b-1)/(b+1) sin^2(q/2)
        - |a_b(q)|^2.                                  (13)

To prove it, set v_(-1)=v_b=0. Component r of sum_s exp(iqs)h_s v is

    exp(iq(r-1/2)) [epsilon_b cos(q/2)v_r
                         - i sin(q/2)(v_(r+1)-v_(r-1))].

The real and imaginary terms have zero cross term in its norm. Also
sum_r(v_(r+1)-v_(r-1))^2=4 sin^2(theta_b)(b-1)/(b+1).
Orthogonal projection away from v subtracts |a_b(q)|^2, proving (13).
In particular F_b(0)=0 exactly. This is not the biased frustration-free
clock or exponentially small terminal tail used by the earlier d=22 work.

### Exact one-vacancy reduction

For V=1 the bath is the N-cycle adjacency operator. Direct diagonalization
by Fourier transform gives

    gamma(0)=0,
    gamma(2 pi l/N)=1/[4t N sin^2(pi l/N)], l!=0.       (14)

For an isolated cell let Y_k=sum_(j=0)^k xi_(b,j-1), k=0,...,b-1,
the first b partial sums of the ordered list of b+1 sources. Define
U_b=sum_k||Y_k||^2 and W_b=||sum_k Y_k||^2.
The inverse cycle Laplacian then yields the exact identity

    R_b = U_b/(tN)-W_b/(tN^2),                          (15)

provided its b+1 contacts do not repeat the same material bond. Thus the
first nontrivial correction in (8) is a finite trigonometric calculation,
not a many-clock approximation. V=1 while N grows approaches full
occupancy; these numbers are NOT a finite-density thermodynamic test.

## 4. The first nontrivial comparisons

Here nu is the exact midpoint (5), and t=1. For E_B-E_A write
J D_1+J^2(D_2,scalar+D_2,clock)+O(J^3). The following evaluations are
[numerical]; the formulas defining them are exact.

| Target -> competitor | N,V | D_1 | D_2,scalar | D_2,clock |
| --- | ---: | ---: | ---: | ---: |
| (10,10) -> (9,11) | 20,2 | 0.00360722618 | -0.00000540520 | -0.00016772686 |
| (10,10) -> (20) | 20,2 | 0.08103222951 | +0.00206290471 | +0.00757527689 |
| (10,10) -> (6,7,7) | 20,2 | 0.09974839075 | -0.00178435480 | -0.01580460037 |
| 9 length-10 cells -> 10 length-9 cells | 90,1 | 0.01981342371 | -0.00014401670 | -0.00123127918 |
| 11 length-10 cells -> 10 length-11 cells | 110,1 | 0.01985390057 | +0.00009807454 | +0.00076196778 |

The balanced defect and shorter-cell head change are helped by the
vacancy correction. The longer-cell head change is hurt by it. Thus the
favorable all-length clock-response sign of the earlier biased d=22
clock does not transfer to this model. None of the signs in this table
alone gives a finite-J or thermodynamic inequality.

## 5. An integer certificate at n=9, with a finite J/t window

The primary calculation is canonical: N=20, V=2, physical L=22,
occupied density rho=10/11. Choose t=1, J=1/25, nu=177/100. The matrices
5000 H_w have integer entries: vacancy hops -5000, head hops -200,
and +177 for each occupied contact half incident on a head. Their sizes
are 21,000 for (10,10), 20,790 for (9,11), 4,200 for (20), and 61,740
for (6,7,7). These are four bounded one-worker calculations, not a scan.

For any strictly positive integer vector u and real stoquastic H,

    min_i (Hu)_i/u_i <= E_0(H) <= max_i (Hu)_i/u_i.       (16)

For example, shift -H by a sufficiently large scalar and apply the
Perron comparison inequalities; irreducibility is unnecessary for the
two bounds. Floating eigensolvers supply only u. The certificate checks
the bounds using integer matrix-vector multiplication and exact integer
division, rounded outward to multiples of 10^-10. It checks positivity
and a row-sum overflow bound first. Thus eigenvector residuals and
floating-point eigenvalues are NOT the certification step.

The actual certified lower bounds on E_B-E_(10,10), in units t, are:

| Competitor | At nu=1.77 | Uniformly for 1.77<=nu<=1.7703 |
| --- | ---: | ---: |
| (9,11) | 0.0001440227 | 0.0001200227 |
| (20) | 0.0032671549 | 0.0032431549 |
| (6,7,7) | 0.0039518859 | 0.0039278859 |

The exact midpoint lies in this rational nu interval. `checks.py` proves
the enclosure using rational alternating arctangent bounds in Machin's
formula for pi and rational alternating cosine bounds. No rounded cosine
is used in the certificate. To extend from 1.77, increasing nu raises
every competitor's Hamiltonian, while the target energy increases by at
most 2J(1.7703-1.77)=0.000024. This proves the last column.

There is also a small certified open parameter neighborhood. Since
||sum_j A_j h_j|| <= m(nu+2) for a vacuum word of m cells, the last column
loses at most (m_B+2)3.7703 |J-1/25| when J changes. Therefore, throughout

    |J/t-1/25| <= 1/200000,   1.77<=nu<=1.7703,         (17)

the balanced gap remains above 0.0000446167 t, and the other two gaps
also remain positive. This narrow rigorously enclosed neighborhood is a
conservative verification allowance, not evidence of physical fine
tuning. We do NOT claim that every 0<J/t<=1/25 is certified.

The witnesses are saved in `witnesses.npz` and can be checked without
running an eigensolver. This is a finite comparison of three specified
competitors. It neither enumerates all partitions of 20 nor covers a
defect in an arbitrarily long bath, nor establishes half-occupancy order.

## 6. Charged-head domination and a finite-port extension

There is an additional exact zero-port simplification specific to this
alphabet. At f=lambda=0, assume

    D>=0,   Delta/2 >= max(1, |g|, |t_d|).              (18)

For EVERY legal word, any combination of charged/vacuum cell flags has
ground energy at least that word's all-vacuum ground energy, at every
J,t>0 and every N,V. This is a domination statement, not a quantitative
charged gap.

Proof: over each gap configuration and set of head positions, take the
norm of the wavefunction over rail, species, and payload labels. Each
head hop and each vacancy hop induces a unitary matching or contraction
between these fibers, with modulus at most the corresponding vacuum
amplitude. Cauchy--Schwarz bounds its energy below by the vacuum hopping
of the fiber norms. Complex express/parity phases can only weaken that
lowering. The nu contact term is identical in the two problems.

The remaining within-fiber operators are the left-end rail switch (norm
1 times the left contact), the right-end closing switch (norm |t_d|
times the right contact), and the conversion (norm |g| times the right
contact). The latter two cannot occur at the same head position: one is
at the far right, the other requires a data letter to its right. The
extra charged mass is (Delta/2)(A_left+A_right); D is nonnegative.
Condition (18) pays each within-fiber norm separately. Summing the
result proves the quadratic-form lower bound after taking fiber norms.
Boundary switches are identity on the neighboring cell, as required by
the exact legal grammar. In a physical cyclic representation the norm
map discards all phases and lands in the admissible vacuum zero-momentum
comparison. Conversely, the vacuum ground state is an allowed trial.

The certified markerless schedule has Delta=4-2cos(pi/11)>2, D=1000,
|g|=2^-163, and tiny |t_d|, so it satisfies (18). A uniform strictly
positive charged excitation cost does NOT follow: a head can have both
contacts absent in individual vacancy configurations.

As an additional literal check, the 48-state charged cell at b=3 is
compared against its vacuum clock for all eight binary contact masks,
with Delta=4-2cos(pi/11), D=0, g=1/2, t_d=2/5 and nontrivial phases.
Every difference is nonnegative. This bounded check supports the
fiber identification; it is not the proof of the all-word statement.

Restore a port and the external field J f_ext Q/2, with f_ext>=0. There
are at most m_B incident cell boundaries in a whole word B. Each port
has norm sqrt(2)|lambda|. Also Q>=-m_B, since only charged heads have
negative charge. The target vacuum all-zero-bit trial has zero field
and zero port expectation. Consequently the finite certificate extends
to the complete charged/bit spaces of these particular words as

    E_B(full)-E_A(full) >= gap_B(vacuum)
                        -J m_B[sqrt(2)|lambda|+f_ext/2]. (19)

At J/t=1/25 and 1.77<=nu<=1.7703, the simple allowance

    sqrt(2)|lambda|+f_ext/2 <= 1/1000                  (20)

leaves every comparison strictly positive; the balanced lower bound is
0.0000400227 t. The physical band-schedule lambda and field are vastly
smaller. A nonnegative contact field instead obeys the same norm estimate
with f in place of f_ext, although it retains its position-dependent
logical error. These are finite whole-word bounds. Applying the m_B
cost in (19) to a defect with infinitely many unchanged spectator cells
would charge an extensive energy to that defect and is NOT justified.

Illegal squeezed words are not compared, nor are every legal word and
head count. Hence even this finite-port extension is not a theorem that
the unrestricted microscopic ground state has the intended pattern.

## 7. The precise relative-response estimate still needed

Fix a density, for example rho=1/2 with V=N, and target b. Let A be an
all-b word on a compatible ring and B a same-length local replacement,
such as (b,b)->(b-1,b+1), or either head-number patch in (8). Outside the
patch the material sources are unchanged. Put d=sigma_B-sigma_A. Then
the EXACT second-order excess coefficient is

    D_2 = -2 d^T G sigma_A - d^T G d
          +sum_(removed cells) R_b -sum_(inserted cells) R_k. (21)

There is no extensive remainder from unchanged cells in (21). In
particular the scalar response does not disappear merely because the
material and vacancy numbers are equal.

For a useful explicit sufficient coefficient bound, suppose

    sup_(N,q) t gamma_(N,V)(q) <= Gamma_rho < infinity.   (22)

This is an UNPROVED susceptibility hypothesis here. Let p be the b-site
period of sigma_A and let p_hat_r=(1/b)sum_j p_j exp(-2pi i rj/b).
Circulant symmetry gives

    ||G sigma_A||_infinity <= (Gamma_rho/t) sum_r |p_hat_r|.

Using Parseval on (11), define

    C_patch = 2||d||_1 sum_r |p_hat_r| + ||d||_2^2
                   + sum_(changed cells, both words) sum_s ||xi_(k,s)||^2.

Then |D_2|<=Gamma_rho C_patch/t, uniformly in the number of unchanged
cells, if (22) holds. The last sum has the explicit closed form

    sum_s ||xi_(k,s)||^2
      = [epsilon_k^2+4 sin^2(theta_k)(k-1)/(k+1)]/2
                                               -sum_s a_(k,s)^2. (23)

A sufficient ALL-ORDER replacement for this coefficient statement is
the following relative curvature inequality, in the zero-port problem:

    |d^2/dJ^2 [E_B(J)-E_A(J)]|
                          <= 2 Gamma_rho C_patch/t     (24)

for 0<J<=J_0, uniformly in compatible canonical volumes. If there are
level crossings, the corresponding integrated second-difference bound
must be used instead. Integrating (24), with the right derivative fixed
by (12), gives

    E_B(J)-E_A(J) >= J alpha Delta_patch
                             - Gamma_rho C_patch J^2/t. (25)

Thus J/t < alpha Delta_patch/(Gamma_rho C_patch), within the interval
of (24), would be a genuine sufficient allowance for this patch. The
finite coefficient calculation checks only the curvature at J=0; it
does NOT prove (22) or (24). This formulation keeps the gapless fluid,
subtracts the same interacting bulk, and states the missing estimate
without an artificial vacancy gap.

The n=9 constants can already be evaluated without any vacancy-state
diagonalization. The following [numerical] values use unwrapped patches,
so their two exterior contacts remain distinct:

| Replacement | C_patch | 0.29 Delta_patch/C_patch |
| --- | ---: | ---: |
| (10,10) -> (9,11) | 0.561732997 | 0.002068760 |
| 9 length-10 -> 10 length-9 | 3.108260770 | 0.001869358 |
| 11 length-10 -> 10 length-11 | 3.021884595 | 0.001922791 |

At V=N and N>=20, alpha_(N,N)>0.29. For example this follows from
sin x>=x-x^3/6 in (9), with x=pi/(2N)<=pi/40. Thus the last column
would be an illustrative allowance on Gamma_(1/2) J/t if (24) were
proved. It is not an unconditional numerical J/t budget. Keeping the
distinct exterior contacts is essential when using this table for a
patch inside a larger ring rather than an isolated cyclic word.

The dependence on n is at least identifiable. At the selector midpoint
and k=b+O(1), a_(k,s)=O(b^-3), epsilon_k=O(b^-2), and (23)=O(b^-2).
The balanced patch has C_patch=O(b^-2) by this conservative bound and
Delta_patch=2m_b=Theta(b^-4). A head-number patch has O(b) changed
cells, C_patch=O(b^-1), and Delta_patch=b m_b=Theta(b^-3). Both would
therefore admit a sufficient scale J/t=O(b^-2/Gamma_rho), IF a
uniform all-order estimate of the form (24) held with comparable
constants. This is a conditional power count, not an established
mobile schedule. Actual cancellations in (21) can improve it.

The decisive next calculation is to establish or refute (24) for the
balanced and shorter-cell head-number patches at V=N, with fixed positive
J/t. A scalar/contact response estimate alone at J=0 is insufficient.
Even that result would leave other wrong lengths, illegal sectors,
finite-port relative errors, dressed logical bands, and the unitary
compact Virasoro CFT identification as separate tasks.
