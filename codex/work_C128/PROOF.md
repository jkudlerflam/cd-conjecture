# C128 proof and limitations

2026-11-06. Natural logarithms, lattice spacing one, hbar one. Tags apply to the paragraphs and displays that follow them. The interacting crossed majorization is open. A failure of independent column clocks below is not a counterexample to the scalar target.

## 1. Histories, gates, and the two different notions of positivity

[Proved] Consider h=Sx tensor Sx+Sy tensor Sy+Delta Sz tensor Sz. Rotate alternate sites by pi about y. On even rings this gives the translation-invariant, stoquastic representative

    h_y=-Sx tensor Sx+Sy tensor Sy-Delta Sz tensor Sz
       =diag(-Delta/4,Delta/4,Delta/4,-Delta/4)
          -(1/2)(|00><11|+|11><00|).

Delta=0 is XX and Delta=1 is Heisenberg. Write t=beta/m and

    a=exp(t Delta/4) cosh(t/2),
    b=exp(-t Delta/4),
    c=exp(t Delta/4) sinh(t/2).

The gate g=exp(-t h_y), in the order 00,01,10,11, is

    [ a  0  0  c ]
    [ 0  b  0  0 ]
    [ 0  0  b  0 ]
    [ c  0  0  a ].                                                (1)

All entries are nonnegative. In the more usual alternating z rotation, the jump exchanges 01 and 10, and the diagonal energies are Delta/4 on parallel spins and -Delta/4 on antiparallel spins. Both representations describe the same even-ring trace. The y convention is useful because it also makes reflection positivity explicit.

[Proved] A discrete column history is s=(s_0,...,s_(2m-1)), periodic in the time index. For a site of fixed spatial parity, the active bond in layer j is known. Decorate each change s_j to s_(j+1) by L or R according to that bond. This decoration is determined by the checkerboard and the history, so no additional factor of two is inserted. Define

    K_r(s,u)=product over j congruent r mod 2
       g_[s_(j+1),u_(j+1);s_j,u_j],  r=0,1.

Here and below j+1 is taken modulo 2m. Inserting product-basis identities between all layers and regrouping the same nonnegative factors by spatial bond gives

    Z_m(L,beta)
      =Tr[(exp(-t H_even) exp(-t H_odd))^m]
      =Tr[(K_0 K_1)^(L/2)].                                      (2)

Every time layer touches every site once. Thus every nonzero vertex change is assigned to one bond, and multiplication of the column factors counts every vertex once. This is the explicit decorated-history kernel, with two alternating one-column factors. For general two-site stoquastic h the same formula applies with its own nonnegative gate, without the special entries (1).

[Proved] Entrywise nonnegativity and operator positivity must be distinguished. Realignment of (1), with row index (output-left,input-left), has the matrix

    [ a  0  0  b ]
    [ 0  c  0  0 ]
    [ 0  0  c  0 ]
    [ b  0  0  a ],

with eigenvalues a+b,a-b,c,c. For Delta>=0 these are nonnegative. Each K_r is a tensor product of these positive matrices, after a permutation of its indices. Consequently

    U_m=K_0^(1/2) K_1 K_0^(1/2)>=0,
    T_m=U_m^(1/2),
    Tr T_m^L=Z_m(L,beta)                                          (3)

on even rings. The matrix B_m=K_0 K_1 is entrywise nonnegative and has the same nonzero spectrum as U_m. The square roots in (3) need not preserve entrywise nonnegativity. We do not identify T_m with a nonnegative one-site history kernel. The safe statements are the nonnegative two-column contraction (2) and the positive operator (3). This repairs a premise in the question without changing R_L or any even-ring trace.

[Proved, inherited from C125] At fixed beta, spectral compactness reconstructs a common positive spectrum for the limiting traces at every even L>=4. It does not by itself reconstruct a pointwise continuous-history kernel or prove trace convergence at L=2. XX has an explicit positive spectrum including L=2. All physical traces themselves converge by the finite-dimensional Trotter formula. At finite m, (2) is Z_m, not the exact Z. The numerical checks distinguish contraction error from Trotter error.

[Assumed scope] For the zero-field spin-1/2 models here, vacuum-compatible lengths are the even rings, which admit the alternating rotation and carry no forced odd-ring defect. L=2 uses the convention of two identical periodic bonds. General positive-transfer chains require their own common reflection-compatible length subsequence; no classification is assumed here.

## 2. Continuous-time weights and what a Poisson bound actually proves

[Proved] In the y convention put D_b(s,u)=-Delta z(s)z(u)/4, where z(0)=1 and z(1)=-1. A configuration consists of periodic site states s_i(tau), and bond events (i,tau) that flip 00 to 11 or 11 to 00 on (i,i+1). Its differential weight is

    exp[-sum_i integral_0^beta D_b(s_i(tau),s_(i+1)(tau)) d tau]
         product over events [(1/2) d tau].                       (4)

Sum over initial site states and integrate ordered event times for each bond. Equivalently integrate unordered times and divide by the factorial on each bond. A decorated column records its piecewise constant state and the incident bond of each jump. Adjacent columns must have identical times and consistent before/after states for their shared events. The bond factor is the exponential in (4), the activity 1/2 for each common event, and the compatibility constraint. Assign the activity once to the bond, or its square root to each endpoint, but not once independently to both endpoints. The general stoquastic formula replaces 1/2 by -h_(a'b';ab)>=0 and D_b by h_(ab;ab). Events that change only one site still carry a bond label.

[Proved] Equation (4) follows by the time-ordered expansion of exp[-beta(D-W)], with W the nonnegative off-diagonal part of -H. Diagonal time intervals give the exponential action and each W insertion gives the displayed event factor. Nonnegativity permits regrouping all terms without cancellation. The d-valued state supplies finitely many marks per vertex. Locality in imaginary time supplies the action integral and event factors. None of these steps supplies a frequency-dependent spatial energy cost.

[Proved count bound, XXZ at every Delta>=0] Let N_i count all jump events on the two bonds incident to site i. Let W_i be the sum of their off-diagonal flip matrices, each of norm 1/2, so ||W_i||<=1. For fugacity z>=1 the normalized event generating function is

    F_i(z)=sum_config weight(config) z^(N_i) / Z
          =Tr exp[-beta(H-(z-1)W_i)] / Tr exp(-beta H).

The min-max principle gives lambda_j(H-(z-1)W_i)>=lambda_j(H)-(z-1)||W_i||. Summing exponentials gives

    F_i(z)<=exp[beta(z-1)].                                      (5)

In particular the mean count is at most beta. Since M=(1+max(1,Delta))/4>=1/2 for these bonds, (5) is bounded by exp[2 beta M(z-1)]. This is a Poisson upper bound on the generating function of the count. It does not assert a coupling to independent temporal clocks, or a bound on every conditional jump intensity. It is valid on the closed interacting ring, not just on an isolated vertex.

[Proved general version] For any fixed stoquastic representative the same proof gives F_i(z)<=exp[beta(z-1)||W_i||]. If J is its centered operator norm, each bond's off-diagonal part has norm at most 2J, and ||W_i||<=4J. A gauge that minimizes the bond half-spread need not preserve the chosen stoquastic expansion, so J is not silently identified with the gauge-optimized M. This issue does not occur for the XXZ representatives above. Even granting a count bound of order beta M for general chains would not close the spatial estimate.

[Proved obstruction to independent column clocks] On a short time interval t, divide the 00 to 11 gate entry by the 00 to 00 entry. From (1) the relative weight is

    c/a=tanh(t/2)=t/2+O(t^3).

Both neighboring columns jump. Two independent column events with rates bounded by a constant r have a simultaneous-bin weight of order r^2 t^2. The ratio is unbounded as t tends to zero. In continuous time the paired-jump measure is supported on equal-time diagonals in the product of neighboring event-time spaces. Independent diffuse column measures give those diagonals zero measure. Therefore no finite pointwise density majorant by such independent column measures exists, already for the explicit XX chain. The same statement holds for Heisenberg and both requested XXZ anisotropies. It persists at arbitrarily small nonzero coupling.

[Proved qualification] Independent bond clocks do generate a valid reference expansion, since one bond clock creates both endpoint changes. That construction avoids the singularity. It retains compatibility between adjacent columns and the action weight in (4). Removing these constraints gives count or volume estimates, with exponent of order beta times the number of bonds. It does not give beta/L in a normalized crossed trace. The preceding obstruction excludes independent column-clock kernel domination, not every possible bond-polymer estimate or every scalar Fock majorant.

[Proved conditional-intensity obstruction] The flip block on one active bond is a sigma_x. Its unequal-endpoint propagator is sinh(a t), a=1/2. In a bridge whose present state differs from its prescribed endpoint at remaining time r, the instantaneous jump rate is

    a cosh(a r)/sinh(a r)=a coth(a r)~1/r.                        (6)

Such bridges occur by conditioning an interval of a ring world line; the complementary interval restores the closed history. Hence uniformly bounded conditional intensities are false even though (5) holds. Average count control and history-by-history conditioning are different statements.

[Proved explanation of the missing scale] On a time mesh with 2m layers, an N-event column has up to d times binomial(2m,N) times a fixed mark-count to the power N possible records. A count estimate controls N at order beta M. It does not make this count uniform as m tends to infinity, and does not replace N by beta M/L. Obtaining the latter requires a spatial suppression of temporal structure, not just a bound on its total activity.

## 3. XX: a positive free product and a sharper scalar bound

[Proved, inherited input] Use the exact crossed spectrum in C125 PROOF section 3, rather than assuming linear dispersion for the lattice column. Restore exchange J>0, so the physical fermion velocity is v=J and M=v/2. Set b=beta v. Define

    k_n^o=asinh((2n+1)pi/b), n>=0,
    k_n^e=asinh(2n pi/b), n>=1,
    P_o^+=(product_(n>=0)(1+exp(-L k_n^o)))^4,
    P_o^-=(product_(n>=0)(1-exp(-L k_n^o)))^4,
    P_e^+=(product_(n>=1)(1+exp(-L k_n^e)))^4.

With p_c=(2pi)^(-1) integral log[2 cosh(b cos k/2)] dk and the positive sector offset delta_b of C125, the normalized trace is

    Ztilde_L=exp(-L p_c) Z_L
            =(P_o^++P_o^-)/2+2 exp(-L delta_b) P_e^+.             (7)

All spectral multiplicities in (7) are nonnegative. The first term projects four fermion modes per odd frequency to even total occupation; the second has two vacua and four fermion modes per positive even frequency. This is a free factorization specific to XX, not a consequence of count bound (5).

[Proved] The XX vertex weights satisfy a^2=b_gate^2+c^2, since Delta=0 in (1). At nonzero Delta the difference is

    a^2-b_gate^2-c^2=2 sinh(t Delta/2).

The free factorization is lost. In world-line terms, the diagonal interaction in (4) depends jointly on adjacent histories. No independent-mode identity for that weight has been established here. Thus the proof below uses the positive XX product as an extra solvable-model input. It is not advertised as a derivation of that product from Poisson counts alone.

[Proved] Since P_o^-<=P_o^+, delta_b>0, and k_n^e>=k_(n-1)^o,

    Ztilde_L<=3 P_o^+,
    log Ztilde_L<=log 3+4 sum_(n>=0) f((2n+1)pi/b),
    f(x)=log(1+exp[-L asinh x]).                                 (8)

The function f is convex on x>=0. Indeed

    f'(x)=-L/[sqrt(1+x^2)(1+exp(L asinh x))],

whose denominator increases. Convex midpoint quadrature on each interval [2n pi/b,2(n+1)pi/b] gives

    sum_(n>=0) f((2n+1)pi/b)
        <=b/(2pi) integral_0^infinity f(x) dx.                    (9)

The integral converges for L>1. With x=sinh u it is

    I_L=integral_0^infinity cosh(u) log(1+exp(-L u)) du
       =L integral_0^infinity sinh(u)/(1+exp(L u)) du.

The boundary term in integration by parts vanishes for L>1. Put a=1/L and t=exp(-L u). Then

    I_L=(1/2) integral_0^1 (t^(-a)-t^a)/(1+t) dt
       =(1/2)[integral_0^infinity t^(a-1)/(1+t) dt-1/a]
       =pi/[2 sin(pi/L)]-L/2.                                   (10)

The middle equality uses t^a/(1+t)=t^(a-1)-t^(a-1)/(1+t) and inversion on the interval (1,infinity). The last integral is the Euler beta integral for 0<a<1. Combining (8)-(10),

    log Ztilde_L<=log 3+b[csc(pi/L)-L/pi],  L>1.                  (11)

For physical traces use even L>=2. This is the main new XX estimate.

[Proved lower denominator and ratio bound] If the largest transfer eigenvalue is lambda_0, set q_j=lambda_j/lambda_0 in [0,1], with at least one q_j=1. Then

    1<=Ztilde_(2L)<=Ztilde_L,
    log Ztilde_L<=R_L<=2 log Ztilde_L,
    Z_(2L)>=lambda_0^(2L).                                      (12)

For XX lambda_0=exp(p_c). Equations (11)-(12) provide both the normalized numerator bound and the required lower bound on the denominator. In particular

    R_L<=2 log 3+2 beta v[csc(pi/L)-L/pi].                        (13)

No extensive vacuum normalization is lost in this step.

[Proved explicit constants] The positive-coefficient expansion of x/sin x implies that L[csc(pi/L)-L/pi] decreases for L>1. It is 2-4/pi at L=2 and 4sqrt(2)-16/pi at L=4. Therefore

    R_L<=2 log 3+(8-16/pi) M beta/L,        even L>=2,
    R_L<=2 log 3+(16sqrt(2)-64/pi) M beta/L, even L>=4.             (14)

At fixed L, beta tending to infinity gives the same constants for L(E_(2L)-2E_L)/M. The second coefficient is about 2.256, compared with C119's 10pi/3. The additive constant 2 log 3 is slightly larger than C119's 3 log 2. This is an improvement of the displayed XX bound, not of C119's theorem for all its other free families.

[Proved failure of unchanged XX constants for Heisenberg] The XX estimate (11) with v=2M=1 cannot extend unchanged to Heisenberg. Its exact four-site ground energy is E_4=-2, while the inherited Bethe bulk energy is e_0=1/4-log 2. The largest crossed eigenvalue has logarithm equal to the infinite-chain pressure. Thus, at fixed L=4,

    lim_(beta->infinity) log Ztilde_4/beta
       =4e_0-E_4=3-4 log 2=0.2274112778... .

The proposed XX right side has slope sqrt(2)-4/pi=0.1409740176... . The additive log 3 cannot repair that strict slope discrepancy. This is an actual lattice failure of the unchanged XX majorant, not a failure of all choices of S, V', or K. It also agrees with Heisenberg's inherited A=pi/2, versus XX's A=1 at the same M=1/2.

[Proved literal linear Fock majorant at the trace level] A product with linear frequencies can also majorize the entire trace, despite failing mode by mode. Let x=b/L, q=1-exp(-1), and c=-log q. For every integer k<=x,

    -log(1-exp(-k/x))>=c.

Consequently

    log product_(k>=1)(1-exp(-k/x))^(-2)>=2c floor(x)>=2cx-2c.

Since 2c>2-4/pi, (11) gives, for L>=2,

    Ztilde_L<=exp(K) product_(k>=1)
                 [1-exp(-2pi k L/(V' beta))]^(-S),
    K=log 3+2c, S=2, V'=2pi v=4pi M.                            (15)

This is an aggregate trace comparison obtained after (11). It is not an injection of lattice modes into same-frequency linear modes, and is not a positive-kernel inequality. The standard decreasing-function integral bound gives log of the product <=S pi V' beta/(12L). Thus (15) has the requested form with explicit constants, although it is weaker than (11).

[Proved consistency checks] As L tends to infinity, csc(pi/L)-L/pi=pi/(6L)+O(L^(-3)). The coefficient in (11) matches the exact XX normalized crossed pressure, with A=v, and does not discard slow velocities. At beta tending to zero, Ztilde_L tends to one, while the right side remains finite. At L=2, (10) gives I_2=pi/2-1, finite. Real L<=1 is excluded, as the XX product diverges there. Independent XX stacks multiply (7), hence add (11)-(14), with sum v_a=2M and additive constant 2n log 3 for d=2^n.

## 4. Frequencies, reflection, and the interacting missing input

[Proved, inherited lattice obstruction] C125 gives locally detected XX states at Matsubara frequency omega_k=2pi k/beta with gap

    e_k=asinh(pi/(beta v))+asinh((2k-1)pi/(beta v)).               (16)

At fixed beta, e_k grows logarithmically with k, so e_k>=omega_k/V' fails for every finite V'. This is a counterexample inside the explicit nearest-neighbor, stoquastic, reflection-positive class. Neither a coupling of histories nor a reflection estimate can prove the false all-frequency linear statement. It does not refute the aggregate trace comparison (15).

[Proved, inherited corrected decay] C125's Laplace-transform Lieb-Robinson argument gives V=4eM and the asymptotic spatial exponent

    alpha(omega)=omega/V,           0<omega<=V,
                 log(omega/(4M)),   omega>=V.                    (17)

At the matching point a polynomial prefactor is allowed. For 0<omega<V the explicit two-point estimate is

    |G_AB(i omega,x)|<=2||A||||B|| V/[omega(V-omega)]
                          times exp(-omega x/V).

Only a positive reflected spectral measure allows this to exclude gaps below alpha in the detected sector. General off-diagonal correlators can cancel, and unseen sectors are not bounded by this argument. Applying a Fourier transform to nonnegative histories introduces oscillatory factors; entrywise nonnegativity alone is not a frequency-by-frequency comparison. This run supplies no independent world-line coupling proof of (17).

[Proved limitation of reflection] Realignment positivity in section 1 proves (12), spatial log-convexity, and Schwarz inequalities for reflected insertions. These are useful for the denominator and for positive spectral measures. They do not assign an additive gap to a product of several insertions. A chessboard application would still need to bound the resulting repeated-defect partition functions, including their vacuum normalization. Taking that estimate to be exp[O(beta M/L)] would insert the missing target. The classical rank argument works when histories collapse to finitely many separator labels, as in C124's commuting proof. Noncommuting XX already has infinitely many crossed modes despite its finite set of local marks.

[Assumed sufficient additional lattice input] Here is a precise replacement for linear-frequency domination. Suppose, uniformly in a convergent Trotter approximation, one can prove

    Ztilde_L<=exp(K_0+c_0 beta M/L)
       product_(k>=1)[1-exp(-L alpha(2pi k/beta))]^(-S(d)),        (18)

for even L>=2, with alpha from (17). The prefactor includes all zero-frequency sectors. The new content is a bound on the number of independent modes and their many-insertion partition function, with additive decay costs after grouping events into interacting clusters. A bound on one-column counts, a bound on single-insertion spans, or a bound on a two-point function does not establish (18). Interactions can mix clusters and produce low states not described by independent insertion costs. No theorem excluding that effect quantitatively is proved here.

[Proved conditional closure] Put f_L(omega)=-log(1-exp[-L alpha(omega)]). It is decreasing. Hence

    sum_(k>=1) f_L(2pi k/beta)
       <=beta/(2pi) integral_0^infinity f_L(omega) d omega.

The low-frequency integral is at most V pi^2/(6L). On the high-frequency interval expand the positive logarithmic series and integrate each term:

    integral_V^infinity f_L(omega) d omega
       =V sum_(r>=1) exp(-Lr)/[r(Lr-1)]
       <=V/(L-1) [-log(1-exp(-L))]
       <=V/(L-1) /(exp(L)-1)
       <=V/L times 2/(exp(2)-1), L>=2.

Therefore (18) implies

    log Ztilde_L<=K_0+
       {c_0 M+S(d)V/(2pi)[pi^2/6+2/(exp(2)-1)]} beta/L.          (19)

Combining with (12) proves the requested ratio bound. Under the stated CFT limit it also bounds A. Thus logarithmic ultraviolet decay is sufficient; forcing linear decay in the ultraviolet is unnecessary. The missing step is the interacting many-state estimate (18), including zero-frequency weight, not the convergence of its high-frequency tail.

[Open conclusion] We have proved a Poisson count bound for the requested XXZ models and a Fock trace bound for XX, but not (18) for Heisenberg, XXZ at Delta=0.5 or 2, or general stoquastic positive-transfer chains. No actual lattice chain violating every choice of K(d), S(d), V'=O(M) has been constructed. The explicit lattice failures concern independent column-clock domination and linear modewise decay. Calling either of them a counterexample to the scalar C(d) target would be incorrect. Equation (18) is a sufficient missing input, not a claim of a uniquely minimal axiom; the direct effective-rank estimate is a weaker sufficient scalar input.
