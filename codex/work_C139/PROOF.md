# C139 proof: the repulsive U(1) class, 0 <= Delta < 1

2026-11-07. Run by a Claude Opus agent, not Astra. J=1 in all analytic sections (restore J by homogeneity). Notation follows C134 and C136: Delta=cos gamma, now 0<gamma<=pi/2, and

    p(x)=2atan[tanh(x)cot(gamma/2)],   theta(x)=2atan[tanh(x)cot(gamma)],
    a1=p'/(2pi)=sin(gamma)/[pi(cosh(2x)-Delta)],
    K=theta'/(2pi)=sin(2gamma)/[2pi(sinh(x)^2+sin(gamma)^2)],
    e(x)=-a a1(x),  a=pi sin(gamma),  b=pi/gamma,
    q=||K||_1=1-2gamma/pi,  B=1/(1-q)=pi/(2gamma),
    r0(x)=1/[2gamma cosh(bx)]   (C134's rho),
    I_j=j-(n+1)/2,  nu=n/L,  1<=n<=L/2.

The Bethe equations are L p(x_j)-sum_k theta(x_j-x_k)=2pi I_j, the candidate energy is Delta L/4+sum_j e(x_j), and

    mu=(1/L)sum_j delta_(x_j),   sigma_L=a1-K*mu,
    z_L(x)=p(x)/(2pi)-(1/L)sum_k theta(x-x_k)/(2pi),  z_L'=sigma_L,
    z_L(+-infinity)=+-c,  2c-nu=(1-gamma/pi)(1-2nu).

For 0<gamma<pi/2 the kernel K is strictly positive, even, decreasing on x>0, and its transform Khat(w)=sinh((pi-2gamma)w/2)/sinh(pi w/2) is strictly positive. At gamma=pi/2 everything reduces to free fermions (K=0).

Status tags: [Proved], [Proved, conditional on X], [Numerical], [Assumed], [Literature].

## 0. Results in one place

* Theorem A [Proved, conditional on the identification hypothesis (I) of Section 2]. For 0<delta<=1, Delta in [0,1-delta], every L>=2 and 0<=n<=L/2,

      |F_L(n)-L e_bulk(Delta,n/L)| <= K_rep(gamma)/L <= (1086/delta)/L,

  and on the repulsive U(1) class (J_z/J in [0,1-delta], any field, any DM angle, every L>=2)

      L(E_(2L)-2E_L) <= C(delta) M,   C(delta)=5400/delta.

* Hypothesis (I) is [Literature] for even L (Yang and Yang 1966, as restated by Kozlowski 2015; existence of a ground-state real solution continuous in Delta, symmetry not stated) and [Proved, conditional on (N)] for all L, where (N) is nonsingularity of the Yang-Yang Hessian along the continuation branch. (N) follows from sigma_L(x_j)>=0 at the roots. [Numerical] sigma_L>0 on the whole real line in every case computed.
* Theorem C [Proved, conditional on (M)]: if the ground-state counting function is monotone on its Fermi interval (sigma_L>=0 on [-Q_L,Q_L]), all constants become uniform in gamma in (0,pi/2], and L(E_(2L)-2E_L) <= (2215+pi^2/4) M on J_z/J in [0,1], Heisenberg included by finite-matrix continuity.
* The single remaining estimate is (M) (Section 8.4). Everything else in the repulsive class, including the continuum dressing at a finite Fermi boundary, is proved here unconditionally.

## 1. Task 1: where C134 and C136 use the sign of the kernel

Each line gives the place, what the negative sign gave, and the replacement for K>0 (or the failure).

C134.

1. Section 2, Yang-Yang Hessian. theta'<=0 made Y strictly convex, giving existence, uniqueness, symmetry and an analytic branch. For K>0 this fails (Kozlowski 1508.05741, Section 2: the functional "is not convex for Delta>0"). Replacement: coercivity survives unchanged (Section 2.1), so a strictly ordered real solution exists; the Hessian equals 2pi L diag(sigma_L(x_j)) plus a strictly positive Gram matrix (Lemma 2.2), so it is positive whenever sigma_L>=0 at the roots. Uniqueness and identification become hypothesis (N)/(I). This is the one place where no unconditional replacement was found.
2. Section 2, coordinate Bethe vector, analytic characteristic polynomial, Perron-Frobenius. Sign-free. For distinct real rapidities the two contact factors 1+z_j z_k+2Delta z_(j or k) have equal modulus, so they vanish only together, which forces z_j=z_k.
3. Section 3, Neumann inverse norm <= B. Sign-free (two-sided series, ||K||_1=q<1). The full-line formulas r0=(I+K*)^(-1)a1 and epsilon=-a r0 are Fourier identities valid for all 0<gamma<pi. Better, for K>0 the true inverse is I-J* with J>=0 and ||J||_1<1/2, so ||(I+K*)^(-1)||<=3/2 uniformly (Section 3). The Neumann value B is lossy here.
4. Section 4, global positivity sigma=a1-K*mu>0 and the midpoint representation in a global counting coordinate. Not provable from the sign; numerically true. Replacement: Lemma 4.1, a cumulative root bound ||G||<=B/(2L) that uses only the ordering of the roots (half the C134 constant B/L).
5. Section 4, derivative bounds ||K'||_1=2|K(0)| and |K''|<=|K|(8/s^4+4/s^2). Sign-free (|K| is even unimodal). Sharper here: ||K''||_1=4 max|K'| (Lemma 5.0).
6. Section 5, exact dressed quadrature identity. Sign-free.
7. Section 6, change of variables t=z(x) on the whole counting interval [-c,c]. Needs global monotonicity. Replacement: only on the core where sigma_L>=r/2 is guaranteed, with first-order bounds outside (Section 5.3).
8. Section 7, odd-L outer strips discarded with the favourable sign f<=0. Uses the global coordinate. Replaced by the finite-boundary identity at Q_L, which treats both parities alike.
9. Sections 7-8, tail masses from (4.2) and the small-size spectral branch. Sign-free.
10. Section 11 (Delta=1/2). Independent of the kernel; reused as a check.

C136 (all with T=-K>=0 there).

11. (2.1) positivity of (I-T_Q)^(-1). Fails: (I+K_Q)^(-1) has an alternating series. Replacement: hole duality s_Q=r0+J*(1_out s_Q) with J>=0 (Lemma 3.2), a positive series on the complement.
12. (2.2) a1<=r_Q<=r_infinity and 1<=Z_Q<=B, monotone in |x|. The order reverses: r0<=r_Q<=(1+4/pi)r0, 1/(1+q)<=Z_Q<=1 (Lemma 3.3). Monotonicity in |x| is not needed; Var(Z_Q)<=q_J replaces it (Lemma 3.4).
13. (2.3)-(2.4) m'(Q)=2r_Q(Q)Z_Q(Q)>0, and 1/2-m(Q). Positivity now from Lemma 3.3; 1/2-m=(1-q_J)h_Q with h_Q<=int_(|x|>Q) r0/(1-q_J) (Lemma 3.3).
14. (3.2)-(3.3) midpoint counting with monotone z_L. Replaced by Lemma 4.1 at filling nu, which needs only monotonicity of the continuum counting function (r_nu>=0).
15. (4.1)-(4.2) Q_L, the sawtooth |Cum R_L|<=1/(2L), Var(Z)<=2(B-1). Replaced: Q_L is the first crossing z_L=nu/2 beyond x_n, |Cum R_L|<=3B/(2L) (Lemma 4.3), Var(Z)<=q_J.
16. (5.1)-(5.2) epsilon_Q'>0 from the odd reduction T(x-y)-T(x+y)>=0, and |epsilon_Q'|<=2a r_Q from |a1'|<=2a1. Both fail for K>0: the odd reduction of -K_Q is negative, and |a1'|/a1=2sinh(2x)/(cosh(2x)-Delta) exceeds 2 near x=0 when Delta>0. Replacement: the dual odd reduction J(x-y)-J(x+y)>=0 gives epsilon_Q'>=0 on (0,Q] and |epsilon_Q'|<=ab r_Q (Lemma 3.5).
17. (5.3) 0<=e_bulk''<=a. Replaced by 0<=e_bulk''<=ab(1+q)^2/2 (Lemma 3.5).
18. (6.1)-(6.4) positive resolvent majorant, kappa lower bound, |r_Q'|<=H1 r_Q, Var(v_Q). All rebuilt from duality, with kappa=1, |r_Q'|<=b r_Q, |r_Q''|<=b^2 r_Q and Var(v_Q)<=(2+40/pi)ab (Lemma 3.5). This part is simpler than in C136.
19. (6.4) restricted densities ordered in Q. Holds for the extended density s_Q, which decreases in Q (Lemma 3.6).
20. (8.4)-(8.5) high-field region via the Delta=-1 operator inequality. Attractive-specific, not needed.
21. Section 9 and C129 reductions (flux, field, DM, M>=M_0=J/2 for |Delta|<=1, seam gluing). Sign-free.

So the sign enters essentially only through global monotonicity of the finite counting function (items 1, 4, 7, 14, 15). Everything on the continuum side has a replacement, and on the finite side the replacements cost the factor B=pi/(2gamma), which is the only source of non-uniformity as Delta->1.

## 2. Existence and identification of the ground-state roots

### 2.1 Existence and ordering [Proved]

For 0<gamma<pi the Yang-Yang functional Y(x)=L sum P(x_j)-sum_(j<k) T(x_j-x_k)-2pi sum I_j x_j (P'=p, T'=theta) is coercive. The C134 recession computation does not use the sign of theta: moving the top k coordinates right gives slope k[L(pi-gamma)-(pi-2gamma)(n-k)-pi(n-k)]=k(pi-gamma)(L-2n+2k)>0, with the sign-adapted ray decomposition of the C134 referee, and sorting cannot increase Y. The same bound applied to <grad Y(x),x>, using x p(x)>=(pi-gamma)|x|-C_1(gamma) and (x_j-x_k)theta(x_j-x_k)<=|pi-2gamma||x_j-x_k|, puts every critical point in a ball whose radius is bounded for gamma in compact subsets of (0,pi). A global minimizer exists and is sorted. Every critical point has distinct coordinates, since equal x_j=x_k would give equal left sides with different right sides 2pi I_j. Hence along any continuous family of critical points the strict order is preserved.

### 2.2 The Hessian [Proved]

    Y''=2pi L diag(sigma_L(x_j)) + Theta,   Theta_jk=theta'(x_j-x_k).

Indeed Y''_jj=Lp'(x_j)-sum_(k!=j)theta'(x_j-x_k)=2pi L sigma_L(x_j)+theta'(0). Since theta'=2pi K has strictly positive Fourier transform, Bochner gives u^T Theta u=int thetahat'(w)|sum u_j e^(iwx_j)|^2 dw/(2pi)>0 for distinct x_j and u!=0. Thus Y''>0 whenever sigma_L(x_j)>=0 for all j.

### 2.3 Identification [Proved, conditional on (N)]

(N): along the branch of critical points continued from the free solution at gamma=pi/2, det Y''!=0 for every gamma in [gamma_0,pi/2].

Under (N) the implicit function theorem gives a real-analytic branch on a maximal interval; by the a priori ball of 2.1 the branch has limit points, which are ordered critical points; by (N) the branch continues through them, so it covers [gamma_0,pi/2]. Local uniqueness in the implicit function theorem and the reflection x->(-x_n,...,-x_1) give x_j=-x_(n+1-j). The C134 argument (coordinate Bethe vector near the free point, identity theorem for det[K_L(cos gamma)|_n-E(gamma)], Perron-Frobenius simplicity and the open-closed argument) then gives E(gamma)=F_L(n) at both parities. By 2.2, (N) holds if sigma_L(x_j)>=0 at every root along the branch.

(I), the hypothesis actually used below: for each L,n, gamma there is a real, strictly ordered, reflection-symmetric solution with the ground-state I_j whose Bethe energy is F_L(n). (N) implies (I).

[Literature] Kozlowski (arXiv:1508.05741v2, introduction) states that Yang and Yang proved that the ground state for Delta>=-1 comes from a real solution with these quantum numbers, continuous in Delta, and that they did not prove uniqueness for Delta>0. Kozlowski assumes even length. The YY original was not accessible here. This covers (I) at even L except the symmetry clause; it does not cover odd L at the minimizing flux. The spec's statement that Yang-Yang convexity identifies the ground state on all of (-1,1) is therefore not supported in that form.

[Numerical] Bethe energies equal ED sector minima at minimizing flux for all 1<=n<=L/2, L=2,...,12, Delta in {0.1,0.5,0.9,0.99,0.999} (max error 6.2e-14, ed_check.json). sigma_L>0 everywhere and at all roots in the 572 cases of check_finite.py, the 144 cases of sigma_sign.py and the 65 cases of the half-filled scan heis_scan.py (doubled rings to 8192 sites).

Note. The energy estimates in Sections 4-6 use only (I): any ordered symmetric real solution, unique or not.

## 3. The hole kernel and the continuum at a finite Fermi boundary [Proved]

### Lemma 3.1 (hole kernel)

Let J be the kernel of (I+K*)^(-1)K*, so Jhat=Khat/(1+Khat). With B'=(pi-gamma)/2, C'=gamma/2 (so B'-C'=(pi-2gamma)/2 and B'+C'=pi/2),

    Jhat(w)=[sinh((B'-C')w)/sinh(B'w)] [1/(2cosh(C'w))].

The first factor is the transform of k1(x)=(1/(2B'))sin(phi)/(cosh(pi x/B')+cos phi), phi=pi(B'-C')/B' in (0,pi), positive, even, decreasing in |x|, with ||k1||_1=(pi-2gamma)/(pi-gamma). The second is the transform of r0. Hence

    J=k1*r0>=0, even, nonincreasing in |x| (convolution of symmetric unimodal densities),
    ||J||_1=q_J=(pi-2gamma)/(2(pi-gamma))<1/2,   J(0)<=||k1||_1 ||r0||_inf<=1/(2gamma),
    |J'|<=bJ,  |J''|<=b^2 J   (from |r0'|<=b r0, |r0''|<=b^2 r0).

Also (I+K*)^(-1)=I-J* on bounded functions, of norm <=1+q_J<=3/2 uniformly in gamma. [Numerical check: jkernel.py, minimum of J is zero to 1e-15, monotone, norm matches.]

### Lemma 3.2 (duality)

For Q in (0,infinity] let r=r_Q solve r+K_Q r=a1 on [-Q,Q] and extend s_Q=a1-K*(1_Q r) to R (so s_Q=r on [-Q,Q]). Then

    s_Q=r0+J*(1_out s_Q),   out={|x|>Q}.

Proof: s+K*s=a1+K*(1_out s); apply (I+K*)^(-1). Conversely the unique bounded solution of the outside equation (Neumann, norm q_J) reproduces r, by injectivity of I+K* on bounded functions. The same algebra applied to 1 and to e-h gives Ztilde=1/(1+q)+J*(1_out Ztilde) and epsilontilde=epsilon_0+J*(1_out epsilontilde), epsilon_0=-a r0-h/(1+q).

### Lemma 3.3 (density bounds)

All Neumann terms on the complement are nonnegative, so on R

    r0<=s_Q,   1/(1+q)<=Ztilde<=1,
    h_Q:=int_out s_Q <= int_out r0/(1-q_J) <= (2/(pi(1-q_J)))e^(-bQ),
    m(Q)=1/2-(1-q_J)h_Q.

The last identity is the integral of s_Q=r0+J*(1_out s_Q) over R together with int s_Q=m+h_Q. Inside, U1:=J*(1_out s_Q) obeys U1(x)<=J(0)h_Q, and r0(x)>=e^(-b|x|)/(2gamma), so

    0<=U1/r0<=omega(x):=(4/pi)e^(-b(Q-|x|)),   r0<=r_Q<=(1+4/pi)r0  on [-Q,Q].

Also |r_Q'|<=b r_Q, |r_Q''|<=b^2 r_Q (differentiate r_Q=r0+U1 and use Lemma 3.1), and m'(Q)=2r_Q(Q)Z_Q(Q)>0, so m is a bijection from [0,infinity] to [0,1/2].

### Lemma 3.4 (variation of Z)

On [-Q,Q], Z'=J'*(1_out Ztilde), and

    Var_[-Q,Q](Z_Q) <= ||Ztilde|| int_(-Q)^Q int_out |J'(x-y)| dy dx <= 2 int_0^inf u|J'(u)|du = q_J,

using unimodality of J. [Numerical: VarZ/q_J<=0.975, check_continuum.json.] Also the restricted inverse is uniformly bounded: if (I+K_Q)f=g, then F=g1_Q-K*(1_Q f) solves F=(I+K*)^(-1)(g1_Q)+J*(1_out F), and

    ||(I+K_Q)^(-1)||_(inf->inf) <= (1+q_J)/(1-q_J) <= 3.            (3.1)

### Lemma 3.5 (dressed energy)

Let epsilon_Q=(I+K_Q)^(-1)(e-h), h=h(Q)=-a r_Q(Q)/Z_Q(Q), so epsilon_Q(+-Q)=0. The C136 identities dE/dQ=2r_Q(Q)g(Q), e_bulk'(nu)=h(Q), h'(Q)=epsilon_Q'(Q)/Z_Q(Q) and

    e_bulk''(nu)=epsilon_Q'(Q)/[2r_Q(Q)Z_Q(Q)^2]

hold for either sign (they use only self-adjointness and epsilon_Q(+-Q)=0). Differentiating the dual equation, the boundary terms vanish because epsilontilde(+-Q)=0, so epsilontilde'=epsilon_0'+J*(1_out epsilontilde'), epsilon_0'=-a r0'. Termwise,

    |epsilontilde'|<=ab s_Q,   |epsilontilde''|<=ab^2 s_Q,
    epsilon_Q'>=0 on (0,Q]   (odd reduction J(x-y)-J(x+y)>=0 for x,y>0),
    0<=e_bulk''<=ab(1+q)^2/2=:a2<=2ab<=2pi^2.

Write v_Q=epsilon_Q'/r_Q, ell=r0'/r0=-b tanh(bx), U2=J*(1_out epsilontilde'). Then v_Q=-a ell+w, w=(U2+a ell U1)/r_Q, |w|<=2ab U1/r_Q, and |w'|<=5ab^2 U1/r_Q<=5ab^2 omega. Hence

    |v_Q|<=ab,   Var_[-Q,Q](v_Q) <= 2ab+40ab/pi =: Vv,                (3.2)

uniformly in Q and in gamma (ab=pi^2 sin(gamma)/gamma<=pi^2). [Numerical: Var(v)/Vv<=0.14, e_bulk''/a2<=0.94, r_Q/r0 in [1,1.19], check_continuum.json.]

### Lemma 3.6 (monotonicity in Q)

d s_Q/dQ=-s_Q(Q)(I-J1_out)^(-1)[J(.-Q)+J(.+Q)]<=0, with |d s_Q/dQ|<=(2/gamma)r_Q(Q)<=(2/gamma)(1+4/pi)r0(Q) and the same bound times b for the x-derivative.

## 4. Finite-size root estimates without the kernel sign

Assume (I). Let nu=n/L, Q_nu the continuum Fermi point with m(Q_nu)=nu (Q_nu=infinity at even half filling), s_nu=s_(Q_nu), z_nu(x)=-c+int_(-infinity)^x s_nu (same c as z_L, since c depends only on nu), C_nu=P_nu(z_nu) (projection onto [-nu/2,nu/2]), C_mu(x)=mu((-infinity,x])-nu/2, D=C_mu-C_nu.

### Lemma 4.1 (cumulative bound) [Proved, conditional on (I)]

    ||D||_inf <= B/(2L),   ||z_L-z_nu||_inf <= qB/(2L).

Proof. z_L-z_nu=-K*D (integration by parts; D vanishes at +-infinity because both measures have mass nu). At a root, C_mu(x_j^+)=j/L-nu/2 and z_L(x_j)=I_j/L in [-nu/2,nu/2], so z_L(x_j)=P_nu(z_L(x_j)) and

    D(x_j^+)=1/(2L)+P_nu(z_L(x_j))-P_nu(z_nu(x_j)) in [1/(2L)-q||D||, 1/(2L)+q||D||],
    D(x_j^-)=D(x_j^+)-1/L.

Between roots C_mu is constant and C_nu is nondecreasing (Lemma 3.3: s_nu>=r0>0), so D is nonincreasing; D(+-infinity)=0. Hence ||D||<=1/(2L)+q||D||. Only the ordering of the roots was used. [Numerical: L||D|| is between 0.50 and 0.522 in all 572 cases, against B/2 up to 17.6, finite_summary.json.]

### Lemma 4.2 (density) [Proved, conditional on (I)]

    sup|sigma_L-s_nu| <= ||K'||_1||D|| <= D_0/L,   D_0=B cot(gamma)/pi,
    sup|sigma_L'-s_nu'| <= ||K''||_1||D|| <= D1_0/L,  D1_0=2M_1 B,

with M_1 from Lemma 5.0. [Numerical: L sup|sigma_L-s_nu|/D_0<=0.131.]

### Lemma 4.3 (Fermi boundary) [Proved, conditional on (I)]

Let Q_L be the first x>=x_n with z_L(x)=nu/2 (it exists, z_L(x_n)=nu/2-1/(2L), z_L(infinity)=c>=nu/2; Q_L=infinity at even half filling); by symmetry z_L(-Q_L)=-nu/2. Let R=mu-1_[-Q_L,Q_L]sigma_L dx, total mass 0, Cum R=C_mu-z_L on [-Q_L,Q_L]. Then

    sup_[-Q_L,Q_L]|Cum R| <= 3B/(2L),    |nu-m(Q_L)| <= q_J 3B/(2L) <= 3B/(4L).

Proof: |C_mu-z_L|<=|D|+|z_nu-z_L|+|C_nu-z_nu|, and C_nu-z_nu vanishes on [-Q_nu,Q_nu] while on [Q_nu,Q_L] it is at most z_nu(Q_L)-nu/2<=z_L(Q_L)+qB/(2L)-nu/2=qB/(2L). The mismatch is nu-m=int Z dR=-int Cum R dZ (C136 (4.1), sign-free), bounded by Var(Z) sup|Cum R|.

Consequently, by Lemma 3.6 and m'(Q)>=r0(Q), on [-Q_L,Q_L]

    |sigma_L-r_(Q_L)|<=D'/L,  D'=D_0+(2/gamma)(1+4/pi)(3B/4),
    |sigma_L'-r_(Q_L)'|<=D1'/L,  D1'=D1_0+b(2/gamma)(1+4/pi)(3B/4).

## 5. Quadrature and Theorem A

### Lemma 5.0 (kernel derivatives) [Proved]

With s=sin gamma and t=|sinh x|, |K'|=(sin 2gamma/pi) t sqrt(1+t^2)/(t^2+s^2)^2, and t/(t^2+s^2)^2<=9/(16 sqrt3 s^3), t^2/(t^2+s^2)^2<=1/(4s^2) give

    max|K'| <= M_1 := (2/pi)(0.325/s^2+0.25/s).

K'' has the sign of 2u^2+(3-2s^2)u-s^2 (u=sinh^2 x), which has one positive root, so K' is unimodal on (0,infinity) and ||K''||_1=4 max|K'|<=4M_1.

### 5.1 Exact identity [Proved, sign-free]

With epsilon=epsilon_(Q_L), all roots in [-Q_L,Q_L] (C136 (4.3)):

    F_L(n)-L[E(Q_L)+h(nu-m)] = L int epsilon dR.                     (5.1)

At even half filling Q_L=infinity, h=0, epsilon=-a r0 (C134 (5.1)).

### 5.2 Core [Proved, conditional on (I)]

Let A>=0 solve r0(A)=2D'/L (possible when L>=4gamma D'). On [-A,A] intersected with [-Q_L,Q_L], r_(Q_L)>=r0>=2D'/L, so sigma_L>=r_(Q_L)/2>0 and z_L is increasing. Let x_(j1),...,x_(j2) be the roots there. They are consecutive, z_L(x_j)=I_j/L, and with f=epsilon o z_L^(-1) the trapezoid Peano kernel (sup h^2/8) gives

    |sum_(j1..j2) epsilon(x_j) - [epsilon(x_j1)+epsilon(x_j2)]/2 - L int_(x_j1)^(x_j2) epsilon sigma_L| <= Var(f')/(8L).

Here f'=v(r/sigma_L), v=v_(Q_L), so Var(f')<=||v||Var(r/sigma_L)+2Var(v). With |(r/sigma)'|<=4(bD'+D1')/(L r) and int_(-A)^A dx/(L r0)<=1/(bD'),

    Var(f') <= ab[4+4D1'/(bD')] + 2Vv = 8ab+80ab/pi+4aD1'/D' =: V_f.

### 5.3 Ends [Proved, conditional on (I)]

Integrating by parts on (x_j2,Q_L], with Cum R(x_j2^+)=1/(2L) exactly and Cum R(Q_L)=0, the right end plus the half weight equals -L int_(x_j2)^(Q_L) epsilon'(C_mu-z_L)dx. No monotonicity is used here. Its size is at most L(3B/(2L)) ab int_(x_j2)^(Q_L) r_(Q_L), and

    int_(x_j2)^(Q_L) r0 <= [3B/2+4(1+4/pi)D'/b]/L,

from Lemma 4.1 at x_j2 and at A, int_A^inf r0<=2r0(A)/b, and the bound z_nu(Q_L)-nu/2<=qB/(2L). Same on the left. If no root lies in the core, then (for L>=5B and A>=1/b) n<=B, and the same first-order bound over all of [-Q_L,Q_L] gives a smaller term.

### 5.4 Theorem A [Proved, conditional on (I)]

With c_r=1+4/pi, u=1/gamma,

    K_large = V_f/8 + 3B ab c_r(3B/2+4c_r D'/b) + a2(3B/4)^2/2,
    L_cut = max(4gamma D' cosh(1), 5B),
    K_rep = max((L_cut+1)^2, K_large).

For L>=L_cut, (5.1), 5.2, 5.3 and the curvature bound |e_bulk(nu)-e_bulk(m)-h(m)(nu-m)|<=a2(nu-m)^2/2 give |F_L(n)-L e_bulk(nu)|<=K_large/L. For L<L_cut, |F_L-L e_bulk|<=L (bond spectral diameter 1, e_bulk in the bond interval as in C136). So

    |F_L(n)-L e_bulk(Delta,n/L)| <= K_rep(gamma)/L      (all L>=2, 0<=n<=L/2).

Majorization (gamma in [gamma_0,pi/2], u=1/gamma>=2/pi, sin gamma>=2gamma/pi, cot gamma<=u, ab<=pi^2, a<=pi/u): D'<=5.856u^2, D'>=5.356u^2, D1'<=19.67u^3, V_f<=376.5, K_large<=47.1+2055.3u^2, L_cut<=36.15u, hence K_rep<=2171.5u^2. Since gamma_0=arccos(1-delta)>=sqrt(2delta),

    K_rep <= 1086/delta  on Delta in [0,1-delta].

Pointwise values (constants.py): K_rep=518 at Delta=0, 1537 at 0.5, 9740 at 0.9, 1.02e5 at 0.99; sup_delta delta sup K_rep=1027.

## 6. Assembly on the repulsive U(1) class [Proved, conditional on (I)]

Doubling with matched sectors (2L,2n), and seam gluing for L<L_cut (sector defect <=4M_0, C136 Section 9), give

    L[F_(2L)(2n)-2F_L(n)]/M_0 <= max(4L_cut, 5K_large) <= 5373/delta.

The C129/C136 reductions (DM angle to flux, field cancellation in matched sectors, flux minimization cost pi^2 M(B)/4, M(B)>=M_0=J/2 for |Delta|<=1) are sign-free, so for every bond with J_z/J in [0,1-delta], every field, DM angle and L>=2,

    L(E_(2L)-2E_L) <= [5373/delta+pi^2/4] M <= (5400/delta) M.

Together with C136 (9.2) this covers J_z/J in [-1+delta,1-delta] with C=6*10^15/delta^(25/3), still conditional on (I) for the positive part.

## 7. Consistency checks

* Delta=0 [Proved]: K=0, J=0, q=q_J=0, B=1, r_Q=a1, Z=1, and (5.1) is the free midpoint identity. K_rep(0)=518 versus the exact sector constant 8-4sqrt2 (C129): consistent, coarse. C136's own constant at its boundary delta=1 (which is the single point Delta=0) evaluates to 1.5*10^12; the two constants agree on validity, not size. The attractive and repulsive constants are both finite at Delta=0, so the two compact theorems glue across Delta=0 without a gap.
* Delta=1/2 [Numerical, Proved in C134]: odd-L Bethe energies equal -3L/8 to 4e-15 (L to 257), and the even-L twisted ED energy at total twist 2pi/3 equals -3L/8 to 2e-14 (L<=12). Theorem A at 1/2 gives 5K_rep/2=3842 J against the exact 2pi^2/9 J. Consistent.
* Casimir [Numerical]: L(F_L-Le_bulk) at half filling, L=1024 or 2048, equals -pi v/6, v=pi sin(gamma)/(2gamma), to 2e-6 relative at Delta=0.1, 0.5 and to 2e-4 to 5e-4 at 0.9, 0.99, 0.999 (slow approach near Delta=1).
* Intermediate inequalities [Numerical]: all ratios to the proved bounds are below one (finite_summary.json, check_continuum.json); the largest is VarZ/q_J=0.975.

## 8. Heisenberg (Task 3)

### 8.1 What fails unconditionally [Proved, tracked]

As gamma->0, B=pi/(2gamma), D'~u^2, L_cut~u, K_rep~2171/gamma^2. The divergence enters only through Lemma 4.1 (the factor B from the fixed-point inequality ||D||<=1/(2L)+q||D||) and through Lemma 3.6 (the 1/gamma, which is natural, since r0~1/(2gamma)). All continuum constants (Lemmas 3.1-3.5) are uniform: q_J<1/2, (3.1)<=3, ab<=pi^2, a2<=2pi^2, Vv<=(2+40/pi)pi^2, and the ratios D/b, aD1/D, gamma D are bounded once B is replaced by an absolute constant. So the only non-uniform estimate is the cumulative root bound.

### 8.2 Theorem C [Proved, conditional on (M)]

(M): for every 0<gamma<pi/2, L, n, the ground-state roots have sigma_L>=0 on [-Q_L,Q_L] (at even half filling, on R).

Under (M): (N) holds (Lemma 2.2), so (I) holds at both parities. z_L is increasing on [-Q_L,Q_L], so |Cum R|<=1/(2L) (C136 (4.x)) and |nu-m|<=q_J/(2L)<=1/(4L). Let G=Cum(mu-r_(Q_L)dx). Integrating (I+K_Q)(mu-r)=R gives

    G+K_Q G = S + c_0 - (nu-m)psi,  psi in [0,q],  S=Cum R,

so G=g_1+c_0 Z with g_1=(I+K_Q)^(-1)[S-(nu-m)psi], and G(-Q)=0 forces |c_0|<=2||g_1||. By (3.1), ||G||<=3||g_1||<=27/(4L). Then sup|sigma_L-r|<=[K(0)/4+2K(0)27/4]/L, i.e. D:=55/(4pi gamma), and D1:=27.25M_1. Repeating Section 5 with these (the ends now use |Cum R|<=1/(2L) and int_(x_j2)^(Q_L) r<=[4c_r D/b+27]/L):

    V_f<=400.9,  ends<=391.5,  curvature<=0.62,  K_large^M<=443,

and the core is nonempty for L>=77. With seam gluing for L<77,

    L[F_(2L)(2n)-2F_L(n)]/M_0 <= max(4*77, 5*443) = 2215

for all 0<=Delta<1, all L>=2, 0<=n<=L/2. For fixed L both sides are continuous in Delta (finite matrices), so the bound holds at Delta=1. With the C129 reductions,

    L(E_(2L)-2E_L) <= (2215+pi^2/4) M    on J_z/J in [0,1], all fields, DM, L.

This is a conditional decision of the Heisenberg endpoint: no logarithmic obstruction appears at the level of an O(1/L) bound; the marginal operator only slows the approach of the coefficient to pi^2/4.

### 8.3 Numerics on the endpoint [Numerical]

Half-filled defect L(F_(2L)-2F_L)/M_0 (heis_scan.json), decreasing in L at every Delta:

| L | 0.5 | 0.9 | 0.99 | 0.999 | 1 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 3.255 | 3.864 | 3.987 | 3.999 | 4.000 |
| 4 | 2.282 | 2.696 | 2.782 | 2.790 | 2.791 |
| 8 | 2.100 | 2.471 | 2.550 | 2.557 | 2.558 |
| 64 | 2.0415 | 2.3896 | 2.4658 | 2.4732 | 2.4741 |
| 1024 | 2.040529 | 2.385142 | 2.460654 | 2.468099 | 2.468922 |
| 4096 | 2.040525 | 2.384805 | 2.460071 | 2.467517 | 2.468340 |
| limit pi v/2 | 2.040524 | 2.384593 | 2.459171 | 2.466579 | 2.467401 |

The supremum is 4 at L=2, Delta=1. Through Delta=1 the minimum of L sigma_L over the roots stays bounded below in natural units (1.48 at Delta=1, rapidities of XXX), and sigma_L>0 on the whole line.

### 8.4 The single remaining estimate

(M) Monotone counting function. For 0<Delta<1 (and at Delta=1 in XXX rapidities), every L>=2 and 1<=n<=L/2, the ground-state Bethe roots satisfy

    sigma_L(x) = a1(x) - (1/L) sum_k K(x-x_k) >= 0   for |x|<=Q_L.

[Numerical] True in every computed case, with sigma_L>0 on all of R, Delta in {0.1,0.5,0.9,0.99,0.999,1}, L up to 4096 (doubled 8192), fillings from 1/L to 1/2. The tightest point is the outermost root, where L sigma_L(x_n) is about 1.44 (Delta=1/2), 3.29 (0.9), 33 (0.999, x-units) and 1.48 (Delta=1, lambda-units). The far tail coefficient sin(gamma)-sin(2gamma)(1/L)sum e^(2x_k) of sigma_L is positive and tends to zero at half filling like L^(-q), so global positivity is a statement about cancellation at the edge, not about the sup-norm error, which is of the same order as sigma_L there.

(M) gives, at once, the identification at both parities, the uniform cumulative bound, and Theorem C. Given (I), a weaker version already makes the constant uniform: int_[-Q_L,Q_L](sigma_L)_- dx <= c/L with c absolute, since it bounds the overshoot of the sawtooth Cum R by c/L and the argument of 8.2 then runs with 1/2 replaced by 1/2+c. Neither version is proved.

[Assumed, unresolved] Without (M), Theorem A's constant diverges like 1/(1-Delta), and Heisenberg is not reached.
