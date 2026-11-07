# C140 proof: monotone counting function on the repulsive side

2026-11-07. Run by a Claude Opus agent, not Astra. J=1. Notation of C139 (PROOF Section 0): Delta=cos gamma, 0<gamma<=pi/2, b=pi/gamma, a1, K, r0(x)=1/(2gamma cosh bx), q_J=(pi-2gamma)/(2(pi-gamma)), the hole kernel J=(I+K*)^(-1)K* (C139 Lemma 3.1: J=k1*r0>=0, even, nonincreasing in |x|, ||J||_1=q_J<1/2), and the dual continuum density s_Q=r0+J*(1_out s_Q), out={|x|>Q} (C139 Lemma 3.2). Roots x_1<...<x_n, symmetric, mu=(1/L)sum delta_(x_j), sigma=sigma_L=a1-K*mu. Gaps g_j=[x_j,x_(j+1)], j=1,...,n-1. Tags [Proved], [Numerical], [Assumed].

## 1. Three exact facts

(F1) [Proved] Dual identity on R: sigma = r0 + J*(sigma dx) - (1/L) sum_j J(.-x_j).
Proof: sigma+K*sigma = a1+K*(sigma dx-mu); apply (I+K*)^(-1)=I-J* (C139 Lemma 3.1) and use (I+K*)^(-1)a1=r0 and (I+K*)^(-1)K*=J*. [Numerical: check_lemma1.json, relative residual below 1.5e-3, the J-grid error.]

(F2) [Proved] Gap masses. The Bethe equations give z_L(x_j)=I_j/L, so int_(g_j) sigma = 1/L for every gap.

(F3) [Proved] Tail mass. int_(x_n)^infinity sigma = c-nu/2+1/(2L) >= 1/(2L), since z_L(infinity)=c and 2c-nu=(1-gamma/pi)(1-2nu)>=0.

## 2. Kernel constants [Proved]

Jhat=Khat/(1+Khat) with Khat(w)=sinh((pi-2gamma)w/2)/sinh(pi w/2). For 0<A<C, sinh(Aw)/sinh(Cw)<=e^(-(C-A)w) (equivalent to e^((C-A)w)>=e^(-(C-A)w)), so 0<=Khat<=t:=e^(-gamma w), and since u/(1+u) is increasing, 0<=Jhat<=t/(1+t). With dw=-dt/(gamma t), and J^(*m)(0)=sup J^(*m)=(1/pi) int_0^inf Jhat^m dw,

    J(0) <= (1/(pi gamma)) int_0^1 dt/(1+t) = ln2/(pi gamma),
    sum_(m>=2) J^(*m)(0) <= (1/(pi gamma)) int_0^1 t/(1+t) dt = (1-ln2)/(pi gamma),     (2.1)

using sum_(m>=2) t^(m-1)(1+t)^(-m) = t/(1+t). [Numerical: J(0)=0.0700 b at Delta=0.99 against the bound 0.0702 b; the bound is asymptotically sharp as gamma->0.]

## 3. The two-neighbour lemma

Hypothesis (P): sigma>=0 on [x_1,x_n] (n>=2).

Lemma 3.1 [Proved, under (P)]. Write v=sigma on out:={|y|>x_n}, Q=x_n, J_out the operator with kernel 1_out(x)J(x-y)1_out(y), ||J_out||<=q_J, positive.

(a) For x in a gap g_k, sigma(x) >= r0(x) + (J*(1_out v))(x) - (1/L) max(J(x-x_k),J(x-x_(k+1))).
(b) For x>x_n, r0(x)+(J_out v)(x)-J(x-x_n)/L <= sigma(x) <= r0(x)+(J_out v)(x), and symmetrically for x<-x_n.

Proof. By (F1), sigma(x) = r0(x) + int_out J(x-y)v(y)dy + sum_j int_(g_j) J(x-y)sigma(y)dy - (1/L)sum_j J(x-x_j). J(x-.) is nondecreasing on (-inf,x] and nonincreasing on [x,inf). For (a), pair g_j with x_j when j<k (on g_j, J(x-y)>=J(x-x_j), and with sigma>=0 and (F2) the pair is >=0), pair g_j with x_(j+1) when j>k (same), and keep g_k with x_k, x_(k+1): int_(g_k) J(x-y)sigma >= min(J(x-x_k),J(x-x_(k+1)))/L by (F2), and min+max equals the sum. For the lower bound in (b), pair g_j with x_j for all j<=n-1; x_n is left over. For the upper bound pair g_j with x_(j+1) (on g_j, J(x-y)<=J(x-x_(j+1)), so each pair is <=0), and x_1 is left over with sign minus. The contribution of the mirror part of out to x>x_n is included in J_out v. End.

Corollary 3.2 [Proved, under (P)] (outer Neumann bounds). With psi=1_(y>Q)J(y-Q)+1_(y<-Q)J(y+Q) and phit=(I-J_out)^(-1)psi>=0,

    s_Q - phit/L <= v <= s_Q   on out.

Proof. (b) says r0 1_out - psi/L + J_out v <= v <= r0 1_out + J_out v. The resolvent sum_k J_out^k is positive and v is bounded, so v=(I-J_out)^(-1)(r0 1_out - psi/L + e) with e>=0, and similarly above; (I-J_out)^(-1)(r0 1_out)=s_Q on out (C139 Lemma 3.2). End.

Consequences: ||phit||_1 <= ||psi||_1/(1-q_J) = q_J/(1-q_J) <= 1, so

    int_out sigma_- <= 1/L.                                          (3.1)

Also J_out^k psi <= J^(*(k+1))(.-Q)+J^(*(k+1))(.+Q) by positivity, so for every x, using (2.1),

    (1/L) int_out J(x-y) phit(y) dy <= (2/L) sum_(m>=2) J^(*m)(0) <= 2(1-ln2)/(pi gamma L).   (3.2)

Lemma 3.3 [Proved, under (P)] (location of the outermost root). L r0(x_n) >= b(1-q_J)/4 = pi^2/(8 gamma(pi-gamma)).

Proof. By (F3) and Corollary 3.2, 1/(2L) <= int_Q^inf s_Q = h_Q/2, where h_Q=int_out s_Q. Integrating s_Q=r0+J*(1_out s_Q) over out gives h_Q = int_out r0 + q_J h_Q - int_(|x|<Q) J*(1_out s_Q) <= int_out r0 + q_J h_Q. And int_Q^inf r0 = 2atan(e^(-bQ))/(2gamma b) <= e^(-bQ)/(gamma b) <= (2/b) r0(Q). So 1/(2L) <= (2/b) r0(Q)/(1-q_J). End. [Numerical: L r0(x_n)/b >= 0.442 in every case, minimum at (L,n)=(4,2); the proved value is 0.125 to 0.235. location_scan.json.]

Theorem 3.4 [Proved, under (P)]. For x in [x_1,x_n],

    sigma(x) >= r0(x) - (2-ln2)/(pi gamma L) >= theta(gamma) r0(x),   theta(gamma)=1-8(2-ln2)(pi-gamma)/pi^3.

Proof. In Lemma 3.1(a) bound J*(1_out v) >= -(1/L) int_out J(x-y)phit(y)dy (Corollary 3.2, s_Q>=0) and use (3.2) and J<=J(0)<=ln2/(pi gamma): the loss is (ln2+2-2ln2)/(pi gamma L). Then r0(x)>=r0(x_n) and Lemma 3.3. End.

theta(gamma)>0 iff pi^3 > 8(2-ln2)(pi-gamma), i.e. gamma > 0.1760. At gamma_0=0.177, cos gamma_0 = 0.98437 >= 63/64, so theta>0 on 0<=Delta<=63/64. [Numerical: min over gaps of sigma/r0 is 0.974 to 0.994 (location_scan.json) against theta=0.437 (Delta=0.1), 0.294 (0.5), 0.093 (0.9), 0.008 (0.98).]

## 4. Continuity argument: (M) on [x_1,x_n] and identification

Theorem 4.1 [Proved]. For 0<=Delta<=63/64, every L>=2 and 1<=n<=L/2, the branch of real ordered critical points of the Yang-Yang functional continued from the free solution at gamma=pi/2 exists on [gamma,pi/2], is real-analytic, reflection symmetric, and its roots satisfy

    sigma_L >= theta(gamma) r0 > 0 on [x_1,x_n],     int_(|x|>x_n) (sigma_L)_- <= 1/L,

and the Yang-Yang Hessian is positive definite. Hence (N) and C139's identification (I) hold at both parities on this range.

Proof. n=1: x_1=0 by symmetry, and sigma = a1-a2/L with a2=K; a1/a2=(cosh2x-cos2gamma)/(2cos gamma(cosh2x-cos gamma)) lies between 1/(2cos gamma) and (1+cos gamma)/cos gamma, both >1/2>=1/L for 0<gamma, so sigma>0 on R and the 1x1 Hessian 2pi L sigma(0)+theta'(0)>0. For n>=2 let G be the set of gamma_1 in [gamma,pi/2] such that the branch exists on [gamma_1,pi/2] with (P) strict there. pi/2 is in G (K=0, sigma=a1>0). G is open in [gamma,pi/2]: strict (P) on the compact [x_1,x_n] gives sigma(x_j)>0, hence a positive definite Hessian (C139 Lemma 2.2), and the implicit function theorem extends the branch, with (P) strict nearby by continuity. G is closed: at gamma_*=inf G the roots stay in the a priori ball of C139 Section 2.1, a limit point is an ordered critical point with distinct coordinates, and sigma>=0 on [x_1,x_n] by uniform convergence on compacts. Theorem 3.4 then gives sigma>=theta r0>0 there, so (P) is strict, the Hessian is positive definite, and by local uniqueness the branch reaches gamma_*. So G=[gamma,pi/2]. Symmetry follows from local uniqueness and the reflection. (3.1) gives the negative mass. Identification with the sector ground state is C139 Section 2.3 (coordinate Bethe vector near the free point, Perron-Frobenius, open-closed argument), which needs only (N). End.

Remark. Theorem 4.1 is the statement (M) on [x_1,x_n], with an explicit margin, plus an O(1/L) bound on the negative mass beyond. Positivity strictly beyond the outermost root is not proved: the sufficient continuum inequality obtained by dropping the paired gap terms, s_Q(x) >= 2(int_Q^inf s_Q) phit_Q(x) for x>Q, is false far out for Delta>=0.5 (tail_cont.json: the ratio is 0.32 to 0.81 at x-Q=12, while it exceeds 1.6 at x=Q). The true sigma_L is positive there in every computed case (sigma_L/r0>=0.952 on R, ratio_scan_*.json).

## 5. What this gives downstream

5.1 [Proved] On 0<=J_z/J<=63/64, C139 Theorem A, which was conditional only on (I), is now unconditional: with delta=1/64,

    L(E_(2L)-2E_L) <= (5400/delta) M = 345600 M

for every field, DM angle and L>=2.

5.2 The sawtooth. On [x_1,x_n] z_L is increasing, so |Cum R|<=1/(2L) there (C136). On [x_n,Q_L], 0<=nu/2-z_L<=1/(2L)+int_(x_n)^x sigma_- <= 3/(2L) by (3.1), and at even half filling (Q_L=infinity) also nu/2-z_L(x)=int_x^inf sigma >= -1/L. So sup|Cum R| <= 3/(2L) on [-Q_L,Q_L]. This is the weak form of (M) named in C139 Section 8.4 with c=1. [Proved, bookkeeping not redone line by line] Running C139 Section 8.2 with the sawtooth 1/2 replaced by 3/2 multiplies every term by at most 9 (each is of degree at most two in the sawtooth constant: D, D1 and L_cut linear, the ends and the curvature quadratic, V_f invariant), so its constant becomes at most 9*2215 = 19935, and L(E_(2L)-2E_L) <= (19935+pi^2/4) M on [0,63/64]. The uniformity in Delta of this second constant is what matters for the Heisenberg end.

5.3 The remaining sub-estimate for the repulsive class. All of Section 3 is uniform in gamma except the one comparison r0(x_n) L >= b(1-q_J)/4 > (2-ln2)/(pi gamma). Define

    (L*) along the branch, L r0(x_n) >= lambda_* b with lambda_*=(2-ln2)/pi^2 = 0.13241.

[Proved] Under (L*) at some gamma (and (P) along the continuation, which (L*) then propagates), Theorem 3.4 holds with a positive margin, hence Theorem 4.1 and 5.2 hold on all of 0<=Delta<1 with the uniform constant 19935+pi^2/4, and at Delta=1 by finite-matrix continuity. Lemma 3.3 proves L r0(x_n) >= (pi/(8(pi-gamma))) b, which exceeds lambda_* b exactly when gamma>0.1760; as gamma->0 it tends to b/8=0.125 b, short of 0.1324 b by 6 percent. [Numerical] L r0(x_n)/b >= 0.442 for every L from 4 to 1024, every sampled filling, and Delta up to 0.9999, a margin of a factor 3.3. The loss in Lemma 3.3 is the factor 1/(1-q_J), which ignores the part of J*(1_out s_Q) that returns to [-Q,Q]; since J has a peak of width gamma and a broad part of width about 1/2, a lower bound on that returned mass uniform in Q was not obtained.

## 6. The ferromagnetic side (D)

[Not proved] For -1<Delta<0, K=-T<=0 and the hole kernel J=-(I-T*)^(-1)T* is negative with norm B-1=pi/(2eta)-1, so neither the pairing argument of Section 3 (which uses J>=0 and unimodality) nor the positive Neumann comparisons transfer. Positivity of sigma is automatic there (sigma=a1+T*mu>0), and (D) is an accuracy statement, not a sign statement, so a single lemma covering (M) and (D) was not found. The C138 numerics stand: eta L sup|sigma_L-r_Q| in [0.23,0.32] for L=16 to 8192, 1+Delta from 1e-2 to 1e-6.

## 7. Consistency checks

* Delta=0 [Proved]: K=J=0, sigma=a1=r0, theta(pi/2)=1-4(2-ln2)/pi^2=0.47, the bounds hold trivially.
* L=2, n=1 [Proved]: sigma=a1-a2/2>0, with tail coefficient sin gamma(1-cos gamma)>0 vanishing only at Delta=1, matching the XXX value 0 (C139).
* Delta=1/2 and the Casimir value pi v/6, Delta=-1 zero defect: inherited unchanged from C139 and C138 (no new energy estimate here).
* The proved lower bound sigma >= r0-(2-ln2)/(pi gamma L) is checked directly in check_lemma1.json (bulk margin positive in all 28 cases) together with Lemma 3.1(a) itself.
