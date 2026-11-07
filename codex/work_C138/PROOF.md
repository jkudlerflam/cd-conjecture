# C138 proof notes: the ferromagnetic crossover of the XXZ ring

2026-11-07. Written by a Claude Opus agent, not Astra. J=1 in the analytic sections; restore J by homogeneity. Tags: [Proved] is a mathematical deduction, [Numerical] is floating point without enclosures, [Heuristic] is an argument with an unproved step, [Assumed] is a hypothesis.

Summary of what is and is not established here.

- [Proved] Task 1. The nondegenerate limit near gamma=pi is reached by dividing the equations by eta, not by rescaling rapidities. The rapidities stay of order one, only the kernel narrows to width eta, and the crossover parameter is eta L, the ratio of kernel width to root spacing. The eta->0 limit at fixed L is a strictly convex log-gas whose energy reproduces degenerate first-order perturbation theory exactly (Section 2).
- [Proved] An unconditional uniform bound |F_L(n)-L e_bulk(n/L)| <= (5pi/4) eta for all L, n and -1<Delta<=0 (Section 3). This gives (X) on the crossover side eta L <= A_0 with K_* = (5pi/4)A_0, but not beyond.
- [Proved, conditional] A reduction of (X) at all L to five explicit estimates, one discrete (D) and four continuum ones (D1), (W1)-(W3), with an explicit K_* (Section 4). No length threshold is needed, so the C129 band is not used.
- [Numerical] All five hypotheses hold in every sampled case with margin, uniformly down to epsilon=10^-6 (Section 5). The quantity C136 controlled, the cumulative root error, is not uniformly O(1/L) at a finite Fermi boundary: it grows slowly as eta->0 (Section 5.3). The pointwise density error is the right currency.
- [Assumed, unresolved] (D) is not proved. It is the exact missing estimate.

## 1. Conventions and eta-normalized constants

[Proved, inherited] All notation is that of C136 Section 1: gamma=arccos Delta in [pi/2,pi), the Bethe equations with I_j=j-(n+1)/2, F_L(n)=Delta L/4+G_L(n), G_L(n)=sum_j e(x_j), the kernels a1=p'/(2pi), T=-theta'/(2pi)>=0, the restricted objects r_Q, Z_Q, epsilon_Q, m(Q), E(Q), h(Q), the counting function z_L, counting density sigma_L, the counting boundary Q_L with z_L(Q_L)=nu/2, and R_L=mu-1_[-Q_L,Q_L] sigma_L dx. We write

    eta = pi-gamma in (0,pi/2],    Delta = -cos eta,
    epsilon = 1+Delta = 2 sin^2(eta/2) <= eta^2/2.

In C136's notation alpha=eta. The following are exact.

    q = integral T = 1-2eta/pi,        B = 1/(1-q) = pi/(2eta),
    a = pi sin(eta),                   b = pi/gamma in [1,2),
    integral a1 = eta/pi,              a1(0) = tan(eta/2)/pi,
    e(0) = h(0) = -epsilon,            T(0) = cot(eta)/pi,
    ||T'||_1 = 2T(0) = 2cot(eta)/pi.                              (1.1)

The values of a1(0), T(0) follow from a1(x)=tan(eta/2) sech^2 x/[pi(1+tan^2(eta/2) tanh^2 x)] and T(x)=tan(eta) sech^2 x/[pi(tanh^2 x+tan^2 eta)]. T is even and decreasing on x>0, so ||T'||_1=2T(0).

[Proved, inherited from C136] The following hold for every Delta in (-1,0] with constants independent of Delta. r_Q and Z_Q are even and nonincreasing on [0,Q], a1<=r_Q<=r_0=1/(2gamma cosh bx), 1<=Z_Q<=B (C136 (2.2)). |nu-m(Q_L)|<=(B-1)/L (C136 (4.2)). epsilon_Q is even, epsilon_Q(+-Q)=0, epsilon_Q'>0 on (0,Q], and |epsilon_Q'|<=2a r_Q (C136 (5.1)-(5.2)). e_bulk''>=0 and e_bulk'(nu)=h, with h(0)=-epsilon and h->0 at nu=1/2 (C136 (5.3)). The exact identity

    F_L(n)-L[E(Q_L)+h(Q_L)(nu-m(Q_L))] = L integral epsilon_Q dR_L     (1.2)

(C136 (4.3)), and at even half filling the same with Q=infinity, h=0, m=nu=1/2.

## 2. Task 1: the rescaling and the nondegenerate limit

### 2.1 The kernel symbol

[Proved] With f^(w)=integral exp(iwx)f(x)dx, the transform of C136 (6.2) at beta=iw gives

    T^(w) = sinh((pi-2eta)w/2)/sinh(pi w/2),
    a1^(w) = sinh(eta w/2)/sinh(pi w/2),
    1-T^(w) = 2cosh(gamma w/2) sinh(eta w/2)/sinh(pi w/2).          (2.1)

The last line is the identity sinh A-sinh C=2cosh((A+C)/2)sinh((A-C)/2) with A=pi w/2 and C=(pi-2eta)w/2. Define the normalized symbol

    N_eta(w) = [1-T^(w)]/eta.

Since 0<=T^(w)<=T^(0)=q for real w, N_eta(w)>=2/pi for every w and every eta, and the bound is attained at w=0. Thus the normalized operator eta^(-1)(I-T) is uniformly coercive on L^2 with constant 2/pi, and on any interval eta^(-1)(I-T_Q) has L^infinity inverse norm at most eta B=pi/2 (C136 (2.1)). As eta->0 at fixed w,

    N_eta(w) -> w coth(pi w/2),     a1^(w)/eta -> (w/2)/sinh(pi w/2),

and the limit symbol w coth(pi w/2) is a first-order elliptic symbol, comparable to 2/pi+|w|. Its inverse applied to the limit source gives 1/(2cosh(pi w/2)), the transform of 1/(2pi cosh x)=r_0 at gamma=pi. For |w|>>1/eta, instead, N_eta(w) ~ 1/eta, so the normalized operator is not uniformly elliptic beyond frequency 1/eta.

[Proved, interpretation] The divergence B=pi/(2eta) of the Neumann inverse is matched exactly by the source: a1 and e-h are of order eta and eta^2. Consequently r_Q<=1/(2gamma), eta Z_Q<=pi/2, and |epsilon_Q|<=B sup|e-h|<=(pi/(2eta))epsilon<=pi eta/4 are all uniform in their natural normalization. No rescaling of rapidities is needed: the root density lives on |x| of order one (order log L in the half-filled tails) at every eta, while the kernel has width eta. The losses in C136 come from applying the norm B to inputs that are not of order eta, namely the discreteness errors of the roots, which are of order 1/L.

### 2.2 Step-free Bethe equations

[Proved] For y not zero, theta(y)=-2 atan(tanh(y)/tan(eta)) equals -pi sign(y)+phi(y) with

    phi(y) = 2 atan(tan(eta)/tanh(y)),

because atan(u)+atan(1/u)=pi/2 for u>0. For ordered roots, sum_(k!=j) sign(x_j-x_k)=2j-n-1=2I_j. Hence the Bethe equations are equivalent to

    L p(x_j) = sum_(k!=j) phi(x_j-x_k),        j=1,...,n.           (2.2)

The quantum numbers have disappeared. Both sides are of order eta when the root separations exceed eta, and phi is bounded by pi near the diagonal. This form avoids the O(1) cancellation between L p and the step part of theta, and it is what bethe_ferro.py solves.

### 2.3 The limit eta->0 at fixed L

[Proved] Divide (2.2) by eta. As eta->0, p(x)/eta -> tanh x and phi(y)/eta -> 2coth y locally uniformly on y!=0, and both quotients are real-analytic and even in eta near eta=0 at fixed x, y!=0. The limit system

    L tanh(x_j) = 2 sum_(k!=j) coth(x_j-x_k)                          (2.3)

is the critical-point equation of

    W(x) = L sum_j log cosh x_j - 2 sum_(j<k) log sinh(x_k-x_j)

on the ordered cone. W is strictly convex there (log cosh and -log sinh are strictly convex) and tends to +infinity at the boundary of the cone and at infinity: collisions send -log sinh to +infinity, and on the ray where the top k roots move to +infinity together the asymptotic slope is kL-2k(n-k)=k(L-2n+2k)>0 for n<=L/2, and the same for the bottom k. The slope is linear on each sign sector of ordered directions, which is generated by these top-k and bottom-k rays, so it is positive on every ray (the sign-adapted decomposition of C136 Section 1). So (2.3) has a unique solution x^(0), and the Hessian of W is positive definite there. The implicit function theorem gives a real-analytic branch x(eta)=x^(0)+O(eta^2) of solutions of (2.2) for small eta, ordered, with the correct quantum numbers. By uniqueness of the real Bethe solution (C134, C136 Section 1) it is the ground-state root set.

[Proved] Sum rule. Multiply (2.3) by tanh x_j and sum over j. Using (tanh u-tanh v)coth(u-v)=cosh(u-v)/(cosh u cosh v)=1-tanh u tanh v,

    L sum_j tanh^2 x_j = sum_(j!=k) [1-tanh x_j tanh x_k]
                       = n(n-1) - (sum_j tanh x_j)^2 + sum_j tanh^2 x_j.

Summing (2.3) itself gives sum_j tanh x_j=0 by antisymmetry. Hence sum tanh^2 x_j^(0)=n(n-1)/(L-1) and

    sum_j sech^2 x_j^(0) = n(L-n)/(L-1).                            (2.4)

Since e(x)=-(eta^2/2)sech^2 x+O(eta^4) and epsilon=eta^2/2+O(eta^4),

    G_L(n) = -epsilon n(L-n)/(L-1) + O_L(eta^4).                     (2.5)

[Proved, consistency check] This is exactly the degenerate first-order energy of the ferromagnetic multiplet state used in C129 Section 10, J epsilon[L/4-n(L-n)/(L-1)] above J Delta L/4 after subtracting -JL/4. [Numerical] (2.4) holds to 1e-8 at eta=1e-5 for (L,n)=(8,3),(16,8),(40,13),(101,50), and the roots solve (2.3) to 3e-8.

[Proved, interpretation] At fixed eta the rapidities interpolate between the log-gas (2.3) when root spacing, about 1/L, exceeds eta, and the continuum density r_Q when root spacing is below eta. The crossover is at eta L of order one, which is epsilon L^2 of order one. In the sparse regime the counting density sigma_L is a sum of Lorentzian spikes of height about 1/(pi eta L) at the roots; in the dense regime it is smooth.

### 2.4 Which C136 constants are artifacts

[Proved, bookkeeping] The C136 large-size estimate rests on three bounds whose eta-dependence is as follows.

1. Cumulative error ||G||_infinity<=3B/(2L), which is of order 1/(eta L).
2. Pointwise density error ||sigma_L-r_Q||<=||T'||_1||G||, of order 1/(eta^2 L) (C136 D/L with D of order eta^(-2)).
3. Curvature e_bulk''<=a, of order eta, used against |nu-m|<=B/L, giving a mismatch cost aB^2/(2L) of order 1/(eta L).

[Numerical] The true sizes are, uniformly in the sampled range: L||G|| between 0.67 and 0.76 at half filling (Section 5.2); eta L sup|sigma_L-r_Q| at most 0.32 (Section 5.1); e_bulk''/eta^2 at most 4/pi (Section 5.4). So items 2 and 3 lose a full factor 1/eta. Item 1 loses a factor 1/eta at half filling, and a smaller but apparently unbounded factor at a finite Fermi boundary.

[Proved] Even with item 1 sharpened to O(1/L) at half filling, the C136 density route (item 2) cannot give a uniform constant at any L: the clipping point must satisfy L r_0(x_c) of order 1/eta^2 for item 2 to give a relative error, and the first-order cost of the tail beyond it is 2a integral_(x_c)^infinity r_0, of order 1/(eta L). This is why C136's coefficient diverges even for L>>1/epsilon. A bound on sigma_L-r_Q of order 1/(eta L), not 1/(eta^2 L), is what is needed.

## 3. An unconditional uniform bound

[Proved] Theorem A. For every -1<Delta<=0, every L>=2 and every 0<=n<=L/2,

    |F_L(n) - J L e_bulk(Delta,n/L)| <= (5pi/4) eta J.               (3.1)

Proof. n=0 is exact. Under t=z_L(x) the quadrature cells [t_j-1/(2L), t_j+1/(2L)], t_j=I_j/L, tile [-nu/2,nu/2], so the right side of (1.2) is sum_j [f(t_j)-L integral_cell_j f dt] with f(t)=epsilon_Q(z_L^(-1)(t)). Each bracket is at most the oscillation of f on its cell, so

    |L integral epsilon_Q dR_L| <= Var_[-Q,Q](epsilon_Q) = 2|epsilon_Q(0)|
          <= 2 integral_0^Q 2a r_Q = 2a m <= a = pi sin(eta).          (3.2)

Here we used that epsilon_Q is even, increasing on [0,Q] and zero at Q, and |epsilon_Q'|<=2a r_Q. At half filling Var(epsilon)=2a r_0(0)=a/gamma<=a. For the mismatch, convexity and h in [-epsilon,0] give

    0 <= L[e_bulk(nu)-e_bulk(m)-h(m)(nu-m)] <= L|nu-m| |h(nu)-h(m)|
       <= (B-1) epsilon <= (pi/(2eta))(eta^2/2) = pi eta/4.          (3.3)

Adding (3.2) and (3.3) and using sin eta<=eta gives (3.1). End of proof.

[Proved, corollary] For eta L<=A_0 this is (X) with K_*=(5pi/4)A_0, two-sided, at all densities. In terms of epsilon, eta L<=A_0 contains epsilon L^2<=A_0^2/2. The matched doubling gives L[F_2L(2n)-2F_L(n)]<=(15pi/4)eta L J, so in units of M_0=J/2 the sector constant (15pi/2)A_0 on eta L<=A_0. This is an independent Bethe-ansatz version of the C129 band; it grows linearly in eta L instead of quadratically and is better than C129's epsilon L^2 when eta L>15pi.

[Proved, interpretation] (3.1) shows the size of the room in the problem. The trivial first-order quadrature bound is already O(eta), uniformly in L, and the target O(1/L) asks for one additional factor 1/(eta L), not for the full second-order accuracy (whose natural size is eta/L).

## 4. Reduction of (X) to a pointwise density estimate

### 4.1 A three-region quadrature lemma

[Proved] Lemma 1. Let the cells be split into a set C of core cells and the rest. Then

    |L integral epsilon_Q dR_L|
       <= (1/(8L)) Var_K(epsilon_Q'/sigma_L) + integral_U |epsilon_Q'| dx,     (4.1)

where K is the union of the x-images of the core cells and U the union of the others.

Proof. For a core cell the Peano kernel of the midpoint rule, (h/2-|s|)^2/2<=h^2/8 with h=1/L, gives |f(t_j)-L integral_cell f|<=Var_cell(f')/(8L), and f'=epsilon_Q'/sigma_L in the x variable; variation is invariant under reparametrization. For the other cells use the oscillation bound of Section 3. Sum over cells. End of proof.

This replaces C136's clipping. There are no jumps and no tail-mass estimate; the non-core cells pay only the variation of epsilon_Q over them.

### 4.2 The hypotheses

Fix Q=Q_L (or Q=infinity at even half filling) and write r=r_Q, rho=M/(eta L) with M=4C_D, and K_c={|x|<=Q: r(x)>=rho}, an interval [-x_c,x_c] since r is even and nonincreasing on [0,Q].

[Assumed] (D) There is a constant C_D with |sigma_L(x)-r_Q(x)|<=C_D/(eta L) for x in K_c.

[Assumed] (D1) There is a constant C_D1 with |sigma_L'(x)-r_Q'(x)|<=C_D1/(eta^2 L) for x in K_c.

[Assumed] (W1) Var_[-Q,Q](v_Q)<=C_W eta for every Q in (0,infinity], where v_Q=epsilon_Q'/r_Q.

[Assumed] (W2) For Q>=1 and |x|<=Q-1, r_Q(x)>=kappa r_0(x).

[Assumed] (W3) sup over 0<=nu<=1/2 of e_bulk''(nu)<=C_e eta^2.

All constants are meant independent of eta in (0,pi/2], L, n. At half filling (Q=infinity) (W1) holds with C_W eta replaced by 2ab<=4pi sin eta, since v=ab tanh(bx) [Proved]; (W2) holds with kappa=1 and (W3) is not needed. (D) and (D1) are statements about the Bethe roots; (W1)-(W3) are statements about the restricted continuum equation alone.

### 4.3 The conditional theorem

[Proved, conditional on (D), (D1), (W1)-(W3)] Theorem B. For every -1<Delta<=0, L>=2 and 0<=n<=L/2,

    |F_L(n) - J L e_bulk(Delta,n/L)| <= K_* J/L,

    K_* = (8a/(9M))[C_D + C_D1(1+1/(kappa b))/eta] + C_W eta/6
          + 4aM(1+2/(kappa b))/eta + 16a/3 + C_e eta^2 (B-1)^2/2.      (4.2)

With a<=pi min(1,eta), b>=1, (B-1)eta<=pi/2 and M=4C_D,

    K_* <= (2pi/(9C_D))[C_D + C_D1(1+1/kappa)] + pi C_W/12
           + 16 pi C_D (1+2/kappa) + 16pi/3 + pi^2 C_e/8.             (4.3)

Proof. On K_c, (D) and r>=rho=4C_D/(eta L) give |sigma_L-r|<=r/4, so sigma_L>=3r/4 and g=r/sigma_L lies in [4/5,4/3]. Take as core cells those whose x-image lies in K_c. At most one cell straddles each endpoint of K_c.

Core. epsilon_Q'/sigma_L=v_Q g. By the product rule for variation and |v_Q|<=2a,

    Var_K(v g) <= 2a Var_(K_c)(g) + (4/3)Var(v_Q).

Next, g'=[r'(sigma_L-r)-r(sigma_L'-r')]/sigma_L^2 and sigma_L^2>=(9/16)r^2, so by (D), (D1)

    Var_(K_c)(g) <= (16/9)[(C_D/(eta L)) integral_(K_c)|r'|/r^2
                          + (C_D1/(eta^2 L)) integral_(K_c) dx/r].

Since r is even and nonincreasing on [0,Q], integral_(K_c)|r'|/r^2=2[1/r(x_c)-1/r(0)]<=2/rho. For integral dx/r split K_c at |x|=Q-1. On |x|<=min(x_c,Q-1), (W2) gives 1/r<=1/(kappa r_0) and integral_0^y dx/r_0=(2gamma/b)sinh(by)<=1/(b r_0(y)); since r_0(x_c)>=r(x_c)>=rho this part is at most 1/(kappa b rho). The remaining part of [0,x_c] has length at most 1 and r>=rho there. So integral_(K_c) dx/r<=2(1+1/(kappa b))/rho, and

    Var_(K_c)(g) <= (32/(9M))[C_D + C_D1(1+1/(kappa b))/eta].

With (W1) the core contribution to (4.1) is at most [2a Var(g)+(4/3)C_W eta]/(8L), the first two terms of (4.2).

Tail. The non-core cells lie in the tail T_c={x_c<|x|<=Q} together with the K_c-parts of at most two straddling cells. On a straddling cell's K_c-part, r<=(4/3)sigma_L and the cell has sigma_L-mass 1/L, so the integral of |epsilon_Q'|<=2a r there is at most 8a/(3L); two such cells give 16a/(3L). On the tail r<rho. Per side, the part with Q-|x|<1 contributes at most rho integral of length one; the part with Q-|x|>=1 lies where r>=kappa r_0 by (W2), so r_0<rho/kappa there and its integral is at most integral_(r_0<rho/kappa) r_0 over one side, which is at most 2rho/(kappa b) by integral_y^infinity r_0<=exp(-by)/(gamma b)<=2r_0(y)/b. Hence integral_(T_c) r<=2rho(1+2/(kappa b)) and integral_(T_c)|epsilon_Q'|<=4a rho(1+2/(kappa b)), the third term of (4.2) after the factor L. If Q<1 the second part is empty and the bound still holds. If K_c is empty, the core is empty and every cell is a tail cell; the same tail bound applies, since then r<rho everywhere.

Mismatch. As in (3.3) but with (W3): L[e_bulk(nu)-e_bulk(m)-h(m)(nu-m)]<=L C_e eta^2(nu-m)^2/2<=C_e eta^2(B-1)^2/(2L), using C136 (4.2). This is the last term.

Combine (1.2), Lemma 1 and the three pieces. At half filling the mismatch is zero and the same estimates hold with r=r_0, kappa=1, Var(v)=2ab. For (4.3), bound sin(eta)<=1 in the C_D term, sin(eta)/eta<=1 in the others, and (B-1)eta<=pi/2. End of proof.

[Proved, conditional corollary] With matched doubling (C136 (6.1)) and C129's field, DM and twist reductions, every bond of the U(1) class with -1<=J_z/J<=0, every field and DM angle and every L>=2 satisfies

    L(E_(2L)-2E_L) <= [5K_*+pi^2/4] M(B_field).                     (4.4)

Delta=-1 has zero defect (C129 Section 8). No length threshold A is needed, so C=5K_*+pi^2/4. If the hypotheses are only granted on (-1,-1+delta_0], use C136 on [-1+delta_0,0] and take the larger of 5K_*+pi^2/4 and 6*10^15 delta_0^(-25/3).

[Numerical, constant] Inserting the conservative values C_D=1/2, C_D1=2, kappa=7/10, C_W=8, C_e=13/10 (each above the largest sampled value, Section 5) into (4.3) gives K_*<=125 and C<=628. The dominant term 16 pi C_D(1+2/kappa), about 97, is the first-order tail cost; the core term is under 8.

## 5. Numerical status of the hypotheses and of (X)

All runs are floating point; scripts and JSON output are in this folder. Anisotropies epsilon=1+Delta in {10^-2,...,10^-6}, with spot checks at epsilon in {0.1,0.5,0.9}. Bethe energies agree with ED to 5.8e-15 at L<=12 (bethe_ferro.py). Typhon was not used; nothing on the critical path required it.

### 5.1 (D) and (D1)

[Numerical] diag_half.py (half filling, L=16 to 8192): eta L sup over roots |sigma_L-r_0| is 0.230 to 0.317, increasing toward 1/pi as eta decreases; it is attained at isolated roots in the far tail, where sigma_L has a spike of height about 1/(pi eta L). Between roots it is at most 0.088. The C136 bound on the same quantity, eta D, is 10.8, 33.6, 106, 335, 1061 at epsilon=10^-2,...,10^-6. eta^2 L sup|sigma_L'-r_0'| is at most 0.011.

[Numerical] dens_check.py (finite Fermi boundary, seven cases with epsilon=10^-3 to 10^-6, eta L from 2.9 to 29): eta L sup over [-Q_L,Q_L] of |sigma_L-r_Q| is at most 0.28, attained 1 to 12 kernel widths from the Fermi point. dens_check2.py extends this to epsilon in {0.9, 0.5, 0.1, 10^-2, 10^-4, 10^-6}, L in {16, 64, 256, 1024} and n from 2 to L/2-2 (120 cases): the maximum of eta L sup|sigma_L - r_Q| over [-Q_L,Q_L] is 0.317, 0.302, 0.185, 0.054, 0.009, 0.002 at epsilon = 10^-6, 10^-4, 10^-2, 0.1, 0.5, 0.9. For (D1), eta^2 L sup|sigma_L' - r_Q'| on [0, Q_L - min(eta,Q_L)/4] is at most 0.22 for epsilon <= 0.1. A first pass at epsilon in {0.5, 0.9} gave values up to 3.6, growing linearly in L; this was a mesh artifact (a one-sided difference for r_Q' at x=0, where r_Q'(0)=0 by evenness: the error halves when the mesh is doubled and does not depend on L). With r_Q'(0)=0 imposed the maximum is 0.019 (dens_check2_0.9_0.5.json). The small-epsilon values are upper estimates that may contain the same artifact.

### 5.2 The cumulative root error

[Numerical] At half filling L sup|G| is 0.67 to 0.76 for every epsilon and L, against C136's B/2, which is 5.5 to 555. At a finite Fermi boundary (diag_identity.py, mesh-converged cases) the cumulative error is concentrated at the outermost roots and reaches L sup|G| of about 2 at epsilon=10^-5 and 3 to 4.4 at epsilon=10^-6, while the density mismatch dN=L(nu-m(Q_L)) is 1.6, 3.3 and 4.6 at (epsilon,L,n)=(10^-4,2048,768), (10^-5,2048,768), (10^-6,2048,512) (conv_check.py, mesh-converged to 1 to 3 percent). At fixed eta L between 18 and 29 it grows as 0.32, 0.77, 1.6 at epsilon=10^-2, 10^-3, 10^-4, roughly like eta^(-0.7), still far below the proved bound B-1. This growth is smooth on the kernel scale, since (D) stays bounded. So a uniform cumulative bound L||G||<=C is apparently false near a Fermi boundary as eta->0, and C136's cumulative currency cannot be repaired by sharper constants alone. [Heuristic] The mechanism is the Euler-Maclaurin endpoint correction of the root sum, a source of mass about 1/(L^2 sigma_L(Q)) localized within eta of the Fermi point, amplified by the restricted resolvent, whose low-frequency norm is B.

### 5.3 (W1)-(W3)

[Numerical] cont_diag.py and cont_diag_wide.py (Q from 0.02 to 8, epsilon from 1 to 10^-6): Var(v_Q)/eta is at most 8.0, attained at the free point where it is exact; e_bulk''/eta^2 is at most 1.273=4/pi, again attained at the free point (e_bulk''=pi sin(pi nu) there), and is within 1.05 of eta^2 near Delta=-1; min over |x|<=Q-1 of r_Q/r_0 is at least 0.762 (Q=1, epsilon=10^-6) and is 0.76 to 0.82 across epsilon, converging as eta->0. The edge ratio r_Q(Q)/r_0(Q) scales as eta^(1/2) (0.43, 0.256, 0.147, 0.083, 0.047 at epsilon=10^-2,...,10^-6), so a bound of the C136 form kappa r_0<=r_Q up to the edge cannot be uniform; (W2) away from the edge suffices. Two cases at epsilon=10^-6, Q>=4 show r_Q/r_0 up to 1.05, which violates the proved r_Q<=r_0 and marks them as underresolved by the mesh.

### 5.4 (X) itself and the doubling defect

[Numerical] scan_defect.py: for every sampled (epsilon, L, n), with L=4 to 4096 paired to 2L, 0 < L[F_2L(2n)-2F_L(n)]/M_0 <= pi v/2, the limiting half-filled Casimir value, with v=(pi/2)sin(gamma)/gamma. The ratio to pi v/2 increases monotonically in L at fixed density, crosses one half near eta L of order 3 to 10, and reaches 0.9993 at epsilon=10^-4, L=4096. Below L=8 the ratio is at most 0.83 at epsilon=0.1 and 0.004 at epsilon=10^-6. No C129 band or dilute bound is violated. The ladder identity L(F_L-L e_bulk)=-(1/2) sum_j Lambda_j/4^j along doublings, with Lambda_j=2^j L[F_(2^(j+1)L)-2F_(2^j L)], gives |L(F_L-L e_bulk)|<=0.0386, 0.0118, 0.0036, 0.0010, 0.00018 at epsilon=10^-2,...,10^-6, at most pi v/6. Independently, the identity (1.2) evaluated with the continuum solver agrees with the ladder to 4e-7 at epsilon=10^-3, L=64. So (X) holds numerically with K_* of order eta, three orders of magnitude below the conditional (4.3).

## 6. Where a proof of (D) would have to come from

[Proved] Exact representation. On (-Q_L,Q_L), sigma_L-r_Q=(I-T_Q)^(-1)(T*R_L), where T*R_L(x)=integral T(x-y)dR_L(y) is the midpoint-quadrature error of the function y->T(x-y) against the roots. This follows from sigma_L=a1+T*mu, mu=sigma_L dx+R_L on the interval, and r_Q=a1+T_Q r_Q.

[Heuristic] In the dense region the high-frequency part of T*R_L is exponentially small in eta L sigma_L, since z_L is analytic in a strip of width of order eta and the cell structure is periodic in t; the restricted resolvent acts on it with norm of order one, not B. The low-frequency part is an Euler-Maclaurin edge term of size 1/(L^2 sigma_L) localized at the Fermi points, together with a bulk term of relative size (eta L sigma_L)^(-2). Applying the crude norm B to the bulk term gives (D) only for eta L>=C log(1/eta); avoiding the log requires the Fourier-multiplier structure of the full-line resolvent at frequencies above 1/eta and its Wiener-Hopf form near the edges. In the sparse region the spike structure gives (D) directly with C_D near 1/pi. A uniform version of these two statements, with explicit constants, is the missing input. (D1) would follow from the same argument with one more derivative.

[Heuristic, sketch at half filling; the tail teeth are not written out] A discrete maximum principle reduces the half-filled cumulative bound to a local sawtooth defect. With g_j=z_L(x_j)-z_0(x_j) and Gamma the interpolation of the g_j linear in the continuum counting coordinate, G=Gamma+L^(-1)S with |S|<=1/2, and sup_j g_j<=B sup_j |(T*S)(x_j)|/L. The local defect T*S at a root is of order (local spacing) in the dense region and of order eta/(spacing) times a logarithm in the sparse region. [Heuristic] Neither is uniformly of order eta, so this route alone gives the half-filled cumulative bound only up to a factor log(1/(eta L)) in the sparse regime, although numerically the cumulative error stays below 0.76/L.

## 7. Task 3 and the window

[Proved] In the window epsilon^(-1/2)<<L<<epsilon^(-1), that is 1<<eta L<<1/eta, the first-order perturbative energy of Section 2.3 is the eta->0 limit of the Bethe energy at fixed L, and Theorem A bounds the full deviation by (5pi/4)eta. [Heuristic] Second-order perturbation theory around the multiplet has energy denominators of order 1/L^2 and is useful only for epsilon L^2 small, so it does not reach the window; the Bethe reduction of Section 4 covers the window directly, conditional on (D). This is why Task 3 was not pursued further.
