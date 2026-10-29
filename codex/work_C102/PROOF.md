# C102 derivations

2026-10-29

Claims marked [Proved] are finite-dimensional statements or deductions with their hypotheses stated. [Proved conditional] uses a specified Bethe/CFT identification. [Numerical] denotes a finite computation, not a thermodynamic theorem. [Assumed] identifies continuum or phase input.

## 1. What energy concavity actually gives

[Proved] Write H_L(h)=sum_j h_{j,j+1} on a periodic ring. For every normalized vector psi, its energy is linear in h. Taking the infimum proves

    E_L(t h1+(1-t)h2) >= t E_L(h1)+(1-t)E_L(h2),  0<=t<=1.

The thermodynamic limit e(h)=lim E_L(h)/L inherits concavity. Both functions are positively homogeneous, hence superadditive. Positive homogeneity is meant for nonnegative scale factors. Multiplication by a negative number exchanges ground and highest energies and is not covered by it.

Define the finite-ring diagnostic

    A_L(h)=(6/pi)L[L e(h)-E_L(h)].

[Assumed] When used as a CFT coefficient, A=lim A_L is taken along a compatible periodic vacuum sequence. Incommensurate fillings or twists can introduce order-1/L shell terms. All comparisons below either use even zero-field XXZ rings, spin-1 vacuum lengths, or explicitly commensurate hopping phases. A gapped ring with exponentially small corrections has A=0. Assigning it this value extends the Casimir function beyond nontrivial CFTs; the main XXZ counterexample needs no such extension.

Let D_e=e(h_t)-t e(h1)-(1-t)e(h2), and define D_E analogously. Both are nonnegative. Direct subtraction gives

    A_L(h_t)-t A_L(h1)-(1-t)A_L(h2)
       =(6/pi)L[L D_e-D_E].                         (1)

The two nonnegative defects have no prescribed ordering. In the additive version, with G_e=e(h1+h2)-e(h1)-e(h2) and G_E defined analogously,

    A_L(h1+h2)-A_L(h1)-A_L(h2)
       =(6/pi)L[L G_e-G_E].                         (2)

This is the missing sign in the proposed concavity argument. It cannot be supplied by concavity of either energy separately.

[Proved] On any convex cone on which A exists and is positively homogeneous, convexity and subadditivity are equivalent. Apply subadditivity to t h1 and (1-t)h2 for one direction, and apply midpoint convexity to h1,h2 and then multiply by two for the other. The same argument exchanges concavity with superadditivity.

## 2. A counterexample entirely inside the critical XXZ cone

Use spin operators S=Pauli/2 and put

    h(J,K)=J(Sx tensor Sx+Sy tensor Sy)+K Sz tensor Sz,
    J>0, -J<K<J,  delta=K/J=cos gamma, 0<gamma<pi.

[Assumed, standard Bethe/CFT input] The complete zero-field continuum has c=1 and

    A(J,K)=v(J,K)=(pi J/2) sin(gamma)/gamma.          (3)

The normalization is Eq. (13) of [Sirker and Bortz](https://arxiv.org/pdf/cond-mat/0511272), with their J restored. Their Eq. (1) is the spin convention used here. This input is also recorded in C101. The endpoint delta=-1 has quadratic modes and is excluded. The limit delta=1 is the Heisenberg CFT, with marginal corrections that do not change A.

[Proved conditional] Take h_-=h(1,-1/2), h_0=h(1,0), h_+=h(1,1/2). Since h_0=(h_-+h_+)/2, Eq. (3) gives

    A(h_-)=3 sqrt(3)/8,
    A(h_0)=1,
    A(h_+)=3 sqrt(3)/4,
    A(h_0)-[A(h_-)+A(h_+)]/2
       =1-9 sqrt(3)/16=0.0257214207425...>0.        (4)

Positivity follows without a decimal approximation from 256>243. Equivalently,

    A(h_-+h_+)=2 > 9 sqrt(3)/8=A(h_-)+A(h_+).

All three points, their sums, and the entire connecting segment have one gapless boson. This is not a phase-boundary or varying-central-charge effect.

[Proved] Their optimized scale is M=1/2. For a reflection-symmetric bond h, exchange of the two sites fixes h and sends D(o)=o tensor I-I tensor o to -D(o). Convexity and unitary invariance of the half spectral width w imply

    w(h)<= [w(h+D(o))+w(h-D(o))]/2=w(h+D(o)).

Thus M=w(h). The XXZ eigenvalues are K/4, K/4, -K/4+J/2, -K/4-J/2. On |K|<=J the extreme pair is the last two, and M=J/2. The violation therefore occurs on a constant-M face as well as at fixed J.

## 3. The full XXZ curvature has the opposite sign

[Proved conditional] Define f(delta)=sin(gamma)/gamma, delta=cos(gamma). Differentiate using d gamma/d delta=-1/sin(gamma):

    f'(delta)=(sin gamma-gamma cos gamma)/(gamma^2 sin gamma),
    f''(delta)=[2-gamma cot gamma-gamma^2 csc^2 gamma]
                 /(gamma^3 sin gamma)
              =-q(gamma)/(gamma^3 sin^3 gamma),
    q(gamma)=gamma^2+gamma sin gamma cos gamma-2 sin^2 gamma.

Let r(gamma)=sin gamma-gamma cos gamma. Since r(0)=0 and r'=gamma sin gamma>0 on (0,pi), r>0 and f'>0. To determine the second derivative, compute

    q'=2 gamma+gamma cos(2 gamma)-3 sin gamma cos gamma,
    q''=4 sin gamma(sin gamma-gamma cos gamma)=4 sin gamma r>0.

Together with q(0)=q'(0)=0, two integrations give q>0, hence f''<0 everywhere in the open interval. This proves strict concavity, rather than checking its sign on a grid.

Now A=(pi/2)J f(K/J) is the perspective of this concave function. Its Hessian is

    Hess A=(pi/(2J)) f''(delta)
              [[delta^2,-delta],[-delta,1]],        (5)

which is negative semidefinite, with the radial direction as its only null direction. More explicitly, for h_i=h(J_i,K_i), set w_i=J_i/(J_1+J_2). Then

    A(h1+h2)=(pi/2)(J1+J2) f(w1 delta1+w2 delta2)
              >= A(h1)+A(h2).                      (6)

Equality holds if and only if delta1=delta2. The same statement applies along segments. Consequently no convex subcone of this XXZ cone containing two different anisotropies can make A subadditive. This includes the antiferromagnetic subcone 0<=K<=J. For example, delta=0,1/4,1/2 give a midpoint defect 0.00433760379358...>0.

## 4. What monotonicity survives, and how it proves the XXZ bound

[Proved conditional] At fixed J, A strictly increases with K throughout the critical cone. It also increases with J at fixed K, while remaining in the cone, since

    partial_K A=(pi/2) f'>0,
    partial_J A=(pi/2)(f-delta f')
       =(pi/2)[gamma-sin gamma cos gamma]/(gamma^2 sin gamma)>0.

For the last sign, the numerator vanishes at zero and has derivative 2 sin^2 gamma>0. These are statements about specified physical couplings, not operator ordering.

Because M=J/2,

    A/(pi M)=f(delta),
    0<f(delta)<=1,                                  (7)

and the ratio increases strictly towards Heisenberg. Thus comparison at fixed M reduces this particular family to its Heisenberg endpoint. That endpoint has h=S.S=P/2-I/4, so it is a two-level permutation bond. At M=1 its class is [P], with A=pi. This establishes A<=pi M log_2 d for this family, with equality only at delta=1 and its scale/unitary equivalents. The mechanism is monotonicity, not maximization of a convex A.

[Proved] The class [P] is also an exposed extreme point of the full quotient ball. With Pi_s and Pi_a the symmetric and antisymmetric projectors for d=2, put B=Pi_s/3-Pi_a. It is traceless, has equal partial traces, has trace norm two, and sign(B)=P. Thus ell([h])=Tr(Bh)/2 descends to the quotient, has value at most one on every contraction, and reaches one only at the representative P. The quotient-ball construction in Section 6 then proves exposure.

There is a stronger supporting-plane inequality. Expansion at delta=1 gives

    f(delta)=1-(1-delta)/3-(1-delta)^2/45+O((1-delta)^3).

Concavity places the graph below its limiting tangent there:

    f(delta)<=(2+delta)/3,
    A(J,K)<=pi(2J+K)/6<=pi J/2=pi M.               (8)

For delta<1 the first inequality is strict. This uses a single supporting linear functional, not a uniform bound on the gradient. It is consistent with the round-one obstruction. Near delta=-1, A is asymptotic to J sqrt((1+delta)/2), whereas partial_K A diverges. Euler's pairing J partial_J A+K partial_K A=A remains finite.

[Proved] There is no nonconstant gauge-invariant function monotone under the full positive-semidefinite order. Given h,k, one can always choose a scalar c with k+cI>=h. Scalar invariance and monotonicity would force A(k)>=A(h) for every pair, and reversing the pair forces equality. Positive homogeneity would then force the constant to be zero.

A concrete critical counterexample to increasing operator monotonicity is

    h(1,-1/2)+I/8-h(1,0)=-Sz tensor Sz/2+I/8>=0,

but A decreases from 1 to 3 sqrt(3)/8. Its eigenvalues as an increment are 0 and 1/4. The increment with the opposite SzSz sign is also positive semidefinite and increases A. Neither order direction survives in general.

## 5. Restricted convexity that does hold

[Proved] Suppose C is a convex cone of bonds with a common compatible length sequence, A exists there, and e is linear on C. Then A is convex and subadditive on C. Indeed D_e=0 in (1), so its right side is nonpositive at every L. Taking limits preserves that inequality. Equivalently, G_e=0 in (2).

A sufficient condition for e to be linear is a common infinite-volume translation-invariant state omega that minimizes every h in C. Then e(h)=omega(h). Frustration-free cones with compatible periodic minimizing states satisfy the stronger identity E_L=L e, hence A=0. A common bulk minimizer alone does not assert that stronger finite-ring identity.

This restricted criterion supplies a possible extreme-point argument only for C intersect {M<=1}, and only when the decomposing bonds also have the required Casimir limits. Its extreme points need not be extreme points of the full quotient ball. No nontrivial interacting critical cone meeting this criterion was identified here.

[Proved, explicit free example] Another restricted class is complex hopping with no field or interaction:

    h(t)=t |01><10|+conjugate(t)|10><01|.

Jordan-Wigner gives a dispersion 2|t| cos(k+arg t). The two Fermi velocities have magnitude 2|t|. Thus c=1 and A(t)=2|t| for t!=0. At t=0 the finite-ring coefficient is zero, although the flat Hamiltonian is not a CFT.

The local spectrum is -|t|,0,0,|t|. Gauge optimization cannot reduce its half-width. Compression of h(t)+D(o) to {|01>,|10>} is

    [[a,t],[conjugate(t),-a]],  a=o_00-o_11,

with eigenvalues +/-sqrt(a^2+|t|^2). A scalar shift does not change their separation, and the min-max principle bounds the full width below by this compressed width. Therefore M=|t|. It follows that A=2M is convex and subadditive on this two-real-dimensional subspace. Its unit ball is the disk |t|<=1. Its extreme circle obeys A=2<pi, proving the d=2 proposed bound throughout this restricted class. This is consistent with, and weaker in scope than, the general quasi-free result already recorded in REPORT section 9.71.

The velocity statement holds for every hopping phase. Reading it as a periodic vacuum Casimir requires phase-compatible lengths or an appropriate boundary twist. No generic untwisted finite-size limit is asserted for an irrational phase.

[Proved] This class also shows that A is not globally concave or superadditive, even when all compared bonds are critical. Choose

    t_+=exp(i pi/3)/2, t_-=exp(-i pi/3)/2, t_mid=1/4.

Then A(t_+)=A(t_-)=1 and A(t_mid)=1/2. The whole straight segment has Re(t)=1/4 and never reaches t=0. For the three displayed bonds, L divisible by six removes the endpoint hopping twists and gives a common even vacuum sequence. Thus the reverse Jensen sign is an actual periodic-ring comparison at d=2, not just an appeal to dispersions with incompatible fluxes.

## 6. Gauge envelopes and the extreme-bond proposal

[Proved] On a ring, sum_j D(o)_{j,j+1}=0 as an operator. Adding alpha I changes E_L by alpha L and e by alpha. Hence A_L, A, and M are invariant under the stated gauge. In particular,

    sup_{alpha,o} A(h+alpha I+D(o))=A(h).            (9)

An upper envelope over these gauges cannot repair convexity. Maximizing A/w over representatives just minimizes w, giving A/M when A>=0. Taking a supremum over on-site unitary conjugations also leaves A unchanged.

[Proved] There is a valid geometric statement behind the proposed reduction. Let V be Hermitian two-site bonds modulo scalars and D(o). Since half-width w(X)=min_alpha ||X-alpha I||,

    M([h])=inf_{alpha,o} ||h+alpha I+D(o)||.

This is a finite-dimensional quotient norm. Its unit ball is the image of the Hermitian operator-norm unit ball. A minimizing representative exists because the gauge subspace is closed and the problem is finite dimensional.

The extreme points of the Hermitian contraction ball are the self-adjoint unitaries F=2P-I. To see this, an eigenvalue strictly between -1 and 1 can be changed slightly in either direction, giving a nontrivial average. Conversely, if F=(B+C)/2 with -I<=B,C<=I and F has only eigenvalues +/-1, each extremal eigenspace forces B and C to agree with F there, so B=C=F. Every Hermitian contraction is a finite convex combination of these sign matrices by diagonalizing it and expressing its diagonal vector as a convex combination of cube vertices. Therefore each extreme point of the quotient ball has a two-level representative. Images of extreme points need not all remain extreme.

If A were convex on this entire ball, existed on all the decomposing bonds, and obeyed the desired bound on those two-level representatives, the bound would extend to the whole ball. Each of these hypotheses matters. In particular, critical bonds alone do not form a convex ball.

[Proved] Failure is visible on exposed extreme points themselves. For F_x=X tensor X and F_y=Y tensor Y with Pauli X,Y, M(F_x)=M(F_y)=1. The functional

    ell_a([h])=Tr(F_a h)/4

is well-defined on V because F_a is traceless and has vanishing partial traces. For every contraction B, ell_a(B)<=1, with equality only for B=F_a. Indeed equality saturates every eigenvector expectation of B in the +/-1 eigenspaces of F_a. Thus [F_a] is an exposed extreme point of the quotient ball.

Each ring made from F_x alone or F_y alone is a commuting antiferromagnetic Ising chain. On even rings its ground energy is -L, so e=-1 and A=0 identically. But

    h_mid=(F_x+F_y)/2=2(Sx tensor Sx+Sy tensor Sy),
    M(h_mid)=1,  A(h_mid)=2.                        (10)

For this XX bond E_L=-2 csc(pi/L) on even rings. Expanding gives E_L=-2L/pi-pi/(3L)+O(L^-3), proving A=2 directly. Both endpoints obey the proposed bound, as does the midpoint, but their actual A values cannot be propagated by a convexity argument. In fact no convex majorant agreeing with A on all exposed extreme bonds can dominate A at this midpoint. Any useful convex majorant must already overestimate A on some such extreme bonds.

This does not prove that a global maximum of A/M cannot occur at a permutation bond. Heisenberg still saturates the conjectured d=2 maximum. It rules out deducing that conclusion by convex interpolation of the physical endpoint values. General projector bonds remain unclassified.

## 7. Spin-1 and XYZ interpretation of the checks

[Assumed] The critical zero-field XYZ sheets are XXZ sheets and their appropriate unitary/staggered images. This standard identification is described after Eq. (2) of [Fendley and Hagendorf](https://arxiv.org/pdf/1006.0237); their overall sign and Pauli conventions differ from ours. Consequently Eq. (3) transfers to each chosen critical sheet, with the compatible ring sequence. Off those sheets a generic XYZ bond is massive. We did not fit a small nonzero coefficient and call it a new CFT.

[Proved] For spin 1, put X1=S_i.S_{i+1}. Its total-spin 0,1,2 eigenvalues are -2,-1,1. Therefore h_b=X1+b X1^2 has eigenvalues

    -2+4b, -1+b, 1+b,

and reflection makes their half-width gauge optimal. In particular M=3 at b=-1 and M=1 at b=1. At b=1, h=P+I, fixing the ULS normalization.

[Assumed, inherited Bethe/CFT data from C35 and C101] For the literal h_b convention, the TB point b=-1 has e=-4, c=3/2, v=2pi, and A=3pi. This follows from h_-1=4 sum_{j=0}^2 H_j Pi_j-6I and the inherited velocity pi/2 in harmonic-projector units. The ULS point b=1 has e=2-log 3-pi/(3sqrt(3)), c=2, v=2pi/3, and A=4pi/3. The open interval -1<b<1 is the Haldane phase. At b=1/3 the AKLT energy E_L=-2L/3 is exact on the rings used, hence A=0 without needing an extrapolation.

[Assumed] The b>1 portion sampled here is the extended massless phase. [Itoi and Kato](https://arxiv.org/abs/cond-mat/9605105) derive the SU(3)_1 description near ULS and its marginal corrections. We use the established phase identification and do not claim to prove it, or a velocity formula at b=1.5,2, from our short rings.

[Numerical] At b=1,1.5,2 the L=6,9 two-size coefficients have positive midpoint defect 0.08460897. A quadratic fit in 1/L^2 using L=3,6,9 gives 0.07322004. These are evidence for the same Jensen direction as XXZ on this particular segment, not a proof of concavity in the thermodynamic spin-1 critical phase. The TB and ULS coefficients provide integrable calibrations. The Haldane samples show why short-size Casimir fits are not reliable phase diagnostics. Detailed values and residuals are in NUMERICS.md.

## 8. Scope of the conclusion

[Proved conditional] Global subadditivity and convexity fail inside a single ordinary critical XXZ phase, even on its antiferromagnetic part and at constant optimized M. [Proved] The gauge envelope is identical to A, and interpolation from extreme-bond A values fails. A has the opposite Jensen sign on an explicit free hopping segment, so global concavity is not a replacement.

[Proved conditional] The useful surviving comparison is monotonicity and concavity on the XXZ cone, which reaches the Heisenberg bound by a supporting plane. [Proved] Convexity does hold when the bulk energy is linear, and for the zero-field complex-hopping subspace, where its restricted extreme-point proof gives A<=2M. None of these statements proves the proposed inequality for arbitrary interacting bonds or arbitrary projectors.
