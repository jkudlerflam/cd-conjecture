# C126. Permutation twists and the doubling defect

2026-11-06. All new work is confined to this directory. Claims are tagged below. There was no independent referee review.

## 1. Setting and vacuum compatibility

[Definition] For L>=2, write H_L=sum_x h_(x,x+1), E_L=min spec H_L, and M=inf_(alpha,o)||h+alpha I+o tensor I-I tensor o||. A replica site is a rung containing N copies of C^d. Its permutation representation is the ordinary tensor-factor representation of S_N, unless a fermionic representation is explicitly specified.

[Assumed, infrared input] A vacuum-compatible sequence for a critical model means that the unrestricted physical ground state is in the conformal identity sector, with the correct fermionic spin structure when needed, and

    E_L=e0 L-pi A/(6L)+o(1/L),    A=sum_a c_a v_a.

Whenever a fixed N is used, both L and NL must belong to this sequence. This is a sector prescription, not a uniform bound on the onset of the expansion. For the computations it means even L for the three spin-1/2 antiferromagnets and L divisible by 3 for the fundamental SU(3) chain. For a classically ordered model it means that a minimizing periodic word fits on the ring. For the ferromagnetic Ising example every L fits. These explicit choices do not classify vacuum compatibility for every Hermitian bond.

[Assumed, ordered quantum phases] Where an exponential finite-size correction is invoked, the separate assumption is |E_L-e0 L|<=C L^p exp(-L/xi) along a sequence fitting the ordered period. No uniform C or xi across Hamiltonians is assumed. Vacuum compatibility by itself is not used as a theorem giving this estimate.

## 2. Exact seam identity

[Proved] Label sites (x,a), 0<=x<L and a in Z_N. Move the seam to (L-1,0). Put

    B=sum_a h_((L-1,a),(0,a)),
    B_sigma=sum_a h_((L-1,a),(0,a+1)),
    H_sigma=H_L^(N)-B+B_sigma.

If sigma_0 moves copy a to a+1 on rung 0, then B_sigma=sigma_0 B sigma_0^dagger. The ordered list

    (0,0),...,(L-1,0),(0,1),...,(L-1,1),...,(L-1,N-1)

identifies H_sigma with H_(NL), retaining the orientation of every bond even when h is not reflection symmetric. Therefore

    min spec H_sigma-N E_L=E_(NL)-N E_L.                 (1)

This modifies one bond sum. Conjugating the whole H_L^(N) by sigma_0 would also modify the other bond incident on rung 0 and would preserve the original spectrum. That operation would not give (1).

[Proved] Scalar shifts cancel from (1). A bond coboundary also cancels on the twisted ring. Its seam contribution is sum_a(o_(L-1,a)-o_(0,a+1)), and the sum on rung 0 is unchanged by the cycle. Thus M is the correct gauge-invariant scale.

[Proved] The thermodynamic limit of the ground-energy density exists for a bounded finite-range chain, by cutting and joining open intervals with a bounded boundary error. Taking N to infinity at fixed L in (1) gives

    lim_N (E_(NL)-N E_L)/N=e0 L-E_L.                    (2)

This identity does not require a CFT or an exchange of limits. At N=2, (1) is the target doubling defect.

[Numerical] Small oriented complex bonds were checked with the actual on-site seam conjugation and with the long-ring matrix. See algebra_checks.py and algebra_checks.json.

## 3. CFT check

[Proved within CFT] The covering map w=z^(1/N) has Schwarzian

    {w,z}=(1-1/N^2)/(2 z^2).

In the covering vacuum, each sheet contributes (c/12){w,z} to its stress tensor. Summing the N sheets gives the twist-field weights

    h_sigma=hbar_sigma=(c/24)(N-1/N).

This is the branch-point construction of [Calabrese and Cardy](https://arxiv.org/abs/hep-th/0405152). The replica theory on circumference L has vacuum subtraction -N pi c v/(6L). Exciting its twist ground state adds (2pi v/L)(h_sigma+hbar_sigma), so

    E_sigma-E_untwisted=(pi c v/(6L))(N-1/N).

Adding independent sectors yields

    E_(NL)-N E_L=(pi A/(6L))(N-1/N)+o(1/L)             (3)

for each fixed N along the compatible sequence. At N=2 this is pi A/(4L). Dividing by N and comparing the separate exact limit (2) gives the Casimir coefficient pi A/6 in variant (a). The lattice remainder in (3) is not asserted to be uniform in N.

## 4. Free copies and a spread twist

### 4.1 The symmetry and its boundary condition

[Proved] For a number-conserving quadratic chain, identical fermionic copies have a U(N) symmetry on their copy indices. In a paired quadratic chain a general complex copy rotation does not preserve the pairing tensor. Real O(N) copy rotations always preserve identical Majorana quadratic forms. An odd N-cycle has determinant +1 and belongs to SO(N), so odd N suffices for the spread construction, including the large-N subsequence. An unsigned even cycle need not belong to the identity component of O(N).

[Proved] An ordinary spin swap and a fermionic mode swap are different operators. For two one-mode factors, the mode swap acts by -1 on |11>, while the tensor swap acts by +1. Jordan-Wigner strings and parity-dependent closing bonds therefore matter. A mode rotation is not automatically a local on-site symmetry of the tensor product of spin chains. Even the replicated XX spin Hamiltonian does not acquire the claimed physical on-site U(N) by ignoring these strings.

[Definition, scope of the free construction] Work either in a CAR ring with a fixed periodic or antiperiodic boundary condition, using odd N, or in a spin chain whose standard Jordan-Wigner representation is quadratic in each global parity sector. Choose a translation eigenvector that is a Gaussian ground vector of definite parity. These are the nearest-neighbor quasi-free realizations considered here. An additional non-parity projection, if present in another model, requires a separate admissibility check.

[Proved, parity matching in the spin realization] Write p=+1 or -1 for the parity of the selected L-ring ground vector. In the standard Jordan-Wigner convention its closing quadratic bond has sign b=-p. Use a real signed cyclic matrix R whose successive link signs have product s=(-1)^(N-1). Its determinant is (-1)^(N-1)s=1, so it lies in SO(N) even when N is even. Unwrapping the replica seam bR gives a single fermionic ring with total boundary sign

    b^N s=(-p)^N (-1)^(N-1)=-p^N.

The trial has total parity p^N, since the copy rotations are even fermionic unitaries. This is the physical closing sign for that parity sector of the NL spin ring. The intermediate link signs can be moved to the final closing bond by real sign changes of the fermion modes, which preserve both hopping and pairing forms. The trial is therefore a valid state of the physical long ring. The long-ring ground state need not have the same parity as the trial, since the variational principle minimizes over both physical sectors. This matching establishes admissibility, rather than assuming that the trial already occupies the long-ring vacuum sector.

[Proved] For odd N, R can be the unsigned cycle. For even N, it must contain an extra minus sign. Its eigenphases are the centered phases in section 4.2. The signed CAR construction produces a trial on the same physical NL spin ring, but is not the ordinary on-site tensor permutation on the original replica rungs.

### 4.2 Trial energy

[Proved] Let epsilon_hop=<h_hop>_L be the hopping energy of one bond in the Gaussian vacuum. Its N-copy covariance is diagonal and identical in the copy indices. For odd N choose the real antisymmetric logarithm K of the cycle with eigenvalues i phi_r, where

    phi_r=2pi r/N,    r=-(N-1)/2,...,(N-1)/2.

The rotations O_x=exp(x K/L) distribute the seam holonomy into an equal relative rotation exp(K/L) on every link. On-site quadratic terms are invariant. The cross-copy two-point function vanishes before the rotation, so the expectation of each hopping bond is epsilon_hop times the trace of its relative copy rotation. Consequently the trial cost is

    T_N(L)=L epsilon_hop [sum_r cos(phi_r/L)-N].        (4)

For number-conserving fermions the same calculation uses U(N). The sine terms cancel because the eigenphases occur in opposite pairs, even if the original bond carries a current. For paired chains the real trace supplies this cancellation directly. No derivative of an interacting ground energy is being identified with a free expression.

[Proved] Use 1-cos t<=t^2/2 and

    sum_r phi_r^2=(pi^2/3)(N-1/N).

The inherited bond-scale theorem of REPORT section 9.71 gives ||h_hop||<=M, both for complex hopping, ||h_hop||=||B||_1, and for Majorana hopping, ||h_hop||=||B||_1/2. Therefore

    E_(NL)-N E_L <= T_N(L)
                     <= (pi^2 M/(6L))(N-1/N).         (5)

Using |epsilon_hop| in the upper estimate makes it valid without a sign assumption. The parity argument above establishes the first inequality for the stated spin realization. For a CAR ring with fixed boundary condition, odd N gives that inequality without altering its boundary condition. Comparing (5) with (3), for any compatible odd N>1, gives

    A <= pi M.                                        (6)

The construction is available for every odd N at fixed L in these classes. Thus (2) also gives E_L>=e0 L-pi^2 M/(6L). In the parity-resolved spin realization the signed cycle extends (4)-(5) to every N, so D(L)<=pi^2 M/(4L) for every length with the stated Gaussian ground vector. This is a finite-size zero-temperature bound, not merely its CFT limit. The infrared conclusion A<=pi M still presumes the physical vacuum expansion when A is identified from ring energies.

[Proved, convention check] The signed cycle has phases

    phi_r=2pi [r-(N-1)/2]/N,    r=0,...,N-1,

and the same squared-phase sum for every N. At N=2 these phases are -pi/2 and pi/2. An unsigned fermion swap instead has phases 0 and pi. In the parity-resolved spin construction its total boundary sign would be +1 although the trial has even parity and requires -1. Using it without this correction would not reproduce the vacuum CFT check.

[Proved, comparison with inherited results] Equation (6) is the diamagnetic, or f-sum, constant. It is weaker than the established sharp bound A<=2M for all nearest-neighbor quasi-free chains, including pairing. That sharper theorem uses the velocity-weighted Fermi-vector trace lemma of section 9.71, which is not a consequence of the norm step in (5). The free property used here is a continuous copy symmetry that preserves the product Gaussian vacuum and makes the bond expectation depend only on the relative rotation. This is the property that fails for the interacting trials below.

[Proved, constants] On the parity-resolved doubled free spin rings, (5) gives C=pi^2/4 and hence A<=pi M. C119's all-temperature estimate has coefficient k=10pi/3, giving C<=10pi/3 at zero temperature and A<=40M/3 if only that estimate is used. Its finite-temperature and parity-resolved scope is stronger than a trial-energy statement. The sharp free infrared slope is C_asymptotic=pi/2, from A<=2M. These three constants concern different strengths of estimate and must not be equated.

[Proved, XX check] At half filling, with h=SxSx+SySy and even L, E_L=-csc(pi/L) and M=1/2. Since

    sum_r cos(phi_r/L)=sin(pi/L)/sin(pi/(NL)),

(4) equals E_(NL)-N E_L. Thus the spread free trial is optimal for XX after the correct sewing. Its thermodynamic hopping energy is -1/pi=-2M/pi, and retaining this value in (4) gives A=2M instead of the coarser norm bound. algebra_checks.py checks odd N and the signed N=2 formula through L=30.

[Numerical] An independent paired-XY check constructs the signed two-copy rotation as an eight-spin vector, then moves its internal fermion signs to the final physical closing bond. With Jx=1, Jy=0.3, field=0.71 and L=4, its direct spin energy cost is 0.05297882357309, agreeing with (4) to 7.51e-15. The exact doubled defect is 0.02968411229538. The base and trial parities are both even. This checks the nontrivial internal sign change in the even-N sewing, beyond the XX energy identity.

[Proved, normalization caution] The inherited stiffness budget sum_alpha D_alpha<=4M assumes Hilbert-Schmidt normalized on-site charges. On a replica rung the dimension is d^N and the raw bond norm scales as NM. Copy charges have their own representation index. Applying the budget with the original M while leaving these charges unnormalized is invalid. The direct contraction (4) avoids that mistake.

## 5. Interacting on-site paths

### 5.1 Exact preservation of the relevant ground sector is obstructed

[Proved] A logarithm of a symmetry need not be an on-site conserved generator. Although G=product_x S_x is the global two-copy swap and G Omega=Omega for Omega=psi tensor psi, this does not imply that Q=sum_x S_x annihilates Omega or commutes with H_L^(2). The logarithm of G can preserve the ground space but is nonlocal in the spatial coordinate. The sum of local logarithms is a different operator.

[Proved, all N on the Heisenberg singlet sequence] No continuous path u(t) on the rung joins the identity to a nontrivial copy permutation while u(t)^(tensor L) preserves the product ground vector up to phase at every t, for even L>=4.

Here is a reduced-state proof. Let e=E_L/L, and use the normalized Pauli basis I,sigma^1,sigma^2,sigma^3 with inner product tr(A^dagger B)/2. The singlet two-site marginal has correlation matrix

    T=diag(1,t,t,t),    t=4e/3.

The nearest-rung marginal of N copies has correlation matrix T_N=T^(tensor N). The eigenvalue t occurs only on

    V=span{sigma^i acting in copy a : 1<=i<=3, 1<=a<=N}.

Indeed all its other eigenvalues are powers t^k, and -1<t<0. Conjugation by any rung unitary acts orthogonally on the Hermitian operator space. Preservation of the two-rung marginal implies O T_N O^T=T_N, so it preserves V. The commutator Lie algebra of V is the direct sum of N simple su(2) ideals. An automorphism of this algebra permutes its simple ideals. That permutation is locally constant along a continuous path. The identity fixes each ideal, while a nontrivial copy permutation moves ideals. Such a path is impossible.

[Proved] The range of t needed above follows without a thermodynamic extrapolation. A product of dimers on an even ring gives e<=-3/8. The three-site operator S1.S2+S2.S3 has minimum -1. Summing its translates gives 2H_L, hence e>=-1/2 for L>=4. Thus -2/3<=t<=-1/2, and all powers used above are separated. The standard singlet ground-sector property of the antiferromagnetic Heisenberg ring is the model input.

[Numerical] At L=6, e=-0.467129272955..., t=-0.622839030607.... The Pauli correlation matrix agrees with diag(1,t,t,t) to 3.34e-16. Its two-copy eigenvalues are 1, six copies of t, and nine copies of t^2. The on-site stabilizer obstruction is not inferred from a failed numerical optimizer.

[Proved limitation] This rules out exact ground-preserving on-site paths. It does not rule out a path that leaves the ground sector by an L-dependent amount, or a unitary with support growing with L. Merely checking P Q P on a one-dimensional ground space gives no control. The leading leakage cost is

    <exp(-itQ) H exp(itQ)>-E
       =t^2 <Q(H-E)Q>+O(t^3).

The complementary matrix elements of Q must therefore be controlled in energy norm. Projecting Q first deletes the quantity that needs a bound.

[Proved] An S_N isotypic decomposition diagonalizes the discrete symmetry representation. It does not imply invariance under arbitrary continuous rotations inside its representation blocks. For N=2, paths central in the two isotypic blocks have the form

    u_x=phase times exp(i theta_x S_x).

Their entire variational family is treated next. The proof above also rules out any other exact ground-preserving rung path, not only this central one.

### 5.2 Exact fractional-swap energy

[Proved for Heisenberg; numerical contraction identity for the other listed models] Let

    |Phi(theta)>=product_x exp(i theta_x S_x) |psi_L>^(tensor 2),
    phi_x=2 theta_x,    z_x=cos phi_x,    w_x=sin phi_x.

Use the fixed twisted Hamiltonian H_sigma. The seam can be encoded by z_L=-z_0, w_L=-w_0. Put p=tr[h(sigma tensor sigma)], e=tr(h rho), and a=p-e, where rho and sigma are the two-site and one-site marginals of psi_L. For the real symmetry-preserving ground states in the tests, the cost is

    F(phi)=sum_(x=0)^(L-1) [a(1-z_x z_(x+1))+b w_x w_(x+1)].    (7)

For Heisenberg, sigma=I/2, p=0, a=-e and b=4e^2/3. To derive this, write S=2 S_A.S_B+1/2. Conjugation of a spin by exp(i theta S) gives cos^2(theta) S_A+sin^2(theta) S_B plus a term proportional to sin(2theta)(S_A cross S_B), with the opposite sign in the other copy. In the product of singlets, <S_x^i S_y^j>=(e/3)delta_ij and one-point functions vanish. The two direct-copy terms give

    2e [cos^2(theta_x)cos^2(theta_y)+sin^2(theta_x)sin^2(theta_y)]
       =e[1+cos(2theta_x)cos(2theta_y)].

The cross-product contractions give 12(e/3)^2 sin(2theta_x)sin(2theta_y). Subtract 2e and obtain (7). In XXZ the same contraction gives b=4 sum_i J_i g_j g_k, where (i,j,k) is cyclic, J=(1,1,Delta), and g_i=<S_x^i S_y^i>.

[Numerical] For SU(3), b is evaluated by an exact d^4-dimensional contraction of rho tensor rho and the two endpoint swaps. checks.py verifies (7) at four unrelated angle pairs. The largest error across all models and sizes is 8.89e-16. The script does not use (7) to obtain its coefficient and then omit the independent contraction check.

[Proved] The uniform path phi_x=pi(x+s)/L has

    F_uniform=L [a-(a-b)cos(pi/L)/2]
             =L(a+b)/2+O(1/L).                         (8)

Its bulk energy density is nonzero. At a uniform halfway rotation, phi_x=pi/2, an ordinary untwisted bond costs a+b. The missing continuous symmetry creates a potential energy along the path, not merely a gradient cost.

[Proved] Even optimizing every angle cannot remove the order-one cost on Heisenberg. More generally, if a>|b|, (7) obeys

    F >= (a/2) sum_x (z_x-z_(x+1))^2
             +(a-|b|) sum_x w_x^2.                     (9)

This uses sum(1-z_x z_(x+1))=sum w_x^2+(1/2)sum(z_x-z_(x+1))^2 and |sum w_x w_(x+1)|<=sum w_x^2, including the antiperiodic seam. If some |z_x|<=1/2, the second term is at least 3(a-|b|)/4. Otherwise every |z_x|>1/2, and the antiperiodic boundary forces an adjacent sign change, giving a squared difference greater than 1. Hence

    F >= min{a/2, 3(a-|b|)/4}.                         (10)

There is a sharper bound. Since |w_x w_(x+1)|<=1-z_x z_(x+1), each summand in (7) is at least (a-|b|)(1-z_x z_(x+1)). The sum of 1-z_x z_(x+1), with z_L=-z_0, is multilinear on [-1,1]^L. Its minimum occurs at vertices, where antiperiodicity forces at least one disagreeing edge. Hence F>=2(a-|b|). For Heisenberg, 3/8<=a<=1/2 and b=4a^2/3, so a-b>=1/6. Every profile, including nonmonotone profiles and every branch of the local logarithm, satisfies

    F >= 1/3,    L even, L>=4.                         (11)

The exact doubling defect is pi^2/(8L)+o(1/L), using the inherited Heisenberg Casimir coefficient A=pi/2. The ratio of the best fractional-swap trial to the exact defect therefore diverges at least linearly with L. This is an analytic obstruction to the whole two-isotypic-phase trial family on a vacuum-compatible interacting critical chain.

## 6. Large N does not supply a continuous symmetry

[Proved] The preceding stabilizer argument already holds for every N in the Heisenberg product vacuum. More copies do not connect the disconnected components of the exact on-site stabilizer. Approximate preservation needs an additional estimate, which that argument neither supplies nor excludes.

[Proved, ordered counterexample to a proposed automatic bulk cancellation] Consider h=-J Z tensor Z, M=J, and L>=3. Choose the valid ring ground vector

    psi_L=(|0...0>+|1...1>)/sqrt(2),    Omega=psi_L^(tensor N).

All lengths are compatible and E_(NL)-N E_L=0. Take N an odd prime and let P be the cyclic permutation of the N spins on one rung. Let u(t)=exp(i t arg P), with the principal eigenphases 2pi k/N, -(N-1)/2<=k<=(N-1)/2. Expand

    u(t)=sum_r a_r(t) P^r,
    a_r(t)=(1/N)sum_k exp[2pi i k(t-r)/N].

The two-rung marginal of Omega is the uniform diagonal sum of |w,w><w,w| over binary words w of length N. A nonconstant word has orbit size N because N is prime. Its transformed one-spin expectation is

    m_a(w,t)=sum_r |a_r(t)|^2 z_(a+r)(w).

Constant words obey the same formula since sum_r |a_r|^2=1. Averaging the products z_(a+r)z_(a+s) over all words gives delta_rs. The energy per copy of an ordinary bond after the same rotation on both rungs is therefore

    e_N(t)/N=-J sum_r |a_r(t)|^4.                      (12)

At t=1/2,

    |a_r(1/2)|^2=1/[N^2 sin^2(pi(r-1/2)/N)],
    sum_r |a_r(1/2)|^4=(N^2+2)/(3N^2).

For completeness, the last sum follows by differentiating
sum_r csc^2(x+pi r/N)=N^2 csc^2(Nx) twice and using (csc^2 x)''=6csc^4 x-4csc^2 x, then setting x=pi/(2N). Thus the excess per copy per bond at halfway is

    (2J/3)(1-1/N^2) -> 2J/3.                          (13)

[Numerical] Direct spin matrices give energies per copy -0.4074074074074074, -0.36, and -0.3469387755102041 at N=3,5,7 and J=1, agreeing with (12) to 1.67e-16. This tests the actual on-site cyclic logarithm, including the short orbits.

[Proved] This example disproves the inference that a large copy number makes fractional cyclic permutations preserve the product vacuum at negligible energy per copy. The interval around this barrier can be chosen uniformly in N. Since ||arg P||<=pi and the replicated bond has norm JN, differentiating the common on-site rotation bounds the derivative of its energy per copy by 4pi J. Equation (13) is at least 16J/27 for N>=3, so an N-independent neighborhood retains a positive excess. A uniform spatial spread traverses a number of such bonds proportional to L. Its neighboring-angle error per copy per bond is O(J/L), by the same commutator norm estimate. Its total excess is therefore of order JNL, not a gradient cost of order JN/L. It does not disprove variant (a). Here variant (a) holds with zero cost, by choosing the same polarized phase in every copy. The failure concerns the proposed construction from the symmetric cat product.

[Proved structural obstruction] Fourier transformation in copy labels does not decouple a replicated interaction. For fermionic density interactions it gives

    V sum_a n_(x,a)n_(y,a)
      =(V/N) sum_(q-p+s-r=0 mod N)
          f^dag_(x,p) f_(x,q) f^dag_(y,r) f_(y,s).

There are three free momentum sums. The factor 1/N is not a small operator-energy bound per copy. An exchangeability or de Finetti statement for the original copies does not replace this interaction by independent Fourier modes. Equation (13) is a direct spin example of the surviving cost.

## 7. Finite-group LSM and ordered phases

[Proved] The on-site permutation representation is linear, U(g)U(h)=U(gh), and has invariant vectors |a>^(tensor N). There is no projective on-site obstruction supplied by S_N alone. The twisted translation can obey T_tilde^L=U(sigma) without producing a nontrivial projective commutator. This kinematics does not estimate the ground-energy difference between the twisted and untwisted Hamiltonians.

[Proved, scope of the literature comparison] Discrete-symmetry LSM results use a projective representation per unit cell to force degeneracy or exclude a unique symmetric gapped phase. They do not provide an arbitrary finite-group twist-energy bound of order 1/L. The twisted-boundary argument of [Yao and Oshikawa](https://arxiv.org/abs/2010.09244) makes this projective hypothesis explicit. Applying that result to the linear replica permutation representation supplies no missing energy estimate.

[Proved] The usual ordered domain wall and the permutation twist must be distinguished. In the Ising chain above an ordinary Z2 spin-flip twist forces one frustrated bond, costing 2J on every ring, including compatible lengths. A permutation twist among identical copies does not flip their order parameters. The diagonal assignment of the same phase to every copy is fixed by all permutations. For the classical Ising example it gives zero energy cost. More generally a minimizing periodic classical word that fits at length L can be repeated N times, so no positive permutation defect is forced on such a compatible length.

[Proved] Nevertheless a chosen vacuum sector can suffer an order-one cost. For two Ising copies fixed to opposite polarizations at reference sites, the permutation-twisted long ring must contain two domain walls, costing 4J. The product of symmetric cat vectors has seam expectation 2J above 2E_L. Both original rings have vacuum-compatible lengths. These are a constrained-sector energy and a trial-state cost, respectively. Neither is the unrestricted D(L), which is zero. Selecting or correlating the phase labels is necessary before using an ordered-phase argument.

[Proved conditional on the ordered-phase assumption in section 1] If |E_L-e0 L|<=C L^p exp(-L/xi), then for each fixed N,

    |E_(NL)-N E_L|
      <=C (NL)^p exp(-NL/xi)+N C L^p exp(-L/xi).

The actual unrestricted permutation defect then cannot approach a nonzero domain-wall tension. This is not a uniform dimension-only bound, since C and xi are uncontrolled. For arbitrary interacting bonds no theorem eliminating all finite-size order-one defects from a sector label alone was proved here.

## 8. Numerical trial comparison

[Numerical] checks.py diagonalizes balanced conserved-charge sectors, using deterministic starting vectors, tolerance 2e-12, and an explicit eigenpair residual. It computes two-site marginals of the L-ring vectors and evaluates (7) through an independent small-matrix contraction. The angles are optimized from 17 initial profiles with analytic gradients. This is a best-found value, not a certification of the global minimum. The lower bound (11) supplies the optimizer-independent Heisenberg conclusion.

[Assumed standard model input] The unrestricted antiferromagnetic ground state is in the zero-magnetization sector for the even XXZ rings used here, and in the equal-color sector for the SU(3) rings of length divisible by 3. Full-space ED checks these choices at L=4 for all spin parameters and at L=3,6 for SU(3). Larger sizes use this standard ground-sector property rather than diagonalizing every sector.

[Proved, inherited trial] For comparison, the C114 coherent seam family is included without presenting it as a new route. If q=2(p-e)=2a and rho_r is the interval marginal, its best energy for the real states here is

    T_coherent=q/S_L,    S_L=sum_(r=0)^(L-1) tr(rho_r^2),    rho_0=1.

It coherently sums all prefix-swap seam positions. C114 proves that the energy-defect matrix in this span is q times the identity and that its Gram matrix has largest eigenvalue S_L. It also already proves failure of this entire span on XX. The present computation extends the interacting finite-size comparison, not the conceptual route.

[Numerical] The largest spin calculation is 2L=20, sector dimension 184756. The largest completed SU(3) calculation is 2L=12, sector dimension 34650. The maximum eigenpair residual is 1.06e-11, and full-space comparisons differ by at most 1.34e-15. Exact contractions of the fractional functional differ by at most 8.89e-16. All numerical trial costs exceed the computed exact defect. The detailed values and ratios are in TABLES.md and checks_local.json.

[Numerical] From L=4 to L=10 the optimized fractional/exact ratio grows from 2.2876 to 5.7003 for Heisenberg, from 2.4198 to 6.0957 for Delta=0.5, and from 2.7190 to 8.7458 for Delta=2. The coherent/exact ratios grow from 1.1095 to 1.3728, 1.1460 to 1.4167, and 1.2506 to 1.9055, respectively. SU(3) gives optimized fractional ratios 1.4343 and 3.4203 at L=3,6, and coherent ratios 1.1474 and 1.4959. Two SU(3) sizes are not an asymptotic fit.

[Heuristic, not used in the proof] The usual critical interval-purity scaling would give a coherent/exact ratio growing as L^(c/4), namely L^(1/4) for the c=1 examples and L^(1/2) for SU(3), with finite-size and marginal corrections. The gapped Delta=2 defect eventually decays faster, making a ratio to it an especially stringent test. A divergent ratio alone would not exclude an O(1/L) upper bound when the exact defect is exponentially small. Our analytic obstruction uses critical Heisenberg, where this distinction does not arise.

[Operational limitation] The SU(3) 18-site balanced sector has 18!/(6!)^3=17153136 states. It was not run locally. run_typhon.slurm requests 128 GB and executes the matrix-free extension of checks.py inside SLURM. No Typhon connection or submission was attempted. The L=9 SU(3) trial-to-exact comparison remains pending that computation.

## 9. One missing lemma and its status

[Assumed sufficient lemma, unresolved] A quantitative integrated seam-screening estimate would complete the general case. Here is one precise version which permits states beyond on-site paths. Let

    H(s)=H_L^(2)+s V,    V=B_sigma-B,    0<=s<=1,
    epsilon(s)=min spec H(s),    q_*=epsilon'_+(0).

Let mu=-epsilon'' be the positive curvature measure of this concave function on (0,1). The required lemma is

    integral_(0,1) (1-s) dmu(s) >= q_*-C(d) M/L         (14)

for every vacuum-compatible L, with a finite dimension-only C(d). In a smooth nondegenerate interpolation its left side is

    2 integral_0^1 (1-s) sum_(n>0)
       |<n_s|V|0_s>|^2/[E_n(s)-E_0(s)] ds.             (15)

Degenerate points are handled by the curvature measure, including its derivative jumps. If the s=0 ground state is unique, q_* is the product-ground seam expectation q_L.

[Proved conditional closure] Integrating the first derivative of epsilon gives

    D(L)=epsilon(1)-epsilon(0)
        =q_*-integral_(0,1)(1-s)dmu(s).

Thus (14) implies D(L)<=C(d)M/L and A<=(4/pi)C(d)M. It asks that the response cancel an order-M bare seam cost to accuracy M/L. Local bond norms, finite-group representation theory and exact ground-sector projections do not provide this cancellation. Formula (14) is equivalent to the target after this exact identity, rather than a proved reduction to a weaker estimate. Its virtue here is to state the missing quantitative cancellation without the false on-site-symmetry premise.

[Proved conclusions and unresolved scope] The general C(d) remains unproved, and no counterexample to it was found. The spread construction yields A<=pi M in the boundary-matched quadratic setting, with the inherited sharp all-quasi-free theorem A<=2M remaining stronger. Exact on-site ground-preserving copy paths are impossible on the Heisenberg product vacuum, and every central two-copy fractional-swap profile has a positive L-independent cost there. The large-N cyclic logarithm also has a nonvanishing bulk cost per copy in an ordered example. Finite-group LSM does not remove either obstruction.

[Interpretation] The disconnected on-site stabilizer and the explicit fractional-swap potential are additional concrete obstructions to this proposed route. Their local-energy manifestation is related to lesson 2 of _Cd_target.md, but here a massless interacting vacuum already suffices, without a nearly critical spectator. Recasting the missing cancellation as (14) encounters lesson 4, the remaining reconnection estimate is the target itself. No crossed-spectrum assumption, sector-counting bound, crossover cutoff, or eventual convexity assumption was used, so lessons 1, 3 and 5 are not bypassed or contradicted.


## 10. Consolidated arbitrary-unitary obstruction and independent checks

[Proved] The stronger Ising theorem in [COPY_ROTATION_PROOF.md](COPY_ROTATION_PROOF.md), section 3, applies to every N>=2 and L>=3. For the symmetric cat ground state in each copy,

    inf_(U=tensor_x U_x) [<U Omega|H_sigma|U Omega>-N E_L]=JN,
    U_x in U(2^N).

This includes every spatially varying on-site unitary, not just fractional cyclic permutations. The actual permutation defect is zero. Regrouping the product cats gives 2^(-N/2) sum_a |a>^L. Its two-rung reduced states after U make the Ising expectation multilinear in W_x(b,a)=|<b|U_x|a>|^2. Relaxation to doubly stochastic matrices is tight because a multilinear minimum can be chosen at permutation matrices, which are allowed unitaries. For the resulting bijections f_x, the Hamming triangle inequality bounds each path from f_0(a) to P f_0(a) below by their direct distance. Averaging gives N/2 mismatches and therefore cost JN. Identity maps attain it. The minimum cost per copy remains J at large N, on compatible lengths. Selecting an aligned ordered ground state or coherently combining different unitary products evades this obstruction. In particular the inherited coherent seam trial costs 4J/(L+1) for the cats.

[Numerical] The independent script [copy_rotation_checks.py](copy_rotation_checks.py) constructs the signed Gaussian trial vectors explicitly in the physical spin Hilbert space, including both parities, pairing, and two modes per site. It also tests arbitrary rung unitaries in the Ising theorem and the Heisenberg correlation eigenspace. Its outputs are in copy_rotation_validation.json. The overlapping independent interacting energies and trial costs agree to 3.78e-14. These checks supplement the main numerical scripts and do not change the unresolved integrated seam lemma in section 9.
