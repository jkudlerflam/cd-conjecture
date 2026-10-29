# C115. Velocity-resolved infrared counting

2026-10-29.

## 1. Scope and result

[Assumed] Use the nearest-neighbor, one-site translation-invariant setting and the optimized bond half-spread M of ../../specs/_bound_background.md. Each individual chain has finitely many unitary infrared sectors. Write c_a for the nonchiral central charge and v_a for the mean of the positive right and left speeds. Set

    A = sum_a c_a v_a,
    N(u) = sum_{v_a >= u} c_a,
    V = 4 e M.

[Inherited, proved conditional on dynamical matching] C109 gives each chiral propagation speed at most V, hence v_a <= V. This needs a nonzero local real-time realization of the conformal fields. A Casimir fit alone is insufficient. The general free-fermion bound A <= 2M is inherited from REPORT 9.71 and contrib/fable_bound_freefermion, not rederived here. M=0 has no positive-speed sector under the dynamical matching assumption, so take M>0 below.

[Open] No general bound A <= C M f(d) is proved here. The new results are a counterexample to the unfiltered finite-size counting premise, an explicit correction to the tower inclusion, and precise sufficient velocity-window spectral estimates. The dynamical capacity argument still lacks a quantitative infrared-to-microscopic channel identification. No counterexample to the target bound is constructed.

## 2. The energy-window direction

[Assumed, finite-level spectral matching] On vacuum-compatible ring lengths, a state with conformal labels (h_a, hbar_a) has

    E(L)-E_0(L) = (2 pi/L) sum_a [v_R,a h_a+v_L,a hbar_a]+o(L^-1).

For equal chiral speeds and a state supported in one sector, put Delta=h+hbar. Then

    E(L)-E_0(L) = 2 pi v_a Delta/L+o(L^-1).

[Proved from this matching] In a window E-E_0 <= 2 pi u x/L, the accessible dimensions obey Delta <= u x/v_a, with a small slack at a threshold to absorb the remainder. Thus v_a >= u does not imply inclusion of the tower through dimension x. It gives the opposite inequality. For example, v_a=2u puts the level Delta=x at twice the proposed window energy. The correct inclusion through x holds for v_a<=u. For all fast sectors together, the velocity ceiling only guarantees inclusion through dimension u x/V. One can instead enlarge the energy window to 2 pi V x/L, but slow sectors remain in it.

[Proved] With unequal speeds the window is an inequality on v_R h+v_L hbar, not on the mean velocity times h+hbar. An average-speed threshold cannot select both chiral towers by energy alone. None of the negative results below requires this complication, since their chains have equal chiral speeds.

[Inherited, proved] One Virasoro stress tensor is not c independent states. Its level-two vacuum state has norm squared c/2, whereas the dimension of its span is one. At any fixed descendant level, the number of Virasoro monomials is bounded independently of c. C106 Section 6 proves this obstruction. Counting the vacuum towers of sectors therefore does not, by itself, count their central charges. Counting a direct sum of towers also produces a sum of degeneracies, whereas an extensive lower bound on log of the count generally needs independent products of states.

## 3. An unfiltered ring trace cannot have the proposed bound

### 3.1 A fixed-d, fixed-M family with a fast sector

[Definition] Put two spin-1/2 factors on each site. On each factor use the XX bond

    h_XX(t) = (t/2)(X tensor X+Y tensor Y),   t>0.

Take their sum on the d=4 site and choose

    t_f=1/(1+epsilon),
    t_s=epsilon/(1+epsilon),   0<epsilon<1/3.

[Proved] A single bond has eigenvalues {-t,0,0,t}. The two factors commute and their extremal eigenvalues add, so the bare half-spread is t_f+t_s=1. This is also the optimized M. Indeed physical exchange of the two entire sites fixes the bond and negates every coboundary o tensor 1-1 tensor o, including o which couples the two factors. Convexity and swap invariance of half-spread show that averaging a representative with its swap cannot increase its half-spread and returns the original bond up to a scalar. Thus no gauge decreases M.

[Proved, exact free infrared] Each factor has dispersion 2t cos k, two Fermi points, c=1 and v=2t. Consequently

    d=4,   M=1,
    v_f=2/(1+epsilon)>1,
    v_s=2epsilon/(1+epsilon)<1,
    N(1)=1,
    A=v_f+v_s=2=2M.

The exact full counting function is 2 below v_s, 1 between v_s and v_f, and 0 above v_f, with the stated >= convention at endpoints. This family saturates, rather than challenges, the free bound.

### 3.2 Divergence after the infrared limit

[Definition] Fix tau>0 and u>0, and define the excitation trace

    Z_L(tau,u) = Tr exp[-(tau L/u)(H_L-E_0(L))].

This uses inverse temperature beta=tau L/u. The expression Tr exp[-tau L(H-E_0)] in the question is recovered by fixing energy units or absorbing 1/u into tau.

[Proved using finite-level free matching] The slow XX sector contains the oscillator states (J_-1)^m|0>, m=0,1,..., with dimensions m. Hold the fast sector in its vacuum. For every fixed integer K and every fixed epsilon>0, there are K+1 orthogonal ring eigenstates approaching these labels on compatible lengths. Therefore

    liminf_(L->infinity) Z_L(tau,u)
       >= sum_(m=0)^K exp[-2 pi tau v_s m/u].

Take K to infinity only after the L limit. Positivity of every summand gives

    liminf_(L->infinity) Z_L(tau,u)
       >= 1/[1-exp(-2 pi tau v_s/u)].                 (3.1)

At fixed tau and u=1 the right side diverges as 1/(2 pi tau v_s) as epsilon tends to zero. This proof requires convergence of finitely many levels at a time, not exchange of an infinite trace with the scaling limit. A complete compact-boson partition function would give a stronger divergence, but is unnecessary.

[Proved using the same states] Let D_L(u,x) be the number of eigenstates with E-E_0<=2 pi u x/L, for x>0. All integers m>=0 with v_s m<u x occur in the window for sufficiently large L. Hence

    liminf_(L->infinity) D_L(u,x)
       >= number of integers m>=0 with m<u x/v_s.    (3.2)

This number diverges as epsilon tends to zero. The strict inequality avoids a finite-size endpoint assumption. Both (3.1) and (3.2) hold with d=4, M/u=1 and N(u)=1. Thus no finite function of d, M/u and fixed positive tau or x can bound these unfiltered quantities uniformly over the class, even if the bound is only imposed after taking L to infinity. Ground-state degeneracy is not the cause. Compatible even XX rings have a unique vacuum.

[Proved scope] These are counterexamples to an upper bound on the total number of low-energy states or the total positive excitation trace. They do not exclude a functional that removes slow states. They also explain why a lower bound on energies of "most states" is insufficient unless it controls the exceptional states at the proposed 1/L scale. A small fraction of the microscopic Hilbert space can already violate a uniform infrared counting estimate.

## 4. The spatial transfer matrix has no fixed d^2 dimension

[Proved tensor accounting] On an even ring split H into its even and odd bonds. A finite Trotter approximation to exp(-beta H) is

    [exp(-beta H_even/n) exp(-beta H_odd/n)]^n.

Each gate acts on d by d sites and admits an operator Schmidt decomposition with at most d^2 terms. This is a bound on one gate. Multiplying n time layers crosses a spatial cut n times, and the uncompressed bond index has dimension at most (d^2)^n, not d^2. Equivalently, the spatial transfer matrix acts on an entire imaginary-time column. Its dimension grows with the Trotter number. Taking n to infinity is part of defining the generic noncommuting quantum partition function.

[Proved limitation] Local gate dimension is not a proof of a uniform exact or approximate compression of this column. Such a compression would need its own error and low-temperature bounds. At beta of order L/u these controls would need to remain adequate as L grows. Classical commuting nearest-neighbor interactions can have a fixed finite transfer matrix, but generic quantum interactions do not inherit that construction. Even a hypothetical compression of some special chain would still have to accommodate the trace divergence of Section 3.

[Source context] The quantum transfer matrix and its infinite-Trotter-number limit are explicitly distinguished in Faulmann, Goehmann and Kozlowski, [Low-temperature spectrum of the quantum transfer matrix of the XXZ chain in the massless regime](https://arxiv.org/abs/2305.06679). The tensor-index argument above does not use integrability or a theorem from that paper.

## 5. A sufficient counting lemma and its explicit constants

### 5.1 Dyadic windows

[Definition] Set s_j=V 2^(-j) and

    I_j = {a: s_j/2 < v_a <= s_j},
    C_j = sum_(a in I_j) c_a,   j=0,1,... .

These disjoint windows cover every positive speed once. They also respect equality at the upper endpoint.

[Proved] Integrating each indicator gives the identity

    A = integral_0^V N(u) du.

Moreover,

    (1/2) sum_j s_j C_j <= A <= sum_j s_j C_j.       (5.1)

Thus a weighted summability bound on C_j is the relevant information. An unweighted bound on the number of sectors is unnecessary.

[Assumed sufficient lemma W_gamma] For some fixed 0<=gamma<1 and an explicit function K(d), every chain in the task satisfies

    C_j <= K(d) 2^(gamma j),   for all j>=0.          (5.2)

The constants are independent of the Hamiltonian, minimum speed, ring length and crossover scale. No internal symmetry is required in this statement. Establishing (5.2), rather than merely fitting it on examples, is the missing theorem.

[Proved conditional bound] Substitute (5.2) into (5.1) and sum the geometric series:

    A <= V K(d) sum_(j>=0) 2^[-(1-gamma)j]
      = [4e/(1-2^(gamma-1))] M K(d).                (5.3)

For the concrete choice gamma=1/2 this is

    A <= 4e(2+sqrt(2)) M K(d).                       (5.4)

This is an explicit constant multiplying the still-unproved function K(d), not a completed bound for all chains.

[Proved conditional tail] Put J=floor(log_2(V/u)). Every sector with v_a>=u belongs to one of I_0,...,I_J, although part of I_J may lie below u. For gamma>0,

    N(u) <= K(d) [2^(gamma(J+1))-1]/[2^gamma-1]
          <= [2^gamma/(2^gamma-1)] K(d) (V/u)^gamma. (5.5)

This gives the requested integrable power with delta=1-gamma. For gamma=0 it gives N(u)<=K(d)[1+log_2(V/u)]. Direct summation in (5.3) gives a better constant than integrating this rounded upper envelope.

[Assumed weaker sufficient lemma; proved implication] A bound C_j<=K(d)(j+1) instead yields

    N(u) <= (K(d)/2)(J+1)(J+2),
    A <= 4 V K(d) = 16e M K(d).                     (5.6)

A squared logarithm is integrable too. The identity sum_(j>=0)(j+1)2^(-j)=4 proves the last line. This variant accommodates many comparable sectors concentrated at one extremely slow scale. A uniform bound on C_j, gamma=0, is stronger and is not suggested by those constructions.

[Proved direct-tail constants] If N(u)<=g(d)(M/u)^(1-delta), 0<delta<=1, direct integration gives A<=[(4e)^delta/delta] M g(d). If N(u)<=g(d)[1+log(V/u)], it gives A<=8e M g(d). The borderline estimate N(u)<=g(d)M/u does not integrate uniformly at zero.

### 5.2 A spectral realization which measures c rather than tower number

[Assumed for this construction only] Suppose a velocity window has a compact unitary full CFT Hilbert space with central charge C_j and an S-invariant torus partition function. This is automatic for a tensor product of individually modular-invariant full sectors. It is not automatic for a chiral vacuum module, a projected subspace, or an arbitrary collection of stress tensors inside a nonfactorized modular invariant. In a general chain the corresponding spectral extraction or completion must be supplied as part of a new lemma.

[Definition] Let D_j be the sum of h+hbar on this full window Hilbert space, with vacuum eigenvalue zero. Define the dimensionless trace

    Q_j(t)=Tr exp(-2 pi t D_j),
    Z_j(t)=exp(pi C_j t/6) Q_j(t).

Here D_j is a conformal-dimension operator after the infrared limit, not the microscopic Hamiltonian divided by one common speed.

[Proved within these assumptions] S-invariance gives Z_j(t)=Z_j(1/t). The vacuum term alone gives Z_j(1/t)>=exp(pi C_j/(6t)). Hence

    log Q_j(t) >= (pi C_j/6)(1/t-t),
    log Q_j(1/2) >= pi C_j/4.                       (5.7)

This avoids a large-c or Cardy asymptotic estimate. It uses full state content and modular invariance. A single vacuum Virasoro character has no such guarantee. The partition-function conventions and S transformation can be checked in Dey, Pal and Qiao, [A universal inequality on the unitary 2D CFT partition function, Section 2.1](https://arxiv.org/html/2410.18174v2). Only the vacuum term and S transformation are used here, not their sparse-spectrum estimates.

[Assumed sufficient spectral lemma S] Construct, for each window in every microscopic chain, a positive infrared spectral functional Q_j satisfying both

    log Q_j >= pi C_j/4,
    log Q_j <= B(d) 2^(j/2),                        (5.8)

with an explicit B(d) independent of the chain and crossover. For tensor-product full CFT sectors, (5.7) supplies the first inequality with Q_j=Q_j(1/2). The second is the new lattice input. For general sectors, even obtaining the first inequality from an admissible lattice spectral selection is an additional part of the lemma, not an assumption silently included in the task.

[Proved conditional implication] Equations (5.8) imply W_(1/2) with K(d)=4B(d)/pi. Therefore

    A <= [16e(2+sqrt(2))/pi] M B(d).                (5.9)

An upper bound B(d)(j+1) instead gives A<=64e M B(d)/pi. These are explicit crossover-free implications. No estimate of B(d) for the general interacting class is established here.

[Proved relation to a physical ring trace, equal speeds] For sectors in I_j, take beta=L/(2s_j) and hold other sectors in their vacua. Their individual aspect ratios t_a=v_a/(2s_j) lie in (1/4,1/2]. Each excitation weight is at least its weight at t=1/2. Thus, assuming the full tensor-product matching, the projected physical trace has log at least pi C_j/4. An upper bound on that projected trace would also suffice. The unfiltered trace fails by Section 3. No microscopic projector with the required properties is supplied by locality alone.

[Proved obstruction to a simple momentum filter] A cutoff on the ratio of total energy to total momentum does not supply that projector. In one slow boson sector, take oscillator levels h=N+1 and hbar=N. Their conformal momentum is one, while their energy in the same units is v_s(2N+1). The apparent ratio can exceed any proposed fast threshold, although both chiral constituents belong to the slow sector. Equal left and right levels even have zero total momentum. Velocity selection must resolve excitation branches or independent sector generators, rather than identify a many-body state's E/|P| with a propagation speed.

[Proved lack of a reflection assumption] The abstract dimension trace treats unequal chiral velocities without identifying them. For one isolated factor, write rescaled energy e=L Delta E/(2 pi) and unwrapped conformal momentum p=L Delta P/(2 pi). Then

    e=v_R h+v_L hbar,   p=h-hbar,
    h=(e+v_L p)/(v_R+v_L),
    hbar=(e-v_R p)/(v_R+v_L),
    h+hbar=[2e+(v_L-v_R)p]/(v_R+v_L).              (5.10)

This requires sector-resolved conformal momentum, including the correct reference momentum. Total lattice energy and momentum alone cannot in general resolve many factors. Weighting by (5.10) is not a fixed-temperature Gibbs state. Alternatively beta=L/(4s_j) gives each chiral aspect ratio at most 1/2 for an average speed in I_j, and the two-variable modular vacuum bound still lower-bounds log of a full factor trace by pi c_a/4. A very small speed of one chirality can make this physical trace much larger. A uniform upper estimate is therefore a further issue, not repaired by averaging speeds.

[Proved qualification] Even a full dimension trace is sensitive to compactification radii and light primary multiplicities beyond c. For a compact boson, letting its radius tend to infinity makes Q(1/2) unbounded while c=1. Therefore S is stronger than W, and a universal upper bound on the full Q is not a consequence of bounding c alone. A useful lattice proof might have to select a smaller set of states which still satisfies the lower bound in (5.8). No general such selection is constructed here.

## 6. Dynamical information capacity

[Inherited, proved] Lieb-Robinson bounds constrain the support and norm of commutators outside an upper propagation cone. They do not give a lower bound on the number of operationally distinguishable channels, their amplitude, or their bandwidth inside that cone. Multiplying an infrared stress field by a small microscopic overlap does not alter its conformal central charge or propagation speed. C106 also explains why overlapping with every stress field is weaker than independent access to them.

[Inherited, proved] For evolution across one nearest-neighbor cut, with arbitrary ancillas on either side, small incremental entangling gives

    |dS/dt| <= 18 M log d.                          (6.1)

This bound uses a centered gauge representative and then its infimum. A finite ring bipartition has two cuts and correspondingly two contributions. Local operations within either side do not contribute to the entanglement derivative. The ancilla-independent estimate is established by Van Acoleyen, Marien and Verstraete, [Entanglement rates and area laws](https://arxiv.org/abs/1304.5931).

[Assumed packet premise; proved consequence] Suppose C_w=sum_(w<v_a<=2w)c_a can be converted into an entropy gain at least eta C_w across one cut in time at most b ell/w, with constants eta,b>0 and packet length ell in microscopic sites. Integrating (6.1) gives

    eta C_w <= 18 b M ell log d/w,
    C_w <= (18b/eta)(M/w) ell log d.                (6.2)

Even ell bounded by a constant gives only a 1/w window bound. Summing s_j C_j over infinitely many windows then diverges. Reusing (6.1) independently at different scales does not allocate a single finite capacity budget across those scales.

[Proved analogous space-time counting] During time t a Lieb-Robinson cone contains O(Mt) sites per side. Their Hilbert-space entropy capacity is O(Mt log d). Allowing a wavepacket flight time t of order ell/u therefore gives capacity of order M ell log d/u, with the same unfavorable scaling. A light cone bounds the number of sites available, but does not force one distinct microscopic site per conformal channel.

[Proved obstruction for the proposed IR packet mechanism] In C106, d=4 and v/M=r>0 stay fixed while the two Fermi pockets narrow on momentum scale epsilon/t. To remain in a fixed-accuracy linear neighborhood of an individual crossing, a packet must have momentum width O(epsilon/t). The Fourier uncertainty relation then requires spatial width at least of order t/epsilon. This statement concerns packets confined to those linear neighborhoods. It does not prove that every possible communication protocol is slow. It does refute deriving a uniform microscopic IR packet width from the velocity alone, which is the extra premise needed by that direct construction.

[Assumed sufficient dynamical lemma, explicit missing input] A simultaneous one-cut protocol with all selected sectors, whose entropy rate obeys

    dS/dt >= eta A/ell_0,
    ell_0 <= F(d),   eta>0 universal,

would give A<=18 F(d) M log d/eta by (6.1). Preparing the state, establishing independence of the channels and controlling ell_0 are part of this lemma. The packet construction just discussed does not prove it. Merely postulating it in a dynamical form does not remove C69/C99's preparation-length obstruction. A new multi-scale capacity inequality that avoids such a microscopic IR packet width could still work, but is not obtained here.

## 7. Tests and limitations of the sufficient lemmas

### 7.1 General nearest-neighbor free fermions

[Inherited, proved] For d=2^m and a pure conformal quasi-free infrared without a flat zero band, the free counting theorem gives c_total<=m. One way to see its degree count is that for an m by m number-conserving Laurent symbol h(z)=B^dag/z+X+Bz, z^m det h(z) is a polynomial of degree at most 2m. Linear crossings counted with multiplicity contribute one half each to c_total. In a 2m by 2m Majorana symbol the polynomial degree is at most 4m and each crossing contributes one quarter. A nonzero flat zero band is excluded by the stipulated finite set of positive-speed CFT sectors. See contrib/fable_bound_freefermion Section 3 for the stronger rank count.

[Proved consequence of the inherited bounds]

    N(u) <= min(log_2 d, 2M/u).                      (7.1)

The second inequality is u N(u)<=A<=2M. Thus free chains obey W_gamma, for example with K(d)=log_2 d and gamma=1/2, while their known direct A<=2M theorem is much stronger than (5.4). This verifies the central-charge counting lemma in the free class. It does not by itself verify every proposed full-trace upper bound S.

### 7.2 Stacks of XX chains

[Proved] For m independent XX factors with positive hoppings t_a,

    d=2^m,   M=sum_a t_a,
    c_a=1,   v_a=2t_a,
    N(u)=number of a with 2t_a>=u,
    A=2 sum_a t_a=2M.                               (7.2)

The swap argument of Section 3 proves gauge optimality for arbitrary m and arbitrary on-site coboundaries. Formula (7.2) checks the layer-cake identity, additivity and velocity thresholds. Equal hoppings have N(2M/m)=m. Therefore a universal power tail N(u)<=g(d)(M/u)^(1-delta) must allow

    g(2^m) >= 2^(1-delta) m^delta.                  (7.3)

A proposed logarithmic tail N(u)<=g(d)[1+log(V/u)] must allow g(2^m)>=m/[1+log(2em)]. These are necessary tests on candidate functions, not violations of the relaxed target.

[Proved compact-boson check] In the convention

    Delta_(n,w)=n^2/(2R^2)+w^2 R^2/2+N+Nbar,

the full dimension trace is

    Q_R(t)=product_(k>=1)(1-exp(-2 pi t k))^(-2)
            times sum_n exp(-pi t n^2/R^2)
            times sum_w exp(-pi t w^2 R^2).          (7.4)

Poisson summation and the oscillator product transformation give the modular relation for exp(pi t/6) Q_R(t). At t=1/2, log Q_R>=pi/4 as in (5.7). At the XX radius R=sqrt(2) in this convention, log Q_R is a fixed finite number b_XX. Independent XX windows therefore have log Q_j=b_XX C_j<=b_XX log_2 d. The spectral lemma S holds for these stacks with B(d)=b_XX log_2 d. Equation (7.4) also shows the radius qualification in Section 5: as R grows, its first theta factor grows proportionally to R.

### 7.3 The C106 fixed-speed crossover family

[Inherited exact example; proved test] C106 has two c=1 sectors at v=rM with 0<r<1 fixed, independent of the small avoided-crossing parameter. Thus

    N(u)=2 if 0<u<=rM, and 0 if u>rM,
    A/M=2r.

The velocity window containing rM and its C_j are unchanged as the thermal matching scale collapses. Taking L to infinity separately for each member does not require a uniform finite L of CFT validity. Thus W_gamma is compatible with the C106 counterexample. A proof that instead selects a uniformly bounded ring length from r and d would reintroduce the excluded crossover premise.

### 7.4 Campaign slow-mode constructions

[Assumed campaign infrared interpretation] The fixed-alphabet constructions have n payload bosons, and in the mobile version one additional length mode. Their general lattice-to-CFT transfer remains conditional. The following tests retain that qualification. Velocities must first be expressed per physical site: if a fixed cell has b sites and a logical-cell group velocity v_cell, then v_physical=b v_cell, because k_cell=b k_physical. This factor must not be omitted.

[Proved conditional envelope test] Suppose the payload sector speeds of the n-th construction satisfy

    max_a v_(a,n)/M_n <= K exp(-kappa n),
    c_payload=n,

with K,kappa>0 independent of n. For any fixed u=r M_n, a payload can be counted only if n<=log(K/r)/kappa. Consequently, uniformly over this family,

    N_n(u) <= c_length + max(0,log(K M_n/u)/kappa), (7.5)

where c_length=0 for frozen cells and 1 for the mobile length mode. If n exceeds the threshold, no payload is counted. At u of order M this bound is independent of n, even though c_total grows without bound. Also A_payload/M_n<=K n exp(-kappa n), with maximum over positive real n equal to K/(e kappa). Polynomial cell-length factors can be absorbed into a smaller positive exponential rate when the cell length grows polynomially.

[Proved distinction] An exponential decrease of the smallest velocity alone does not prove (7.5). One needs an upper bound on the largest payload speed or a bound on the ordered speed distribution. If instead the ordered sector speeds within each chain satisfy v_a/M<=K exp(-kappa a), direct counting gives the same logarithmic envelope. These are different assumptions. A statement about v_min must not be substituted for either.

[Proved gauge lower bound for this construction] The vacuum-token hopping move is X_(v,a) R_b -> L_a X_(v,b), with matrix element -t, as specified in C22 RESULT Section 2 and retained in routeA_spatial/ycell.py, bond_moves. Both site letters change. For orthogonal basis letters i!=k and j!=l,

    <ij|[alpha 1+o tensor 1-1 tensor o]|kl>=0.

Thus that matrix element is -t in every gauge representative and after any scalar centering. Every centered representative has norm at least the magnitude of any matrix element, so M>=t. This supplies the needed lower bound on the optimized scale without equating M to the bare norm.

[Inherited numerical evidence with proved conversion, not an all-n theorem] contrib/margins_audit/README.md records both v_min and v_max, not only v_min. For n=1,2,3,4,23 its rounded maximum-velocity logarithms are -13.8,-20.8,-29.5,-38.1,-205.2, in units lambda^2 g^2 alpha_n^4 t per cell. The reported ratio v_max/v_min is at most 10^2.04 over the tabulated n<=23. Restoring b=n+1 and using M>=t gives the exact conversion v_max,physical/M<=b v_max,cell/t. For n=23 the reported right side is approximately

    24 times 10^(-205.2) lambda^2 g^2 alpha_n^4.

For lambda,g,alpha_n<=1, the five quoted n values all have leading payload velocities below 10^(-13) M at the precision reported. They contribute N(M)=0 in their frozen effective infrared description. The printed decimals are not rigorous enclosures. This uses the recorded leading-order velocities and the campaign's perturbative/infrared interpretation; it does not certify the exact interacting spectrum or an all-n exponential bound. A mobile length mode contributes at most one additional unit of central charge at the threshold.

[Inherited heuristic mechanism] contrib/routeA_spatial/README.md gives express-transfer suppression approximately (t_x/D_x)^(n+1) and J_n proportional to lambda^2 g^2 alpha_n^4 times its square. The hierarchy permits further shrinking lambda at each n. The assembly notes give physical payload speeds of order b W and place the length-fluid bandwidth at a comparable small scale. These are the source of the conditional envelope test, not a derivation of a new uniform microscopic estimate. The slow constructions pass the intended check at u of order M under their recorded normalized hierarchy.

[Proved caution for stronger window lemmas] Many comparable payload speeds may occupy only a few dyadic windows even when n is large. If their ratio is uniformly bounded by R, at most ceil(log_2 R)+1 windows are occupied, so one contains at least n/[ceil(log_2 R)+1] central charge. Therefore an all-n comparable-speed realization at fixed d would rule out the stronger uniform C_j<=K(d), while remaining compatible with (5.2) for gamma>0 and with the linear-j bound when the occupied j grows at least proportionally to n. The current finite-n velocity table is evidence for this possibility, not its proof.

## 8. Verification and final status

[Numerical] checks.py performs only small deterministic bond checks, scalar sums, and arithmetic on the documented velocity maxima. It checks the d=4 optimized bond witness, the tower-window direction, growth of the lower trace bound, compact-boson modular normalization, the geometric-series constants, exact XX layer-cake integration and a specified exponential-envelope model. It does not run many-body infrared extrapolations or re-estimate the campaign's spectra. CHECKS.md contains its output.

[Proved consistency checks] Positive rescaling of H multiplies M, velocities and u by the same factor; N at fixed u/M and every dimensionless spectral functional are unchanged. A and the right sides of (5.3), (5.9) scale linearly. The fast-plus-slow example has A=2M for every epsilon. Endpoints of the windows are assigned once. The order L to infinity, then the number of retained oscillator states to infinity, then epsilon to zero is explicit in (3.1). No thermal crossover is presumed in that proof.

[Open final conclusion] The general task reduces to a new velocity-resolved infrared counting inequality such as W_(1/2), or a spectral construction satisfying both parts of S. Locality as used here supplies neither. Raw finite-size traces fail before this step, and the dynamical capacity estimate supplies neither a uniform IR bandwidth nor a summable budget across velocities. The target A<=C M f(d) remains open.

[Proved work record] All created files are inside work_C115. No heavy numerics, SSH command, SLURM job, network shell command, commit or push was needed. Literature access used the browsing tool only. Project STATE.md and REPORT.md were left unchanged because writes outside this folder are prohibited; local STATE.md supplies the handoff.
