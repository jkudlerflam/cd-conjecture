# C127 proof and limitations

2026-11-06. All logarithms are natural. A status tag applies to the paragraph and its displayed formulas. The microscopic setting, gauge-optimized M, and even-ring conventions are those of the shared target. Results inherited from C125 are cited without repeating its transfer construction.

## 1. What is being compared

[Proved, inherited from C125] For a positive spatial transfer operator, put Hc=-log T-epsilon_0 and g(L)=log Tr exp(-L Hc). The support of T is understood. If both traces exist,

    R_L=2g(L)-g(2L),
    g(L)<=R_L<=2g(L).                                           (1)

The first inequality uses g(2L)<=g(L); the second uses g(2L)>=0. The normalized trace includes all crossed states, including any multiplicity of the largest eigenvalue. In the checkerboard construction the physical lengths are even. Its spectral compactness limit reproduces every even L>=4. XX has an explicit limit also at L=2. These are the length conventions used below. The XX and antiferromagnetic spin-1/2 examples have no forced defect on this even sequence. A proof for another model must specify its own vacuum-compatible sequence as well as the lengths on which its transfer construction applies.

[Proved] The following are different statements.

1. The range of a fixed-frequency single-insertion map has small dimension.
2. The invariant subspace detected by that map has small dimension.
3. There is an energy- and momentum-preserving comparison with a finite-species Fock space.
4. There is only a scalar bound on Tr exp(-L Hc).

Statement 1 holds. Statement 2 fails already for XX. Statement 3 with a global linear dispersion fails for XX, even with any fixed multiplicative constant. Statement 4 holds for XX with explicit constants, as proved in section 5. The general statement 4 remains unresolved here.

[Proved] The distinction matters for the target. In a transfer representation reproducing the physical traces, the largest eigenvalue determines the spatial pressure p(beta). Thus

    g(L)=log Z(L,beta)-L p(beta).

Taking beta to infinity at fixed allowed L gives g(L)/beta=e_0 L-E_L. Therefore g(L)<=K+kappa M beta/L implies E_L>=e_0 L-kappa M/L. Under the stated infrared expansion this gives A<=(6/pi)kappa M. Equation (1) also gives the all-temperature ratio bound with constants 2K and 2kappa, and hence C(d)<=2kappa on the selected lengths. No interchange of the large-L and large-beta limits is required for these deductions.

## 2. The insertion space and the detected spectral space

[Assumed representation for this section] Let the reflected column construction have a Hilbert space, a nonnegative Hc, a commuting unitary circle action U(tau)=exp(i tau Pc), and a rotation-invariant vacuum vector Omega. Let J_tau(O) be the vector made by an insertion of the one-site operator O on the chosen reflected boundary. The reflection operation is included in the inner product. It is not legitimate to identify an arbitrary unreflected A,B correlator with a positive spectral measure. At a finite Trotter number the analogous objects use a discrete circle. C125's eigenvalue compactness argument alone does not construct these insertion maps or a continuous circle action in the limit.

[Proved under this representation] At omega_k=2pi k/beta define

    F_k(O)=beta^(-1/2) integral_0^beta exp(-i omega_k tau) J_tau(O) d tau,
    W_k=range F_k,
    C_k=closure span { f(Hc)F_k(O): f bounded Borel, O in M_d }.
                                                                    (2)

The overall normalization can be changed without changing either space. Time covariance puts W_k in the momentum-omega_k eigenspace, with the Fourier sign fixed by the displayed convention. For k!=0, the identity maps to zero. Hence

    dim W_k<=d^2-1,
    dim(W_k+W_(-k))<=2(d^2-1).                                  (3)

At zero frequency remove the vacuum component before counting connected insertions. If multiple vacua are treated separately their detected spaces must also be included. The single-vacuum count does not bound vacuum degeneracy.

[Proved] The second space in (2) is the spectral cyclic span. In general Hc W_k is not contained in W_k. In a discrete spectrum, write

    F_k(O)=sum_E P_(E,k) F_k(O).

For each fixed E the rank of O -> P_(E,k)F_k(O) is at most d^2-1. There can nevertheless be arbitrarily many distinct E with a nonzero image. Even a single vector with nonzero components at infinitely many distinct eigenvalues has infinite-dimensional cyclic span. A rank bound on the matrix of two-point functions at a fixed separation is not a bound on the number of exponential decay rates in that matrix.

### 2.1 The energy statement that follows from C125

[Proved, using the C109/C125 commutator estimate] For disjoint one-site A,B separated by x, positive bosonic omega, and V=4eM,

    G_AB(i omega,x)=i integral_0^infinity exp(-omega t)<[A(t),B]> dt,
    ||[A(t),B]||<=2||A||||B|| min(1,exp(-x+Vt)).                 (4)

Writing C=2||A||||B|| and splitting the integral at x/V gives, when 0<omega<V,

    |G_AB| <= C{[exp(-omega x/V)-exp(-x)]/(V-omega)
                              +exp(-omega x/V)/omega}
           <= C V/[omega(V-omega)] exp(-omega x/V).             (5)

At omega=V the first term is C x exp(-x)/V and the second is C exp(-x)/V. The divergent prefactor in (5) must not be used at that endpoint. Negative frequencies follow by conjugation and exchange of the insertions. At omega=0 this argument supplies no positive gap.

[Proved conditional spectral implication] Suppose the reflected diagonal correlator, possibly with a fixed spatial offset x0 imposed by the reflection plane, is

    C_O(x)=<F_k(O), exp(-x Hc) F_k(O)>
          =integral_[0,infinity) exp(-x E) d mu_(O,k)(E),
    mu_(O,k)>=0.                                               (6)

An offset changes the measure by exp(-x0 E) and does not change its support at finite energy. If mu([0,a])=w>0 for a<omega/V, then C_O(x)>=w exp(-a x), contradicting (5) along arbitrarily large allowed spatial separations. Consequently

    P_[0,|omega_k|/V) F_k=0,
    Hc restricted to C_k >= |omega_k|/V,   0<|omega_k|<V.        (7)

This is a support statement for detected energies, not an assertion that vectors in W_k are eigenvectors. Polarization and the span over O extend it to C_k. It says nothing about orthogonal, undetected states at the same momentum.

[Proved, inherited optimized rate] Retaining the LR parameter mu>0 gives the asymptotic decay rate min(mu,mu omega/(4M exp(mu))). Its optimum is

    alpha(omega)=omega/(4eM),                    0<=omega<=4eM,
                 log(omega/(4M)),               omega>=4eM.    (8)

At the optimized matching point a linear factor in x is allowed. Such a factor still excludes spectral support strictly below alpha. Thus, under (6), Hc>=alpha(|omega|) on C_k. Only the low-frequency branch of (8) is linear. These continuum statements are not estimates for a checkerboard with a fixed nonzero Trotter step.

## 3. A lattice counterexample to counting detected energies by insertion rank

[Proved] Use the actual XX spin chain with h=Sx Sx+Sy Sy, M=1/2. C125 section 3 gives the one-fermion resolvent

    G_x(i nu)=(-i sign nu)/sqrt(1+nu^2)
                [-i sign nu exp(-asinh|nu|)]^|x|,
    nu=(2n+1)pi/beta.                                         (9)

The density Sz is an original one-site spin operator. Its connected Matsubara correlator is the fermion bubble

    chi(i omega_k,x)=-(1/beta) sum_n G_x(i nu_n)G_(-x)(i nu_n-i omega_k).
                                                                    (10)

[Proved] The sign of the reflected correlator is worth keeping explicit. For odd x=2r+1, the phase in each product of (9) is

    [-sign(nu_n)sign(nu_n-omega_k)]^(x+1)=1.

Therefore -chi(i omega_k,2r+1) is the positive sum

    sum_n w_(n,k) exp(-2r E_(n,k)),
    E_(n,k)=asinh|nu_n|+asinh|nu_n-omega_k|,
    w_(n,k)=exp(-E_(n,k)) /
              [beta sqrt((1+nu_n^2)(1+(nu_n-omega_k)^2))]>0.     (11)

This is the bond-reflected density correlator, including the sign reversal of Sz in the staggered reflection convention. The two-site propagation is exp(-2Hc). Using even separation without the reflection sign would instead produce signed coefficients and would not prove the positive-measure statement.

[Proved] Fix any k>=1. For n>=k, both nu_n and nu_n-omega_k are positive, and E_(n,k) is strictly increasing with n and tends to infinity. Equation (11) has a positive mass at every one of these distinct energies. Its cyclic span is infinite-dimensional, although F_k(Sz) is one vector and dim W_k<=3 for the whole one-site algebra. This is a lattice-realizable counterexample, not an abstract spectrum.

[Proved] The obstruction occurs even in a fixed low-energy window as beta grows. For k=1 and n/beta held fixed,

    E_(n,1)=2 asinh(2pi n/beta)+O(beta^(-2)).

Consequently the number of this particular sequence of detected energies below a fixed E>0 grows as

    beta sinh(E/2)/(2pi)+O(1).                                 (12)

For small E it is beta E/(4pi)+O(beta E^3)+O(1). No dimension-only bound on the number of detected energies follows by restricting to E<1. This growth does not violate a bosonic many-body count. Even one compact boson has many multiparticle states at a fixed total momentum. It refutes only the use of (3) as a count of spectral species or eigenstates.

[Numerical] The script enumerates 500 distinct positive masses at k=1 for beta=10,100,1000. The numbers below E=1 are respectively 0,8,82. At beta=100 the first eight energies run from 0.1255195673 to 0.9668165104. The residue formula (9) was independently checked by momentum quadrature, with maximum residual 5.27e-16. The reflected sign and the bound (5), with ||Sz||=1/2, pass all sampled bubble checks. The bubble evaluation has a stated finite frequency cutoff and is a numerical check, while (11) proves the infinite-support assertion.

### 3.1 A different lattice failure, momentum-preserving linear majorization

[Proved, using C125's XX states] The same chain has states at p_k=2pi k/beta with

    E_k=asinh(pi/beta)+asinh((2k-1)pi/beta).                     (13)

For any finite V' and fixed beta, E_k<|p_k|/V' for all sufficiently large k. Every bosonic Fock configuration with integer Matsubara momenta and energies at least |omega|/V' has total energy at least |p_k|/V' in that momentum sector, by the triangle inequality. With finitely many species its sector trace has large-L behavior

    Z_Fock,k(L)=D_k exp(-L |p_k|/V')(1+o(1)),

where D_k is finite. The actual XX sector trace is at least exp(-L E_k). Their ratio diverges along even L. No fixed exp(K) restores that momentum-sector comparison. A separate p=0 factor that has a finite trace and finite ground multiplicity also cannot change the large-L exponential rate.

[Proved limitation] This is not a counterexample to an unlabelled scalar trace bound. In a scalar comparison low-energy reference states at other momenta can supply the needed trace. It also does not rule out a comparison confined to sufficiently low energies, with the ultraviolet part estimated separately. The positive scalar result below makes the distinction concrete.

## 4. What a Fock majorant would prove

[Proved] Let

    B_S(a)=-S sum_(k>=1) log(1-exp(-a k)),
    a=2pi L/(V' beta).

The function f(x)=-log(1-exp(-x)) is positive and decreasing. Right-endpoint quadrature gives

    sum_(k>=1) f(a k) <= (1/a) integral_0^infinity f(x) dx
                       = pi^2/(6a),
    B_S(a)<=S pi V' beta/(12 L).                               (14)

To evaluate the integral, expand f(x)=sum_(n>=1)exp(-nx)/n and integrate the nonnegative terms. Thus there is a uniform upper bound with no positive lower-order correction. In particular a proved majorization g(L)<=K+B_S(a) would give A<=S V'/2.

[Proved normalization] S counts both momentum signs. A single c=1 boson of velocity v has S=2 and V'=v, so the leading oscillator log trace is pi v beta/(6L). It is not an equality at all aspect ratios. The oscillator asymptotic is

    B_S(a)=S[pi^2/(6a)+(1/2)log(a/(2pi))-a/24]
                +O(exp(-4pi^2/a)),  a down to zero.            (15)

The product definition and eta modular transformation give this identity; see [DLMF 23.15](https://dlmf.nist.gov/23.15) and [DLMF 23.18.5](https://dlmf.nist.gov/23.18#E5). Compact zero modes modify the logarithmic term. Only the leading coefficient is used here.

### 4.1 Allowing arbitrary S and K leaves the original scalar problem

[Proved] There is an elementary converse to (14), sufficient for this logical point. Put b=-log(1-exp(-1))>0 and t=V' beta/L. For k<=floor(t/(2pi)), the summand is at least b. Hence

    B_S(2pi/t)>=S b floor(t/(2pi))
                   >=S b t/(2pi)-S b.                         (16)

Suppose the requested scalar target g<=K0+kappa M beta/L is known. Fix any positive numerical c and set V'=cM. Choose an integer S>=2pi kappa/(bc). Then

    g(L)<=K0+S b+B_S(2pi L/(V' beta)).                          (17)

Conversely (14) gives the scalar target from this product bound. Thus existence of a product majorant with arbitrary dimension-dependent S,K and V' of order M is equivalent, up to explicit changes of constants, to the original all-temperature scalar trace target. Replacing an unknown linear coefficient by an unknown species count does not by itself reduce the problem.

[Assumed or unresolved] No proof of that equivalent general statement, and no lattice chain violating it, has been obtained here. The rank counterexample in section 3 is not substituted for such a result.

## 5. A scalar bosonic majorant proved for XX

[Proved, inherited exact spectrum] Write

    f_L(u)=log(1+exp(-L asinh u)),
    ko_n=asinh((2n+1)pi/beta), n>=0,
    ke_n=asinh(2n pi/beta), n>=1.

C125 equation (6), after removing exp(L p_c), is

    Ztilde(L)=(P_o^+ + P_o^-)/2 + 2 exp(-L delta_beta) P_e^+,
    P_o^+ = product_(n>=0)(1+exp(-L ko_n))^4,
    P_o^- = product_(n>=0)(1-exp(-L ko_n))^4,
    P_e^+ = product_(n>=1)(1+exp(-L ke_n))^4,
    delta_beta>0.                                              (18)

The four fermion levels, the even-occupation projection in the first sector, and the two vacua in the second sector are all retained. Products converge for real L>1. On even L>=2 they reproduce the physical XX spin trace.

[Proved] f_L is decreasing and convex on u>=0. Indeed, g(t)=log(1+exp(-Lt)) has g'<0 and g''>0, while t=asinh u has t'>0 and t''<=0. Therefore f_L''=g''(t)(t')^2+g'(t)t''>=0. Let h=2pi/beta and I_L=integral_0^infinity f_L(u) du. Convex midpoint quadrature on each interval [nh,(n+1)h] gives

    sum_(n>=0) f_L((n+1/2)h) <= I_L/h.

Decreasing right-endpoint quadrature gives

    sum_(n>=1) f_L(nh) <= I_L/h.

Consequently P_o^+,P_e^+<=exp(2 beta I_L/pi). Since 0<=P_o^-<=P_o^+ and exp(-L delta_beta)<=1,

    log Ztilde(L)<=log 3+(2 beta/pi)I_L.                        (19)

[Proved] The integral identity also recorded in C125 section 7.1 is

    I_L=integral_0^infinity cosh t log(1+exp(-Lt)) dt
       =sum_(n>=1) (-1)^(n+1) L/(L^2 n^2-1)
       =[pi csc(pi/L)-L]/2,  L>1.                             (20)

For the series step, integration by parts gives L integral sinh t/(1+exp(Lt)) dt, and termwise integration gives the displayed absolutely convergent series. The partial fraction expansion of csc gives the last equality. Combining (19)-(20) proves

    log Ztilde(L)<=log 3+beta[csc(pi/L)-L/pi].                   (21)

[Proved] The positive-coefficient expansion of x/sin x, used in C125, shows that L[csc(pi/L)-L/pi] decreases for L>=2. Its value at L=2 is c0=2-4/pi. Thus

    log Ztilde(L)<=log 3+c0 beta/L,   c0=0.7267604553... .       (22)

The coefficient in (21) is the actual fixed-L zero-temperature limit of log Ztilde/beta. The additive log 3 is a uniform upper bound, not a claim of sharpness. The check L->infinity gives csc(pi/L)-L/pi=pi/(6L)+O(L^(-3)), reproducing A=1.

[Proved] A literal product of the form in the question follows too. Take S=2, V'=2e=4eM, and b as in (16). Since

    2e b/pi > 2-4/pi,

equations (16) and (22) give

    log Ztilde(L)<=K_XX+B_2(2pi L/(2e beta)),
    K_XX=log 3+2b=2.0159625794... .                            (23)

One elementary check on the strict constant inequality uses b>e^(-1)+e^(-2)/2, e<3 and pi<22/7. Then 2e b>2+1/e>7/3>16/7>2pi-4. No fitting is involved. For exchange J>0 replace beta by J beta and V' by 2eJ; M=J/2. Equations (21)-(23) hold for every beta>0 and even L>=2.

[Numerical] Both (21) and (23) pass the exact-ring XX grid L=2,4,6,16,64,256 and beta/L=0.01,0.1,1,4,32. The crossed product is independently evaluated with 20000 modes. At L=2,beta=64 its omitted log-trace tail is at most 0.0051877743, and the observed discrepancy is 0.0051875795. Other comparisons are recorded with their individual tail bounds. This is a controlled product truncation, not a claimed 1e-12 evaluation of that slowly convergent case.

[Proved extension] For n independent XX chains with exchanges J_a>0 on one tensor-product site, d=2^n and M=(1/2)sum_a J_a, as in C119. Their crossed traces multiply. Summing (21) gives

    log Ztilde(L)<=n log 3+2M beta[csc(pi/L)-L/pi].

Alternatively, sum (23) and enlarge each V_a=2eJ_a to V'=4eM. This gives S=2n and K=n K_XX. Microscopic fermionic Wick factorization, including its parity-sector bookkeeping, is a concrete sufficient extra lattice input in this family. It is not a property inferred for the interacting XXZ models.

## 6. The constant and zero-frequency states

### 6.1 A lattice necessity, K(d)>=log d

[Proved] Consider the diagonal ferromagnetic chain

    h=-J sum_(a=1)^d |aa><aa|, J>0,
    W_ab=exp(beta J delta_ab).

Its spatial transfer W is positive, with eigenvalues lambda_0=exp(beta J)+d-1 and lambda_1=exp(beta J)-1 of multiplicity d-1. It is a one-site translation-invariant nearest-neighbor chain in the allowed class. Its gauge-optimized M is J/2. The half-spread J/2 is attained by the displayed bond; no gauge can reduce it because the average bond energy on a uniform product ring is -J, whereas an alternating two-color even product ring has energy density zero. Telescoping gauges cancel in these expectation values.

[Proved] The normalized crossed trace is

    Ztilde(L)=1+(d-1)r_beta^L,
    r_beta=(exp(beta J)-1)/(exp(beta J)+d-1).

Choose even L(beta)=2 floor(exp(beta J/2)/2). Then beta/L->0 and

    log r_beta=-d exp(-beta J)+O(exp(-2 beta J)),
    L log r_beta->0,
    log Ztilde(L)->log d.                                     (24)

Every fixed finite-species oscillator product with V'=c(d)M tends to one on this sequence. Hence any scalar bound of the proposed form forces K(d)>=log d. A bound with K=0 fails on an actual positive-transfer lattice chain, regardless of S(d). These almost-degenerate classical sectors have zero time-circle momentum and are not controlled by a U(1) stiffness. Their inclusion is one concrete role for K. The exact ground energies are -JL and A=0, so this is not a counterexample to the target with an allowed constant.

### 6.2 A useful bound for the proposed Gaussian charge factor

[Assumed spectral comparison] Suppose a separately identified set of zero-frequency primaries is bounded by

    Q(beta/L)=product_a theta(rho_a beta/L),
    theta(u)=sum_(q in Z) exp(-q^2/(2u)), theta(0)=1,
    rho_a>=0, sum_a rho_a<=4M.                                (25)

This requires an energy lower bound for all the relevant charge sectors, with the stated charge normalization and with a controlled multiplicity. The Kohn budget alone is not that lower bound. In a generic Luttinger description stiffness is a response coefficient at small twist; identifying it with a uniform quadratic bound on all crossed charge energies, at all beta and through crossovers, is extra input. Neutral primaries and sectors invisible to those on-site charges are not included in (25).

[Proved under (25)] For any u>0, q^2>=q for integer q>=1 gives

    theta(u)<=1+2 sum_(q>=1) exp(-q/(2u))
             =1+2/(exp(1/(2u))-1)
             <=1+4u<=exp(4u).

Consequently

    log Q(beta/L)<=4(beta/L)sum_a rho_a<=16M beta/L.            (26)

The estimate does not require a bound on the number of charge directions. For countably many, apply it to finite products and use the increasing limit and summability of rho_a. It remains valid when some stiffnesses vanish. This is a deliberately loose bound that avoids an uncontrolled factor per charge. A factor for neutral ground multiplicity would still have to be included separately.

[Proved limitation] A slowly decaying positive-energy primary alone does not refute a scalar majorant with exp(K)>1, since the reference already has its vacuum contribution exp(K). Conversely, setting K=0 is not justified by saying that both traces tend to one at fixed beta and large L. The uniform two-parameter limit (24) prevents that choice. No estimate on all neutral primary multiplicities follows from the present insertion rank argument.

## 7. Where the many-insertion argument stops

### 7.1 Words in insertions and interactions

[Proved algebraic distinction] At n insertions one can form vectors labelled by n operator indices and n time arguments. A finite one-site algebra bounds the number of index choices at fixed times. It does not make the time arguments independent bosonic modes, does not make products of Fourier insertions commute, and does not imply an additive lower bound on the transfer energy of such products. Time ordering, their Gram matrix, and the crossed generator all enter. Even the free microscopic XX example gives a one-site density insertion that is already a two-fermion operator in the crossed description.

[Assumed sufficient input, not established for general chains] An effective Fock comparison would require an energy-filtered spanning or min-max statement. For example, let F_ref be a reference Hilbert space of occupation configurations and a specified zero-mode space. If there is an isometry J from the complete crossed Hilbert space into F_ref such that, as quadratic forms,

    Hc >= J^* H_ref J,                                        (27)

then the min-max principle puts every eigenvalue of Hc above the corresponding eigenvalue of H_ref. Thus Tr exp(-L Hc)<=Tr exp(-L H_ref). Equivalently one can directly prove an appropriate comparison of cumulative eigenvalue counts. A bounded number of additional reference copies supplies an additive K. Completeness, uniformity in Trotter number and beta, the zero-mode comparison, and the energy inequality are all required. They do not follow from (3) or (7). A non-isometric spanning map alone does not imply (27).

[Heuristic, precise point requiring proof] A polymer expansion would have to control the connected weights of arbitrarily many insertions so as to prove such an energy-filtered count. Reflection positivity makes a reflected quadratic form positive, but does not make every ordered operator word positive or control its connected weight by a product of two-point decays. Taking absolute values and bounding each word by the product of operator norms discards the required correlations and typically gives a cost of order beta M. A small-beta cluster convergence argument would leave beta/L large untreated. No all-beta convergence or comparison estimate of the needed form is proved here. The issue is the interacting multipoint weights, not just the number of matrix labels.

### 7.2 Gibbs comparison has a missing expectation in the exact state

[Proved] On a common finite-dimensional space, write H=H0+W, rho=exp(-LH)/Z and rho0=exp(-LH0)/Z0. Direct substitution into relative entropy gives

    D(rho0||rho)=L< W >_0+log Z-log Z0>=0,
    D(rho||rho0)=-L< W >_rho+log Z0-log Z>=0.

Therefore

    -L< W >_0 <= log Z-log Z0 <= -L< W >_rho.                  (28)

The usual computable reference-state variational expectation gives a lower bound on log Z. The upper bound needed here requires an expectation in the unknown crossed Gibbs state, or an independent operator/min-max comparison. Matching two-point functions to a Gaussian reference does not control W when the crossed generator has higher-order terms. Gaussian maximum-entropy statements additionally presuppose an appropriate canonical algebra; a finite one-site matrix algebra does not provide it. Vacuum-energy normalization must also be matched. Thus this proposed route has not supplied an upper trace bound.

### 7.3 Band limitation does not count histories

[Proved consequence of (8)] If E<1, all single-insertion components in the spectral interval [0,E] have |omega|<=VE. High frequencies have alpha>=1 and are excluded in the detected subspace. Without reconstruction of the full Hilbert space from controlled products, this is not a statement about every crossed state.

[Proved obstruction] In a many-insertion history only the sum of the Fourier frequencies is the external momentum. Large positive and negative internal frequencies can cancel. Equation (11) is an explicit instance: fixed external omega_k, arbitrarily many internal pairs, and many distinct low energies as beta grows. A time-bandwidth count of the one-variable insertion map therefore does not bound the cyclic or many-insertion space. Applying such a count to n-variable histories also requires a bound on n and its energy cost. Those are the missing interacting data.

[Proved] If one only had log N(E)<=K+c beta E, integration by parts gives

    Tr exp(-L Hc)=L integral_0^infinity exp(-LE)N(E) dE
                 <=exp(K) L/(L-c beta),                        (29)

and this argument applies only when L>c beta. It says nothing in the required opposite regime. In contrast, an additive Fock reference with frequencies k and S species obeys, by (14) and an exponential cutoff,

    N_ref(E)<=exp[t E+S pi V' beta/(12t)] for every t>0,
    log N_ref(E)<=sqrt(S pi V' beta E/3).                       (30)

The square root is obtained by optimizing t. Producing this additive energy structure, or a comparably strong count, is the substantive missing input.

## 8. A sufficient comparison that tolerates lattice ultraviolet dispersion

[Proved conditional closure] The reference need not have a globally linear dispersion. Suppose (27), including completeness and multiplicities, can instead be proved with at most S oscillator species per positive Matsubara index, counting both signs, and one-particle energies alpha(2pi k/beta) from (8). Include a factor exp(K0) and, if desired, the charge factor (25). The resulting oscillator trace obeys the desired beta/L bound for L>=2.

[Proved] To see this, f(omega)=-log(1-exp(-L alpha(omega))) is decreasing, so

    log Z_ref,osc <= S beta/(2pi) integral_0^infinity f(omega) d omega.

For 0<=omega<=V, use the linear branch and extend that integral to infinity, obtaining V pi^2/(6L). For omega>=V=4eM, exp(-L alpha)=(4M/omega)^L<=exp(-L). The elementary bound -log(1-z)<=z/(1-z) gives

    integral_V^infinity f(omega) d omega
       <= V exp(-L)/[(L-1)(1-exp(-L))]
       <= (V/L) 2/(exp(2)-1),  L>=2.

Therefore

    log Z_ref,osc <= S V beta/(2pi L)
                        [pi^2/6+2/(exp(2)-1)].                 (31)

Adding (26) yields

    log Ztilde(L)<=K0+
       { (2e S/pi)[pi^2/6+2/(exp(2)-1)]+16 } M beta/L,         (32)

with the last 16 omitted when no charge factor is needed. This proves that a logarithmic ultraviolet dispersion is compatible with the target. A false global linear cone is unnecessary. The only unresolved step in this conditional theorem is the substantive all-state comparison, together with any assumed zero-sector estimate, not convergence of its reference product.

[Assumed or unresolved] There is no uniquely ordered notion of a minimal extra lattice axiom. Relative to this route, (27) or a comparable energy-filtered many-insertion count states the missing input. It is stronger than the scalar target and has not been derived from fixed d and nearest-neighbor locality. For free XX, Wick factorization, the explicit spin-structure sum and convex quadrature already supply a weaker sufficient comparison, section 5. Claiming Wick domination or additive insertion energies for an interacting chain would be an additional assumption, not a conclusion of reflection positivity.

## 9. Finite-Trotter checks and their scope

[Numerical construction] The scripts implement C125's reflected XXZ bond h'=-Sx Sx+Sy Sy-Delta Sz Sz, g=exp(-beta h'/m), and the two positive column matrices Ke,Ko on histories of length 2m. A shift by two history entries represents time translation beta/m and commutes with Ke and Ko. Orbit Fourier bases block-diagonalize it with eigenvalues exp(2pi i q/m). Within each block the script diagonalizes U=Ke^(1/2)KoKe^(1/2). For an eigenvalue u, the energy relative to the largest eigenvalue u0 is E=-(1/2)log(u/u0). The factor 1/2 is needed because U propagates two spatial sites.

[Numerical] The full local grid is Delta=0,0.5,1,2, beta=0.25,1,2,8 and m=2,3,4, totaling 48 transfer matrices of dimensions at most 256. Each is compared at L=2,4,8,16,64. Every discrete momentum sector and every retained crossed energy is saved. Momentum is known only modulo m. The aliased comparison uses the reference product with its oscillator momenta reduced modulo m as well. Positive cyclic convolution computes its sector traces without subtracting nearly equal Fourier sums.

[Numerical] All scalar samples pass

    log Ztilde_m(L)<=log 3+B_6(2pi L/(4eM beta)),               (33)

with minimum margin 0.0082541890. Equation (33) is a tested candidate, not an interacting theorem. Sector versions fail, with a sampled log-margin -259.5066468 for Heisenberg at beta=0.25,m=4,L=64,q=2. Exact XX already proves the continuum momentum-preserving obstruction in section 3.1, so this finite-step failure is not used to infer a new continuum theorem.

[Numerical] The two-column contraction agrees with the independently assembled physical four-site Trotter trace to 4.03e-12 in log trace. At beta=2,m=4, errors relative to the exact physical four-site trace are 0.00587127 for XX, 0.01782895 for Delta=0.5, 0.03950913 for Heisenberg, and 0.09952505 for Delta=2. Errors are larger at beta=8. At beta=8,m=2 the Heisenberg sector q=1 violates the continuum low-frequency cone by 0.0724667252. This is a coarse-time-discretization warning, not a contradiction of (5)-(7). The scripts retain it as an expected failure.

[Numerical scope] Eigenvalues below 1e-13 times the largest U eigenvalue are omitted from gap lists. Their total normalized trace contribution is bounded by dimension times (1e-13)^(L/2), and the actual omitted count is used. No sampled scalar pass depends on that truncation. Positivity signs and eigenvalues are floating-point checks, not certified enclosures. The exact XX products, integral estimates, Gaussian sum bound, elementary boson bounds, square-root Fock count and the directions in (28) are separately tested in audit_checks.py.

[Work record] Only small transfer matrices were run locally. run_typhon.slurm and checks.py --cluster prepare m=5,6 checks and refuse that mode without a SLURM job identifier. They were not submitted or run. No network shell command, outside write, or independent agent review was used. A local STATE.md supplies the handoff because project-level writes are prohibited for this run.

## 10. Status of the requested decision

[Proved] The single-insertion rank statement is true, but it is not the required spectral multiplicity statement. The actual XX chain disproves that inference, and disproves global momentum-preserving linear bosonic majorization. The diagonal ferromagnetic chain forces K(d)>=log d. XX nevertheless satisfies an explicit scalar two-species majorant at every allowed length and temperature.

[Assumed or unresolved] The arbitrary interacting scalar majorant with some explicit K(d), S(d), and V'=O(M) has neither been proved nor refuted here. Tasks 2 and 4 are therefore only partially resolved: the failed stronger interpretations have lattice examples, while the scalar target itself remains open. No abstract spectral countermodel is offered as a lattice counterexample. The precise additional many-insertion input needed by this route is stated in sections 7 and 8, but it is not a proved property of the requested interacting chains.
