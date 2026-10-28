# C59. Thermal comparison for weak XXZ

## 1. Status and conventions

[proved partial results; interacting theorem unproved] This note does not establish estimate P for XXZ at a fixed nonzero anisotropy. It identifies the step that prevents the proposed extension of the Ising proof, proves several ingredients of that extension, and states the remaining lemma with its scalar normalization requirement. There is no contradiction to the expected thermal formula.

[read] C57's four-history bootstrap and first-order interaction calculation are used as established local inputs. C58 has closed the physical thermodynamic transfer, conditional on the matched-action infrared estimates. Neither result proves the remaining thermal comparison. The literature search and its access limitations are recorded in SOURCE_AUDIT.md.

Use C57's Hamiltonian

    H = -t sum_x (c_x^dagger c_(x+1)+h.c.)
        +V sum_x (n_x-1/2)(n_(x+1)-1/2)-mu sum_x n_x,
    t>0, Delta=V/(2t).                                      (1)

Fix a compact interval of Fermi points inside (0,pi), excluding pi/2. Smallness constants may depend on that interval and t. Take the spatial thermodynamic limit before T tends to zero, at fixed physical parameters. Counterterms are selected at zero temperature and then kept fixed. Reflection gives equal right and left speeds for (1).

[proved boundary-condition check] An end bond of uniformly bounded norm changes the free energy per site by at most its norm divided by L. Indeed the ordered eigenvalues of H and H+B differ by at most ||B||, so exp(-beta||B||)Z_H<=Z_(H+B)<=exp(beta||B||)Z_H. Open spin XXZ is exactly Jordan-Wigner equivalent to open fermions. Compare each periodic chain with its open chain using this bound. Their thermodynamic pressures agree at fixed T without identifying a periodic spin trace with one particular fermion parity sector. This elementary comparison does not control a finite-size Casimir coefficient.

Write q_h=2^(-h), h>=0, after fixing an ultraviolet energy unit E_*. Thus T and q_h below are dimensionless. A physical free-energy density is restored by multiplication by E_*; an error O(T^(2+theta)) means O(E_* (T_phys/E_*)^(2+theta)). Let T<=q_H<2T. No smallness condition here is uniform as the filling approaches a band edge or half filling.

## 2. What changes in the Ising proof

[read] Giuliani-Mastropietro [1211.3030, Section 3.4](https://arxiv.org/pdf/1211.3030#page=32) first obtain (3.84), an irrelevant gain in each interaction scalar. Equation (3.89) splits the finite-size correction into a fixed-history geometric difference and a history difference. Equations (3.69), (3.90)-(3.91) control the latter and its sum; the paragraph after (3.91) treats the former. The relevant power is 3-epsilon, greater than two. Their Section 1 expressly leaves the interacting c=1 extension open at that time.

[proved distinction] The corresponding XXZ objects must be

    r_h(T)=E_(h,T)(X_h^T)-E_(h,T)(Y_h^T),
    delta_h(T)=r_h(T)-r_h(0),                                (2)

where Y is the complete interacting reference, not its quadratic fermionic part. The terms E include all covariance-change scalar factors. The sought analogs are

    |r_h(0)| <= C q_h^(2+theta),
    |delta_h(T)| <= C q_h^(2+theta)(T/q_h)^p,
    p>2+theta.                                              (3)

An individual interacting reference shell need only be O(q_h^2). Its marginal quartic does not vanish in the infrared. Therefore replacing the Ising interaction scalar e_h in (3.84) by the XXZ interaction scalar is not a valid inference. One must first subtract the reference, then carry that subtraction through the geometric and history estimates.

[proved algebra, using C57] Put

    y=Y_h^0, z=X_h^0-Y_h^0, a=Y_h^T-Y_h^0,
    d=(X_h^T-Y_h^T)-(X_h^0-Y_h^0).

For any scalar map E_T the difference in (2) is the sum

    E_T(y+a+z+d)-E_T(y+a+z),
    E_T(y+a+z)-E_T(y+a)-E_T(y+z)+E_T(y),
    (E_T-E_0)(y+z)-(E_T-E_0)(y).                          (4)

The second term is the integral of D^2 E_T[a,z] over the parameter rectangle. It is present even when the limiting reference beta function vanishes. C57 proves the resulting bootstrap; that proof is not repeated here.

## 3. The zero-temperature data that can actually be extracted

### 3.1 Localized couplings and field normalization

[read] Benfatto-Mastropietro [0907.2837, (81)-(87), (92)-(95)](https://arxiv.org/pdf/0907.2837#page=17) control the constant-history beta function, convergence of the marginal coupling, anomalous field-normalization ratios and matching of the lattice/reference field amplitudes. In their infrared index j<=0,

    |lambda_j-lambda_-infinity| <= C_theta lambda^2 gamma^(theta j).

These are useful zero-temperature inputs. Their Theorem 1.1 is a response/exponent theorem, not a thermal-pressure theorem.

[proved elementary normalization lemma] Suppose two positive normalization sequences obey

    Z_(h+1)^X/Z_h^X=exp(eta+b_h^X),
    Z_(h+1)^Y/Z_h^Y=exp(eta+b_h^Y),
    |b_h^X-b_h^Y|<=B q_h^theta, theta>0.                    (5)

Then there is a constant c>0 such that

    |log[Z_h^X/(c Z_h^Y)]|
         <= B q_h^theta/(1-2^(-theta)).                    (6)

Indeed the logarithm of Z_h^X/Z_h^Y equals its value at h=0 plus the sum of b_k^X-b_k^Y for k<h. The series converges absolutely. Choose log c to be its infinite sum, including the initial value. The remaining difference is minus the tail beginning at h, bounded by the geometric series in (6). Exponentiation gives Z_h^X/(c Z_h^Y)=1+O(q_h^theta). The individual Z_h can diverge. This argument uses equality of the anomalous exponent eta; it does not prove that equality or a finite-temperature version of (5).

The BM matching supplies the relevant zero-temperature normalization mechanism for the instantaneous model. It does not by itself identify the complete Banach-space defect z in (4).

### 3.2 Curvature is a genuine marked irrelevant factor

[proved shell estimate] At a fixed Fermi point p, expand the free lattice dispersion in the right or left patch,

    e(omega p+k)=omega v_F k+r_omega(k),
    v_F=2t sin p, |r_omega(k)|<=C k^2,
    omega in {+1,-1}.                                     (7)

Here the constant part has been removed by the chemical potential. Choose a smooth shell chi_h supported where |k0|+v_F|k| is comparable to q_h, with bounded rescaled derivatives. For sufficiently small fixed initial patch, the interpolating denominators

    D_s=-i k0+omega v_F k+s r_omega(k),  0<=s<=1

have |D_s|>=c q_h on this support. To see this, |D_0| is comparable to q_h and |r|<=C q_h^2; shrink the initial patch until the latter is at most half the former. The resolvent identity gives

    chi_h/D_1-chi_h/D_0
       =-integral_0^1 ds chi_h r_omega(k)/D_s^2.           (8)

After k0=q_h K0, k=q_h K, the symbol on the right and any fixed number of its rescaled derivatives are bounded uniformly. Its Fourier transform therefore obeys, for every fixed integer m,

    |g_h^lat(x0,x)-g_h^lin(x0,x)|
        <= C_m q_h^2 (1+q_h|x0|+q_h|x|)^(-m),
    integral dx0 sum_x (1+q_h|x0|+q_h|x|)^m
                    |g_h^lat-g_h^lin| <= C_m.             (9)

For the second inequality use the first with power m+4. Fourier integration by parts proves the first; the patch cutoff is smooth at its boundary, so there are no boundary terms. In contrast, the corresponding bounds for g_h^lin have prefactors q_h and q_h^(-1). Thus a curvature replacement gains q_h. The same proof works for a uniformly smooth small perturbation of the local coefficients, provided the shell denominator stays uniformly invertible.

[proved limitation] Equation (9) is a covariance estimate, not a bound on the sum of renormalized trees containing that replacement. The Gaussian field factors must be normalized before using it. Propagating this gain through arbitrary nested localizations is an additional step.

### 3.3 Why the full reference cannot be truncated to its local quartic

[proved coefficient-level example] Consider a continuum linear reference with local interaction

    g integral dx rho_+(x) rho_-(x),
    rho_omega=psi_omega^+ psi_omega^-.

At second order, contracting one plus-chirality pair between vertices x,y produces a nonlocal six-field contribution proportional to

    g^2 g_(+,h)(x-y)
          psi_+^+(x) psi_+^-(y) rho_-(x) rho_-(y),          (10)

with its antisymmetrized partner. This is a nonzero Grassmann kernel for distinct points. In the linear continuum shell,

    g_(+,h)(x-y)=q_h G_+(q_h(x-y)).

The pinned absolute integral of this coefficient is proportional to g^2/q_h. Its six-field dimensionless norm multiplies that integral by q_h^(6/2-2)=q_h. It is therefore of order g^2, with no vanishing q_h factor. Antisymmetrization kills the fully local six-field monomial, but not the kernel at separations of order q_h^(-1). Equivalently, the rescaled polynomial in rescaled fields is the same nonzero polynomial at every shell.

This is a check of a specific perturbative kernel, not a lower bound for an arbitrary representation of the complete interacting effective action. It suffices to show why engineering irrelevance alone does not justify discarding all such reference kernels in the full effective-potential norm. They cancel against their lattice counterparts only in a matched comparison. The bound on lambda_h-lambda_* is not a bound on that entire difference.

## 4. Gaussian normalization factors admit the required thermal bound

[proved analytic lemma] This section treats smooth shell logarithms at prescribed running data. It does not assume a bound on the actual thermal histories.

Let a range over a convex parameter domain, and let Phi_h(k0,k;a) be a scalar symbol. It can, for example, be a matrix log determinant. Suppose it and its first two parameter derivatives are supported in a fixed rescaled shell or patch, and, for some integer m>=3, obey

    integral dk0 dk/(2pi)^2
            |partial_(k0)^j D_a^r Phi_h(k0,k;a)|
                  <= C q_h^(2-j),
    0<=j<=m, 0<=r<=2.                                    (11)

Use operator norms for the parameter multilinear forms. Smoothness sufficient for Fourier inversion and integration by parts is understood; compact smooth symbols suffice. Spatial momentum can be on the Brillouin zone. Define

    G_h(T,a)=T sum_(n in Z) integral dk/(2pi)
                   Phi_h((2n+1)pi T,k;a),
    G_h(0,a)=integral dk0 dk/(2pi)^2 Phi_h(k0,k;a).         (12)

Then, uniformly for 0<T<=q_h,

    ||D_a^r[G_h(T)-G_h(0)]||
                    <= C_m q_h^2(T/q_h)^m,  r<=2,
    ||D_a^r G_h(T)||+||D_a^r G_h(0)||<=C_m q_h^2.        (13)

Proof. With the time Fourier transform defined using dk0/(2pi), the antiperiodic Poisson formula gives

    G_h(T,a)-G_h(0,a)
        =sum_(ell!=0) (-1)^ell integral dk/(2pi)
                                      Phi_hat_h(ell/T,k;a).

Integrate by parts m times in k0. The integrated absolute value of the ell-th term, including r parameter derivatives, is at most

    C (T/|ell|)^m q_h^(2-m).

The sum of |ell|^(-m) is 2 zeta(m), finite. This proves the first estimate. The second follows from (11) at j=0 and T/q_h<=1. All differentiations and sums are justified by the displayed absolute majorants. QED.

[proved sufficient symbol condition] Suppose A_h(k;a) is a finite matrix, has a fixed compact rescaled support, has uniformly bounded rescaled momentum derivatives and first two parameter derivatives, and ||A_h||<=r<1. Then

    Phi_h=log det(I+A_h)
         =sum_(n>=1) (-1)^(n+1) tr(A_h^n)/n             (14)

has (11). Each fixed derivative order costs a polynomial in n, and the remaining powers of A_h give a convergent geometric sum. Rescaling the two momenta supplies q_h^2 and each frequency derivative supplies q_h^(-1). This defines the branch of the logarithm continuously from A_h=0, also for complex symbols. Small covariance-update ratios on a smooth shell are an example.

[proved four-history consequence] Apply (4) with E=G. If the state differences satisfy

    ||z||<=C q_h^theta,
    ||a||<=C (T/q_h)^p,
    ||d||<=C q_h^theta(T/q_h)^p,
    2<p<=m,

then (13) bounds the three terms by, respectively,

    C q_h^2 ||d||,
    C q_h^2 ||a|| ||z||,
    C q_h^2 (T/q_h)^m ||z||.

Their sum is at most C q_h^(2+theta)(T/q_h)^p. Thus Gaussian normalization factors of the form (12)-(14) cause no separate loss of thermal power once their actual histories are controlled.

[assumed application boundary] A physical RG implementation must derive its covariance-change factors as (12), or an equivalent controlled expression, from the finite-regulator trace identity. One cannot apply this lemma to an unregularized absolute product of Z_h factors or choose its scalar normalization freely. The lemma establishes the bound for the specified smooth symbols, not that every omitted factor in an unspecified scheme has this form.

## 5. An explicit history norm that tolerates a running marginal coupling

[proved conditional sequence lemma] The following is a concrete version of C57's weighted-history criterion. It does not require the reference coupling to tend to zero. Suppose the mixed distance satisfies

    d_(h+1) <= ell d_h
                 +c epsilon sum_(k=0)^h 2^(-sigma(h-k)) d_k
                 +B q_h^theta (T/q_h)^p,                 (15)

where d_h>=0, sigma>0, ell>=0, epsilon is the small bound on running interactions, and 0<theta<p-2. Define

    w_h=q_h^theta (T/q_h)^p,
    r=2^(-(p-theta))
             [ell+c epsilon/(1-2^(-(sigma+p-theta)))].     (16)

If r<1 and d_0<=D_0 w_0, then d_h<=D w_h, with

    D=max(D_0, B 2^(-(p-theta))/(1-r)).                    (17)

Proof. Since w_k/w_(h+1)=2^(-(p-theta))2^(-(p-theta)(h-k)), the history sum divided by w_(h+1) is at most

    D 2^(-(p-theta)) sum_(j>=0)2^(-(sigma+p-theta)j).

The geometric series gives (16). Equation (15) therefore propagates d_(h+1)/w_(h+1)<=rD+B 2^(-(p-theta))<=D. The initial step is included. This proves (17).

The corresponding pure reference shift uses theta=0. If (16) holds, its analogous contraction condition also holds: replacing p-theta by the larger p decreases both terms in the bracketed normalized expression. A forcing B(T/q_h)^p and initial bound O(T^p) then give ||Y_h^T-Y_h^0||<=C(T/q_h)^p.

[proved feasibility check] For p=3, theta=1/4, sigma=1 and ell=2, the epsilon=0 value is 2/2^(11/4)<1. Hence an O(epsilon) history perturbation is allowed. Relevant one-step amplification by two is compatible with the bootstrap. The required thermal weight grows toward the infrared faster than that amplification. This calculation does not verify (15) for XXZ.

[proved logical distinction] The zero-temperature vanishing beta function is used to match the limiting marginal parameters. It does not assert that DR of the reference map is zero, nor that its vacuum functional is free. A surviving marginal coupling is compatible with (15); proving the joint thermal and defect forcing in its last term is the substantive missing step.

## 6. The reference pressure and the scalar identification

[proved reference extension] C51 already evaluates the compact Gaussian reference with constant speed. The same leading coefficient holds for a Hamiltonian reference whose two oscillator branches have energy

    epsilon(k)=k u(k), k>0,
    u(k)>=u_min>0, u(k)->v>0 as k->0,

and whose oscillator thermodynamic free energy is

    f_ref(T)-e_ref=(T/pi) integral_0^infinity
                               log(1-exp(-k u(k)/T)) dk.  (18)

Set k=Tx. The integrand divided by T^2 tends pointwise to pi^(-1) log(1-exp(-vx)). Its absolute value is bounded by -pi^(-1) log(1-exp(-u_min x)), which is integrable at both endpoints. Dominated convergence gives

    f_ref(T)-e_ref=-pi T^2/(6v)+o(T^2),                   (19)

using integral_0^infinity log(1-exp(-vx))dx=-pi^2/(6v). The last identity follows by expanding -log(1-z)=sum_(n>=1) z^n/n and integrating the positive series termwise. If the branches are cut off at a fixed momentum, extend the rescaled integrand by zero; the same argument applies. The compact charge sectors have subextensive logarithmic partition contributions by C51 and hence do not alter (18) in the thermodynamic-first limit.

[proved normalization distinction] Equation (18) is a Hamiltonian trace statement. A comparison of normalized fermionic correlators or a Ward identity for them does not identify it with an ultraviolet-regularized Grassmann vacuum functional. For example, multiplying every source partition function by exp[-beta L b(T)] leaves Z(J)/Z(0) unchanged while shifting the free energy by b(T). Taking b(T)=aT^2 demonstrates the missing information. This is an ambiguity of source-normalized data, not a claim that an arbitrary b(T) is realizable by the fixed XXZ Hamiltonian.

[assumed, necessary bridge] The reference used in the lattice subtraction must either be the trace-normalized solvable Hamiltonian reference from the outset, or its full scalar must be compared to that reference with o(T^2) error. The identity includes the ultraviolet regulator and every determinant. The zero-temperature density Ward identities do not supply this bridge. No proof of it for a chosen common XXZ/reference RG scheme is supplied here.

## 7. Exact remaining lemma

[assumed, unproved lemma TP] For (1), at a fixed off-half-filled Fermi point and sufficiently small |Delta|, construct the lattice and the matched Hamiltonian Luttinger reference with the same infrared scheme and their physical trace normalizations. The reference has the dressed speed v and limiting forward coupling. Match its finite field amplitude as in (6). Use the four-node thermal localization of the collaborator's THERMAL_LOCALIZATION.md, or another specified scheme with power p>2. Require the following estimates uniformly in the number of integrated shells and thermodynamic regulators.

(TP1) In one common normalized effective-kernel/history norm, including all nonconstant kernels and the local covariance data,

    ||X_h^0-Y_h^0|| <= A q_h^theta, 0<theta<p-2.            (20)

One concrete norm to target, in a specified common slow-field representation, is the sum of the absolute normalized local covariance parameters and

    sum_(l even>=2) sum_(derivative labels alpha)
        B^l q_h^(l/2+|alpha|-2)
        integral |W_(h,l,alpha)(0,X)|
                         (1+q_h diam(X))^5 dX,             (20a)

where B>1 is fixed, the first field coordinate is pinned, and local affine quadratic and constant quartic terms are extracted into their coefficient coordinates. Spatial integrals are replaced by sums when the common representation is on the lattice. Derivative labels retain the Taylor-remainder field derivatives and are included in the displayed dimension. An equivalent fixed polynomial terminal-tree norm dominating (20a) is allowed, with its equivalence constants stated independently of h. The reference has the same labels and normalization. Equation (20) must control this complete difference, not only its localized coefficients.

The construction must supply a fixed stronger field and moment reserve for the marked trees. For the individual integration/localization steps the required additional moment orders are given in WINDING_ESTIMATE Sections 8-9 and THERMAL_LOCALIZATION. At fixed output moment order, those reserves must not grow with h. No claim that an arbitrary kernel controlled only in (20a) has uniformly bounded derivatives of every RG map is being made.

(TP2) The actual four histories stay in the corresponding common domains. The pure reference shift and mixed difference obey (15) and its theta=0 counterpart, with r<1, or the one-step derivative hypotheses (5) of C57. The ultraviolet pure shift is O(T^p), the mixed shift is O(A T^p), with the same fixed bare counterterms at all temperatures. The forcing B in the mixed recurrence is independent of h and T.

(TP3) The interaction part of each vacuum increment, on those domains, obeys

    ||D E_(h,T)^int||+||D^2 E_(h,T)^int|| <= C q_h^2,
    ||D(E_(h,T)^int-E_(h,0)^int)||
                                  <= C q_h^2(T/q_h)^p.   (21)

The derivative and Hessian norms must be compatible with the history/state norm of TP1-TP2. Gaussian updates are retained as smooth symbols satisfying Section 4, or are bounded directly by (21). The exact scalar expansion also has an ultraviolet four-term difference O(T^p), and either a terminal Lipschitz constant C q_H^2 or an exhausted shell partition with no residual covariance. The reference scalar is the physical quantity (18), or has an independently controlled o(T^2) difference from it.

[proved conditional theorem] TP implies

    f_XXZ(T)=e0-pi T^2/(6v)+o(T^2).                       (22)

Proof. TP1-TP2 and Section 5, or C57's one-step bootstrap, give ||d_h||<=C q_h^theta(T/q_h)^p and ||a_h||<=C(T/q_h)^p. Equations (4), (21) and Section 4 give (3). The terminal difference is O(q_H^(2+theta)) by its Lipschitz bound and ||X_H^T-Y_H^T||<=||z_H||+||d_H||. The missing zero-temperature tail is bounded by the sum of C q_h^(2+theta) for h>=H, hence O(T^(2+theta)). For active shells,

    sum_(h<H) q_h^(2+theta)(T/q_h)^p
        =T^p sum_(h<H) 2^((p-2-theta)h)
        <=C T^(2+theta).                                (23)

Together with the ultraviolet difference this proves

    [f_XXZ(T)-e0]-[f_ref(T)-e_ref]=O(T^(2+theta)),

when the reference scalar is exact. An additional o(T^2) reference-normalization error leaves the final conclusion unchanged. Insert (19). QED.

[assumed, precise unproved tree operation] To establish TP from the known XXZ expansion, the GN trees must be compared with both a lattice/reference defect and a thermal replacement retained at once. The defect can be a curvature line, a difference of running local couplings, or a difference of complete effective kernels. The thermal replacement can be a periodized line, a localization difference, or a history variation. The absolute, renormalized tree sum with these two operations must preserve the factor q_h^theta(T/q_h)^p, together with one history derivative or two state derivatives as required above. Scalar roots are included. This is the uniform joint marked estimate still missing here.

The collaborator's WINDING_ESTIMATE.md proves the fixed-data one-shell version, including scalar roots and Cauchy derivatives. THERMAL_LOCALIZATION.md provides p=3 with the stated moments. MOMENT_RESERVE.md supplies those moments initially from the matched action. Its Section 6 also proves that, once the ordinary marked-tree distance routing is controlled, one thermal mark costs only a fixed additional moment and a polynomial in the number of fields/endpoints. That conditional allocation result is already available and is not an additional missing combinatorial lemma. What is not verified here is its application to the full XXZ/reference comparison, the actual four histories and the scalar normalization, with the norms and bounds in TP. The zero-temperature BM estimates and (9) provide defect ingredients, not that complete application. Equations (20)-(21) are not declared consequences of the phrase "vanishing beta function".

## 8. The matched retarded action

[read, imported] C48 and the massive-matching work supply exponentially decaying all-degree endpoints and an exponentially small thermal initialization error. C55 proves the linear one-endpoint all-scale estimate. C55's nonlinear Lemma N and C56's multichannel comparison remain separate inputs. C58 transfers a uniform logical thermal formula to the physical ensemble with activated accuracy.

[proved conditional scope] If TP is proved only for the nearest-neighbor instantaneous Hamiltonian (1), it does not apply to an arbitrary supplied retarded action. Its extension must allow all even endpoint degrees and the fixed field/moment norm from massive matching. Static shifts of the relevant term, kinetic residue and marginal coupling must first be matched, including terms linear in the extra endpoint size. The dressed speed must include the frequency residue. An exponentially small initial thermal error fits TP2 but does not prove the nonlinear tree estimates. These distinctions persist even if the announced manuscript eventually proves (22) for instantaneous chains.

If reflection is retained, the single-channel target is (22) with the rematched v. Without reflection it is

    f(T)=e0-(pi/12)(1/v_R+1/v_L)T^2+o(T^2).

For the original many-channel matched action, both the source-uniform all-scale input (a) and its scalar input (b) are still required. C58's ensemble transfer does not need to be redone.

## 9. Checks and boundary of the result

[proved checks] The curvature gain has the expected one extra power of energy. The Gaussian shell bound has dimension q_h^2 and vanishes identically at T=0. Matching field amplitudes removes only a finite ratio, not the anomalous exponent. The history condition allows a relevant multiplier two and a nonzero small marginal coupling. At V=0 the known free coefficient is recovered. C57's already checked first Taylor coefficient is consistent with (22); it is not recomputed or upgraded to an all-order estimate here.

[numerical] checks.py tests the Poisson identity for a smooth log-determinant symbol, including its first two parameter derivatives, and the explicit history bound (15)-(17). These are controls of the partial lemmas, not computations of the XXZ interacting flow or of its thermal central charge.

[proved work boundary] The all-depth joint marked estimate and the scalar identification in TP were not established. Consequently task (ii) remains incomplete, and (22) is a conditional theorem. The exact remaining lemma is Section 7. No claim of an unconditional XXZ or matched-retarded thermal theorem is made.
