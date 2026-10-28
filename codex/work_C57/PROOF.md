# C57. Thermal vacuum comparison and the interacting XXZ check

## 1. Status and conventions

[proved partial results; requested theorem remains unproved] This note proves a conditional four-history comparison theorem and computes the first Taylor coefficient in the XXZ interaction of the low-temperature free energy. It does not prove estimate P for the full interacting XXZ chain at a fixed nonzero interaction. The retarded and multichannel applications therefore also remain conditional.

[read] The source audit is in [SOURCE_AUDIT.md](SOURCE_AUDIT.md). It includes a directly relevant 2024 announcement by Renzi, Giuliani, Mastropietro and Santini. No full theorem statement from that work was located. There is no basis here for claiming a new unconditional thermal c=1 theorem, or that such a theorem has never been announced.

Use lattice spacing, hbar and k_B equal to one. The spinless Hamiltonian is

    H(V,mu) = -t sum_x (c_x^dagger c_(x+1)+h.c.)
              +V sum_x (n_x-1/2)(n_(x+1)-1/2)
              -mu sum_x n_x,                 t>0.                 (1)

[proved identification] Jordan-Wigner and, if desired, the staggered spin rotation give the XXZ convention J=2t, Delta=V/(2t). A field-dependent scalar caused by replacing n by S^z+1/2 changes the ground energy, not the thermal difference at fixed field. We work with (1) to keep that convention unambiguous.

At V=0 put p=arccos(-mu/(2t)), v=2t sin p. Fix p in the interior of (0,pi), away from pi/2 for the requested off-half-filled application. Every smallness constant may depend on this fixed choice. Thermodynamic length tends to infinity first, then T tends to zero at fixed microscopic parameters. Physical counterterms are never chosen afresh as functions of T.

[read, inherited results] C51 Section 7 proves that P implies the desired coefficient. Its Sections 2–6 evaluate and diagonalize the Gaussian reference. The collaborator's [thermal criterion](../../contrib/thermal_ensemble_20261028/thermal/THERMAL_CRITERION.md), [winding theorem](../../contrib/thermal_ensemble_20261028/thermal/WINDING_ESTIMATE.md) and [localization lemma](../../contrib/thermal_ensemble_20261028/thermal/THERMAL_LOCALIZATION.md) already treat summation, fixed-history shell integration and Matsubara interpolation. Their proofs were read. None is being presented as a new C57 result.

## 2. What the finite-temperature RG has to compare

[assumed representation, not a theorem for (1)] Choose a smooth shell decomposition with q_h=E_* 2^(-h), h>=0. Let H be the integer with T<=q_H<2T. With complete scalar bookkeeping the expansion must have the form

    f_X(T)=a_X(T)+sum_(h=0)^(H-1) E_(h,T)(X_h^T)
                         +Q_(H,T)(X_H^T),
    e_X=a_X(0)+sum_(h>=0) E_(h,0)(X_h^0).                       (2)

X is the lattice theory; Y will denote its dressed Luttinger reference. The scalar E includes the Gaussian normalization produced when changing covariance or field normalization. Q is the entire terminal integration. Equation (2) is a partition-function identity, not just a convention for a normalized source functional.

The effective data X_h must include enough information to perform the next integration. A possible choice includes the complete nonconstant effective potential, normalized quadratic data and shell covariance. If an expansion instead uses histories of marginal couplings, those histories and their norm are part of the data. Absolute Z_h need not converge. One must work with matched normalizations or ratios and retain their scalar Jacobians.

[proved algebraic distinction] The shell error relevant to P is

    delta_h(T) = [E_(h,T)(X_h^T)-E_(h,T)(Y_h^T)]
                  -[E_(h,0)(X_h^0)-E_(h,0)(Y_h^0)].             (3)

It involves four states. An estimate of each individual thermal change with constant proportional to the small interaction gives, at best, C q_h^2 epsilon (T/q_h)^p. It does not yield a factor tending to zero with h. Similarly, separate bounds by q_h^2 eta_h and q_h^2(T/q_h)^p only give their minimum, not their product.

[proved counterexample to that last inference] Set delta_h=q_h^2 min(eta_h,(T/q_h)^p) with eta_h=2^(-theta h). Both separate bounds hold. At T/q_h=eta_h^(1/p), the ratio to q_h^2 eta_h(T/q_h)^p is 1/eta_h, unbounded in h. This is a counterexample to combining two inequalities, not a physical XXZ counterexample.

[assumed reference identification] The fermionic reference in (2) must have its full pressure identified with the Gaussian current theory with the limiting forward coupling. It is not the free fermion model with the bare velocity. In the reflection-symmetric case its thermal contribution is -pi T^2/(6u). In a parity-breaking case it is -pi T^2(1/u_R+1/u_L)/12. The C51 oscillator proof establishes these formulas for the Gaussian theory; identifying an ultraviolet-regularized reference functional including its scalar is an additional step.

## 3. A four-history estimate from one-step derivative bounds

[proved conditional theorem] The following states sufficient local hypotheses and proves their propagation. It identifies the missing RG estimate more finely than postulating P itself.

For simplicity take q_0=E_*=1. General E_* is restored by dimensional rescaling. Let D_h be domains in Banach spaces B_h and let R_h^T:D_h->B_(h+1) be twice continuously differentiable integration/localization/rescaling maps. Both theories use this same map; their different dispersions, endpoints and fixed parameters must be encoded in their states. Suppose all line segments and two-parameter rectangles used below lie in D_h. This is a uniform neighborhood hypothesis, not a conclusion from small marginal couplings alone.

Write

    X_(h+1)^T=R_h^T(X_h^T),   Y_(h+1)^T=R_h^T(Y_h^T),
    z_h=X_h^0-Y_h^0,
    a_h=Y_h^T-Y_h^0,
    d_h=(X_h^T-Y_h^T)-z_h,
    eta_h=A 2^(-theta h),     tau_h=(T/q_h)^p.                  (4)

Assume 0<theta<p-2 and ||z_h||<=eta_h on every zero-temperature
scale h>=0. Assume also, uniformly for h<H,

    ||DR_h^T|| <= L,          ||D^2 R_h^T|| <= M,
    ||R_h^T-R_h^0|| <= b tau_h,
    ||D(R_h^T-R_h^0)|| <= b tau_h,
    L < 2^(p-theta).                                           (5)

The norms of maps and derivatives are supremum norms on the required common domains. The initial conditions are

    ||a_0|| <= a_* T^p,       ||d_0|| <= d_* eta_0 T^p.         (6)

Then there are constants A_*,D_* independent of T,H,h such that

    ||a_h|| <= A_* tau_h,
    ||d_h|| <= D_* eta_h tau_h,                  0<=h<=H.       (7)

[proof] The reference difference obeys

    ||a_(h+1)|| <= L||a_h||+b tau_h.

Since tau_(h+1)=2^p tau_h, choose

    A_* >= max(a_*, b/(2^p-L)).                                (8)

Induction proves the first estimate. Notice that the reference thermal change is not required to have an irrelevant factor eta_h.

For the mixed difference abbreviate y=Y_h^0, z=z_h, a=a_h, d=d_h and R=R_h^T. Add and subtract R(y+a+z) to obtain the exact decomposition

    d_(h+1)
      = R(y+a+z+d)-R(y+a+z)
        +R(y+a+z)-R(y+a)-R(y+z)+R(y)
        +(R_h^T-R_h^0)(y+z)-(R_h^T-R_h^0)(y).                (9)

The first line has norm at most L||d||. The rectangular second difference equals

    integral_0^1 ds integral_0^1 dr
                         D^2 R(y+s a+r z)[a,z],              (10)

so it is bounded by M||a|| ||z||. The final line has norm at most b tau_h ||z||. Thus

    ||d_(h+1)|| <= L||d_h||+(M A_*+b) eta_h tau_h.            (11)

But eta_(h+1)tau_(h+1)=2^(p-theta) eta_h tau_h. Choose

    D_* >= max(d_*, (M A_*+b)/(2^(p-theta)-L)).                (12)

Induction proves the second estimate of (7). No temperature-dependent retuning of a relevant counterterm was used. Its possible one-step amplification belongs in L. For example p=3, theta<1 permits the usual relevant amplification near 2, provided the full map has the stated norm bound. QED.

[proved scalar consequence] In addition suppose

    ||DE_(h,T)||+||D^2 E_(h,T)|| <= C q_h^2,
    ||D(E_(h,T)-E_(h,0))|| <= C q_h^2 tau_h.                  (13)

Require the zero-temperature derivative bound ||DE_(h,0)||<=Cq_h^2
at every h>=0, including scales below the thermal cutoff.

Apply exactly (9) to E in place of R. Its three contributions are bounded by

    C q_h^2 [||d_h||+||a_h|| ||z_h||+tau_h||z_h||]
       <= C' q_h^2 eta_h tau_h.                              (14)

This proves P on the active scales from (5), (6), (13). It exposes the term involving the Hessian and the reference thermal history which is lost if only the lattice/reference distance at T=0 is controlled.

[proved terminal and tail consequences] Suppose Q_(H,T) is Lipschitz with constant Cq_H^2. Equation (7) and tau_H<=1 give

    |Q_(H,T)(X_H^T)-Q_(H,T)(Y_H^T)|
      <= C q_H^2 ||z_H+d_H|| <= C' q_H^2 eta_H.               (15)

The zero-temperature shell difference is bounded by Cq_h^2 eta_h using (13). Use its negative as delta_h for h>H, and combine the h=H zero-temperature difference with (15). Since 1/2<T/q_H<=1, this combined terminal contribution also obeys P up to a fixed factor 2^p. For h>H, q_h<T and the minimum in P is one. If the ultraviolet four-term difference in a_X-a_Y is O(T^p), the exact identities (2) now give C51's full representation. Its theorem yields

    [f_X(T)-e_X]-[f_Y(T)-e_Y] = O(T^(2+theta)),               (16)

because p>2+theta. More general exponents or logarithmic rates can instead use C51's summation theorem. No new proof of that summation is needed here.

### 3.1 A version with history dependence

[proved conditional extension] Replacing a complete state by a list of couplings does not eliminate memory. Suppose the derived mixed recurrence is instead

    ||d_(h+1)|| <= sum_(k=0)^h A_(h,k)||d_k||
                            +B eta_h tau_h,   A_(h,k)>=0.   (17)

Set w_h=eta_h tau_h and assume

    sup_h sum_(k=0)^h A_(h,k) w_k/w_(h+1) <= r<1.            (18)

Then ||d_h||<=D w_h follows by induction with

    D>=max(||d_0||/w_0, B 2^(-(p-theta))/(1-r)).              (19)

Indeed dividing (17) by w_(h+1) makes its right side at most rD+B2^(-(p-theta)). This proves the claim. A corresponding bound with weights tau_h is required for a_h. Equation (18) is an explicit norm condition to verify for the causal history kernel. It cannot be replaced by the phrase "short memory" without checking the weights.

## 4. What is available for the one-step hypotheses

[read, imported local theorem] The collaborator's winding theorem compares the scalar produced by a single smooth shell at fixed running data and fixed covariance. Its hypotheses are a Gram determinant bound and weighted L1 covariance/endpoint bounds of sufficiently high order. The proof unfolds only the BBF spanning tree, retains the loop determinant, and bounds long configurations by a time moment. Its Section 7 proves the derivative version by Cauchy's estimate. Thus the fixed-data interaction part of (13) has a concrete proof under those local norm hypotheses. Its external-field extension supplies the analogous starting point for (5).

[proved limitation] These one-step estimates consume a reserve of field weight and coordinate moments. A fixed loss applied once per shell is not an infinite-scale theorem. To apply Section 3, one must prove that the renormalized tree expansion supplies uniform domains D_h and all the derivative bounds in a compatible family of norms. Also the winding theorem does not itself retain every Gaussian covariance-change scalar. Those factors must separately satisfy (13) and (15).

[read, imported localization lemma] The collaborator's four-node Matsubara projection extracts the quadratic constant/linear terms using frequencies +/-pi T and +/-3pi T. Its error in normalized linear coefficients is O((T/q)^3), with five moments sufficient for the stated spatial coefficient. The quartic projection has a fourth-order thermal error. This provides a possible p=3 scheme at one step.

[proved check of why the scheme matters] The two-node derivative estimator

    B_T f = [f(pi T)-f(-pi T)]/(2 pi T)

gives B_T f-f'(0)=(pi T/q)^2 for f(omega)=omega^3/q^2. Thus high smoothness by itself does not give p>2 for this localization convention. One must improve the projection or prove cancellations between its scalar and coupling contributions. This does not invalidate a physical thermal coefficient, which is independent of localization conventions. It prevents identifying a particular shell bound without specifying the convention.

[assumed, outstanding XXZ application] For ordinary weak XXZ, the BM estimates quoted in SOURCE_AUDIT give control of marginal histories and their limiting coupling. I have not proved their extension to the full common-domain, double-difference bounds (5), (13), with scalar Gaussian factors and the final thermal shell. A statement that the beta function vanishes does not supply those bounds. In particular, the interacting fermionic reference is not Gaussian in the fermion fields. The reference quartic must stay inside the comparison even though its current representation is quadratic.

[proved boundary] Section 3 is consequently a reduction of the missing estimate, not a verification of its hypotheses for XXZ. The single-scale results and zero-temperature flow estimates cannot be combined by multiplying their bounds, as Section 2 demonstrates.

## 5. The retarded correction and what it would change

[read] C55's available progress treats off-half-filled initialization and the linear response to arbitrary-degree retarded endpoints. It does not yet claim the uniform nonlinear all-scale theorem. The collaborator's SINGLE_CHANNEL_RG and UV_AUDIT establish finite-step initialization, with the stated norm qualifications. C48 and the operator-valued massive-matching notes supply candidate thermal kernels and their scalar baseline on their specified ensembles. These are separate imported results, not applications newly verified by C57.

[proved conditional modifications] Suppose the massive matching gives a zero-temperature endpoint R with exponential time rate kappa>0 and a thermal discrepancy of norm at most B exp(-kappa/T). Then for every fixed p and N,

    exp(-kappa/T) <= C_(p,N,kappa) T^(p+N),    0<T<=T_0.     (20)

To prove it set x=kappa/T. The ratio to T^(p+N) is kappa^(-p-N) x^(p+N)e^(-x), bounded on x>=kappa/T_0. Thus the matching discrepancy can meet initial conditions such as (6), and its scalar baseline is o(T^2), provided the retained integration has only the polynomial amplifications already bounded in Section 3. No infrared stability follows from (20) alone.

The zero-temperature R shifts the quadratic coefficients, forward coupling and physical field map. These shifts must be matched before defining z_h. A six-leg endpoint can contract to a quartic at first order in its coefficient, as C50 proves. Its shift is not generally second order. A residual frequency derivative in the quadratic term also changes the physical speed; the collaborator's exact quadratic-memory example proves this without an interacting assumption.

After this rematching, exponential endpoint weights supply any fixed finite number of time moments. The common nonlinear renormalized-tree estimate must still preserve them with the needed field-weight reserves. A bound on two moments alone is not the five-moment hypothesis used by the explicit Matsubara construction. If C55 establishes ordinary marginal/source control but not the mixed thermal comparison (5), (13), estimate (b) remains additional work.

Reflection symmetry of (1) enforces equal left/right speeds. A retarded correction that preserves U(1) and translation but breaks reflection need not preserve that equality. The safe coefficient then uses u_R and u_L separately. A thermal action also needs the correct physical trace normalization, or a matching theorem that supplies it. No scalar can be inferred from normalized correlators alone.

## 6. An unconditional calculation of the first interaction coefficient

[proved perturbative coefficient, not a fixed-V theorem] Define f_1 as the
thermodynamic limit of the first finite-volume Taylor coefficient,

    f_1(T,mu)=lim_(L->infinity) partial_V f_L(V,T,mu)|_(V=0).

At any fixed positive T, differentiation of the finite-volume trace gives the interaction expectation per site. The thermodynamic limit of its free-state Wick expression exists. With

    n_T = integral_(-pi)^pi dk/(2pi) n_F((-2t cos k-mu)/T),
    m_T = integral_(-pi)^pi dk/(2pi) cos k n_F((-2t cos k-mu)/T),
    n_F(x)=1/(1+exp x),

Wick's rule gives

    f_1(T,mu)=(n_T-1/2)^2-m_T^2.                             (21)

No interacting RG or bosonization is used in this identity. At zero temperature n_0=p/pi and m_0=sin p/pi. Denote a=p/pi-1/2 and m=sin p/pi.

### 6.1 The Sommerfeld step with its remainder

[proved] Near energy E=0 the two-root density of states and its cosine-weighted counterpart are

    D_0(E)=1/[pi sqrt(4t^2-(E+mu)^2)],
    D_1(E)=-(E+mu) D_0(E)/(2t).                              (22)

They are smooth on a fixed neighborhood of zero because |mu|<2t. Choose a smooth cutoff supported in that neighborhood and equal to one on a smaller one. Outside it the change of Fermi occupation is exponentially small; the band-edge square-root singularities are integrable. Inside it set E=Tx. The function n_F(x)-1_(x<0) is odd. Its first moment is

    integral_R x[n_F(x)-1_(x<0)] dx
       =2 sum_(l>=1) (-1)^(l+1) integral_0^infinity x e^(-l x)dx
       =2 sum_(l>=1) (-1)^(l+1)/l^2=pi^2/6.                (23)

The integrated absolute terms are summable, justifying the interchange. Taylor expansion through order two of D_i gives a remainder bounded by C|E|^3. The constant and quadratic terms integrate to zero against the odd occupation difference, while its third absolute moment is finite. Consequently

    n_T=n_0+(pi^2/6)D_0'(0)T^2+O(T^4),
    m_T=m_0+(pi^2/6)D_1'(0)T^2+O(T^4).                      (24)

Direct differentiation, using mu=-2t cos p, gives

    D_0'(0)=-2t cos p/(pi v^3),
    D_1'(0)=-2t/(pi v^3).                                   (25)

Insert (24)–(25) into (21):

    f_1(T,mu)-f_1(0,mu)
      = [2 pi t/(3v^3)] [m-a cos p] T^2+O(T^4).             (26)

This establishes the first Taylor coefficient of the thermal pressure correction, including the field-induced density shift. It does not exchange the interaction Taylor series with the T->0 limit.

### 6.2 Matching the dressed speed at first order

[proved coefficient calculation] At fixed particle density rho=p/pi the ground-energy Taylor coefficients are

    e^(0)(rho)=-2t sin(pi rho)/pi,
    e^(1)(rho)=(rho-1/2)^2-[sin(pi rho)/pi]^2.                 (27)

The first is the filled-band integral and the second is the free-state interaction expectation. Their derivatives yield

    kappa^(-1)=pi v+4V sin^2 p+O(V^2)                        (28)

as a formal Taylor expansion of the inverse compressibility. At this order the Fermi point is pi rho: the first-order change of a ground-state occupation expectation is zero at fixed particle count, since that occupation operator is diagonal in the free occupation basis. Equivalently the first-order self-energy only moves the Fermi boundary when the field is held fixed.

[imported first-order identity] Rescaling BM (16) and (18) to (1) gives

    K=1-2V sin^2 p/(pi v)+O(V^2),       u=K/(pi kappa).

Combining with (28) gives the fixed-density sound-speed coefficient

    u(V,rho)=v+V u_1^(rho)+O(V^2),
    u_1^(rho)=2 sin^2 p/pi.                                 (29)

Only the first Taylor coefficients of these identities are used. Their comparison with (26) is a consistency check on the proposed thermal formula, not a derivation of that formula from the susceptibility theorem.

At fixed mu, the ground-state stationarity condition mu=partial_rho e implies

    partial_V p|_(mu,V=0)=-2[a-m cos p]/v.                    (30)

Therefore the fixed-mu speed coefficient is

    u_1^(mu)=2 sin^2 p/pi
                   -(4t cos p/v)[a-m cos p]
             =(4t/v)[m-a cos p].                            (31)

The second equality uses v=2t sin p and m=sin p/pi. It follows that the coefficient in (26) is

    pi u_1^(mu)/(6v^2),                                     (32)

which is the derivative at V=0 of -pi/[6u(V,mu)]. Thus the first interaction correction agrees with thermal c=1 after dressing the sound speed.

### 6.3 Fixed density at positive temperature

[proved] For the canonical free energy, the derivative at V=0 is again (21), evaluated at the free chemical potential mu_T that gives n_T=rho. This follows either from the canonical free-state expectation or by differentiating the Legendre transform; its chemical-potential derivative cancels by stationarity.

Since partial_mu n_0=1/(pi v), (24) gives

    mu_T=mu+[pi^2 t cos p/(3v^2)]T^2+O(T^4).                 (33)

Also partial_mu m_0=cos p/(pi v). Substitution into (24) yields

    m_(T,rho)=m-[pi t sin^2 p/(3v^3)]T^2+O(T^4).

The density term in (21) is now fixed. Consequently

    f_1^(canonical)(T,rho)-f_1^(canonical)(0,rho)
       = T^2/(12t^2)+O(T^4)
       = [pi u_1^(rho)/(6v^2)]T^2+O(T^4).                  (34)

The coefficient's independence of rho at this order follows from the cancellation of sin^2 p against v^2. The expansion is not uniform near an empty/full band, where the permitted temperature neighborhood shrinks.

[proved checks] At half filling a=cos p=0, fixed-field and fixed-density coefficients coincide. At t=1/2, the canonical coefficient is 1/3. The coefficient has units energy^(-2), as required for the derivative with respect to V of a free energy divided by T^2. Replacing p by pi-p leaves (26) invariant, as particle-hole symmetry requires. The computation also recovers the free thermal coefficient from the two simple Fermi crossings, already proved in C46/C53.

[proved limitation] Equations (26) and (34) concern the first derivative at V=0. They do not bound the higher Taylor coefficients uniformly after division by T^2. No O(V^2 T^2) bound for the actual remainder is asserted. A nonperturbative thermal c=1 statement at any fixed V!=0 still requires the all-order comparison.

## 7. Several channels

[proved conditional extension] Section 3 is formulated in Banach norms and is unchanged when the state has finitely many branches. C50's covariance estimates show why fixed unequal velocities are allowed: shell bounds depend on inverse speeds and branch count, not inverse differences between velocities. The positive quadratic current Hamiltonian is diagonalized by C51 and the collaborator's MULTIVELOCITY_REVIEW. Under a complete scalar comparison, the reference contribution is

    f(T)=e_0-(pi/12) sum_a [1/u_(a,R)+1/u_(a,L)]T^2+o(T^2). (35)

The velocities in (35) are the eigenmode speeds of the dressed current Hamiltonian. They are not generically the separate fermionic band slopes or eigenvalues of an unnormalized stiffness matrix. The scalar theorem counts all these modes through their inverse speeds; an unweighted central charge additionally uses their current/stress identification.

[assumed additional RG work] At fixed finite channel count one must verify the full Fermi-point/counterterm map, detuning of every dangerous non-forward quartic, stability of the forward-coupling matrix, and the same mixed thermal estimates for the complete endpoint class. General higher-degree retarded interactions do not fall under the microscopic assumptions of the cited response theorems merely because such kernels are generated in those proofs. Constants may deteriorate with the smallest speed, detuning scale and number of channels. The limits must be taken separately at each fixed microscopic model.

[read, independent ensemble obligation] The physical-length transfer is C58's task. Its current progress distinguishes the positive completion of the bond, the sum over marker numbers, and fixed-density extraction. C57 proves none of these transfers. A conditional rooted-cell scalar result cannot be called a result per physical site without that argument.

## 8. Remaining estimate and deliverable boundary

[assumed, precise unresolved step] Construct one finite-temperature renormalized expansion for the lattice and the fully normalized dressed reference, on domains that persist through every shell, and prove the mixed flow and scalar derivative bounds (5), (13), including the terminal and Gaussian normalization terms. Alternatively prove the corresponding weighted history bounds (17)–(18) and scalar comparison directly in GN trees. The reference's full pressure must also be evaluated with the same normalization. This work has not verified those hypotheses even for instantaneous interacting XXZ.

[proved partial deliverable] The four-history decomposition (9) supplies the needed extra term and a closed bootstrap once the local estimates hold. The first-order calculation checks the interaction correction, field convention, physical velocity and canonical/grand-canonical distinction. The primary-source audit identifies an existing Ising estimate of the correct shape and a directly relevant quantum-chain announcement. These are useful reductions and checks, but estimate (b), the interacting thermal c=1 theorem, and the retarded/multichannel extensions remain unproved in C57.

## 9. Verification and provenance

[numerical] Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py`.
Direct quadrature of the free occupation integrals tests (26) and (34),
with a root solve for the canonical chemical potential. Three parameter
choices include the C53 field mu=1/3 at t=1/2, a second filling and hopping,
and half filling as a check. The final maximum coefficient discrepancies
are 1.35e-5 for fixed field and 2.53e-6 for fixed density. They decrease
quadratically in T, as the proved O(T^4) remainder predicts. These are
finite-temperature approximation errors, not failures of the identities.

[numerical] A nonlinear scalar flow with nonzero reference thermal forcing
tests the exact four-history identity and the normalized mixed bound. At
90-digit precision its largest identity residual is below 1e-91. This
is an algebra fixture, not an RG computation for XXZ.

[work record] The first numerical run stopped because the smallest
temperature in one fixture left a 5.37e-5 coefficient error above the
chosen 5e-5 threshold. Extending the temperature sequence resolved this
at the predicted rate. The initial binary-float history fixture also
failed when subtracting states whose difference was below floating-point
resolution. Decimal arithmetic resolved that cancellation. The final
run passes all assertions and is saved in checks.json. No lattice
simulation, cluster job, or interacting central-charge measurement was
performed. All writes stay in work_C57. No subagents were used and no
person was contacted. The local STATE.md is the handoff because the
parent project STATE.md is outside the permitted write scope.
