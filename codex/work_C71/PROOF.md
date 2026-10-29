# C71. Ultraviolet cancellation with the smooth trace regulator

## 1. Verdict and scope

[assumed, unresolved target] Lemma TP is not proved for XXZ at fixed nonzero small anisotropy. The smooth regulator removes C65's actual moment-norm counterexample. A new, explicit obstruction to applying the available absolute determinant estimate occurs earlier than the infrared joint comparison. The absolute two-line ultraviolet loop grows logarithmically with the regulator. Its integrated coefficient does not diverge. Distinguishing these statements is essential.

[proved partial results] Sections 3 and 4 prove this logarithmic obstruction and the cancellation that repairs the particular second-order density-loop coefficient. Section 5 proves a uniform thermal replacement for that coefficient in every fixed spacetime moment norm, with arbitrary fixed power greater than one, and with two speed derivatives. Thus the next task is an all-order ultraviolet preparation retaining these cancellations, followed by the infrared joint comparison. Merely repeating the smooth single-shell theorem does not perform either operation.

[assumed, imported construction] We use C65 Section 5's sea-preserving attenuation of C62's Hamiltonian reference. Its interaction remains quartic, has a regulator-independent relative form bound, and converges in the full Gibbs trace at fixed length and temperature. C62's matched parameters A,B give dressed speed v = sqrt(A^2-B^2). Their trace normalization and parameter conversion are not rederived. The reference interaction is

    V = g integral dt dx dy v0(x-y) rho_+(t,x) rho_-(t,y),
    g = 2 pi B,   Fourier(v0)(p) = chi(p).

Here chi is real, even, smooth, compactly supported in |p| < kappa, and equals one near zero. Hence v0 is Schwartz. The densities are sea normal ordered. C65 multiplies every fermion leg by r(k/Lambda), where r is smooth, even, lies between zero and one, equals one for |k| <= 1, and vanishes for |k| >= 2. Set w_Lambda(k) = r(k/Lambda)^2. Equivalently use smeared-field covariance

    G_Lambda(omega,k) = w_Lambda(k)/(A k - i omega).        (1)

[proved scope of the calculation] One chirality suffices. Reflection gives the other. The calculation concerns the reference used for off-half-filled XXZ, not a new lattice Hamiltonian. Constants below are uniform for A in a fixed compact subset of (0,infinity). They may depend on the fixed preparation scale, cutoff profiles and moment orders. They are not uniform near an XXZ band edge or the excluded filling.

We use Fourier phase exp(-i omega t + i k x). Bosonic transfer frequency is denoted Omega. The thermal exponent is denoted P to avoid confusing it with the transfer momentum p. The ultraviolet preparation scale is fixed, and can be set to one by rescaling.

## 2. Smooth high covariance and its actual ultraviolet singularity

[definition] Choose a real C-infinity cutoff f(omega,k), compactly supported, even under (omega,k) -> (-omega,-k), and equal to one near zero. Its spatial support lies in |k| <= K. Define

    u = FourierInverse[f/(A k-i omega)],
    a_Lambda = FourierInverse[(w_Lambda-f)/(A k-i omega)],
    a = FourierInverse[(1-f)/(A k-i omega)].               (2)

Take Lambda >= Lambda0 > K+kappa. Then w_Lambda f=f, and (2) is a genuine common low/high covariance splitting. Its thermal version is antiperiodization, with the same temperature-independent cutoff. None of these choices changes a physical counterterm with temperature.

[proved] For t>0, integrating the frequency in (1) gives

    c_Lambda(t,x) = integral_0^infinity dk/(2 pi)
                           w_Lambda(k) exp[-k(A t-i x)],
    c(t,x) = 1/[2 pi(A t-i x)].                           (3)

For t<0 both kernels are fixed by c_Lambda(-t,-x)=-c_Lambda(t,x). This also follows from the odd Fourier denominator. The equal-time prescription is inherited from the Hamiltonian trace. Values on the measure-zero equal-time slice do not affect the estimates below. In particular a_Lambda=c_Lambda-u and a=c-u.

[proved uniform singular bound] Write R=(t^2+x^2)^(1/2), absorbing the fixed speed constants in the bounds. For each fixed integer N, and j=0,1,2,

    |partial_A^j a_Lambda(t,x)| + |partial_A^j a(t,x)|
                  <= C_(N,j) R^(-1)(1+R)^(-N), R>0.     (4)

The bound holds uniformly in Lambda >= Lambda0. It is locally integrable but its square is not.

Proof for a. Decompose (1-f)/(A k-i omega) into smooth ultraviolet annuli with radii Q=2^j, j>=0, plus a fixed crossover. Rescaling and Fourier integration by parts give C_M Q(1+Q R)^(-M) per annulus, also after two A derivatives. When R<=1, sum Q<=1/R directly and Q>1/R with M>1. Each sum is at most C/R. When R>=1, every annulus has Q R>=1 and the sum is at most C_M R^(-M). Increasing M proves (4).

Proof of the cutoff assertion. Put z=A t-i x for t>0. From (3), integration by parts gives

    c_Lambda-c = (1/(2 pi z)) integral_Lambda^(2 Lambda)
                                      w_Lambda'(k) e^(-kz) dk.

All derivatives of w_Lambda' vanish at its support boundaries. Another N integrations by parts imply

    |c_Lambda-c| <= C_N |z|^(-1)(Lambda |z|)^(-N).       (5)

For |z|<=1/Lambda, use |c_Lambda|<=C Lambda and boundedness of u, which follows from local integrability of f/(A k-i omega). For |z|>=1/Lambda use (5) and the bound on a. Since Lambda is bounded below, these estimates give (4). Reflection treats t<0. Differentiating in A inserts t derivatives of the Laplace transform. The resulting factors |t|/|z| are bounded. Factors (Lambda |t|)^j in the edge integrals are paid by j additional integrations by parts in (5), using |t|<=C|z|. This proves the two derivative bounds as well. The compactly supported low symbol and its A derivatives remain locally integrable, so the same argument applies to u. QED.

[proved consequence] C65's smooth annular Gram bounds are compatible with (4). They do not give a uniform bounded Gram constant for the entire high covariance. The physical continuum has arbitrarily large fermion momenta as Lambda tends to infinity, in addition to the unbounded frequency variable. The bounded Brillouin-zone determinant argument in the collaborator's UV_AUDIT does not apply unchanged to this continuum limit.

## 3. A logarithm in the absolute ultraviolet tree bound

[proved] Fix delta>0. Plancherel in x, followed by the positive t integral in (3), gives

    I_Lambda = integral_0^delta dt integral_R dx |c_Lambda(t,x)|^2
      = (1/(4 pi A)) integral_0^(2 Lambda)
                   r(k/Lambda)^4 (1-e^(-2 A delta k)) dk/k
      = (1/(4 pi A)) log Lambda + O(1).                 (6)

At k=0 the numerator divided by k tends to 2 A delta, so there is no infrared divergence in this strip integral. On 0<k<Lambda the cutoff is one. Splitting that integral at any fixed positive k gives log Lambda plus a convergent constant and an exponentially small remainder. On Lambda<k<2 Lambda set k=Lambda s. The remaining integral tends to integral_1^2 r(s)^4 ds/s. This proves (6) without a sharp cutoff.

[proved] The same leading logarithm survives the low subtraction in (2):

    integral_0^delta dt integral_R dx |a_Lambda(t,x)|^2
                   = (1/(4 pi A)) log Lambda + O(1).    (7)

To verify this, expand |c_Lambda-u|^2. The spatial Fourier transform of u(t,x) is supported in |k|<=K and bounded by C(1+|log |k||), uniformly for 0<=t<=delta. This follows by bounding its compact frequency integral by integral d omega/(omega^2+A^2 k^2)^(1/2). The displayed logarithm is square integrable in k. Thus u has a finite strip L2 norm. The cross term, by spatial Plancherel, only uses |k|<=K, where w_Lambda=1. It is finite and independent of Lambda. Both corrections to (6) are therefore independent of Lambda once Lambda>K. This proves (7).

[proved relevance to the effective action] At order g^2, two right-density insertions contract to the chiral loop

    B_Lambda(z) = -a_Lambda(z)a_Lambda(-z)=a_Lambda(z)^2.

The uncontracted left densities give a connected four-field kernel. Before antisymmetrization its coefficient is, up to the cumulant sign and a fixed factor,

    K_Lambda(t,y-y') = integral dx dx'
            v0(x-y)v0(x'-y') B_Lambda(t,x-x')
                     = (phi *_x B_Lambda)(t,y-y'),
    phi = v0*v0,   Fourier(phi)=chi^2.                  (8)

Taking absolute values before the two internal spatial integrations gives, with y pinned,

    integral dy' dx dx' |v0(x-y)v0(x'-y')|
                              |a_Lambda(t,x-x')|^2
                = ||v0||_1^2 integral dx |a_Lambda(t,x)|^2.

Integrating over 0<t<delta and using (7) proves a logarithmic divergence of this majorant. The zeroth moment already diverges, so adding the fifth-moment output weight cannot repair it. A fixed field-weight reserve cannot absorb log Lambda at fixed perturbative order two.

In the two-vertex BBF representation of this contraction, one line is the connecting tree and one is a one-by-one residual determinant. Taking their magnitudes gives the same product. There is no larger determinant left in this contraction whose Gram estimate supplies the missing cancellation. Reorganizing by scales, exploiting Pauli differences, or evaluating a loop first may improve the bound. Those are additional operations, not consequences of (7) being small.

[proved limitation of the obstruction] This is a counterexample to this specific absolute-before-loop-integration bound. It is not an infinite norm of K_Lambda, not a proof of a divergent perturbative coefficient, and not a refutation of TP. Equation (8) includes internal integrations that must precede the coefficient norm. The next section proves that they remove the divergence for this coefficient. Unlike C65's sharp-edge cusp, (7) is an ultraviolet majorant obstruction, rather than a failure of the actual prepared moment norm.

## 4. The same coefficient after the loop is integrated

### 4.1 Exact dynamic bubble on the regulator plateau

[proved] The full finite-temperature chiral density response, before the low subtraction, is

    Pi_Lambda^T(Omega,p)
      = integral dk/(2 pi) w_Lambda(k)w_Lambda(k+p)
              [n_T(k)-n_T(k+p)]/(A p-i Omega),
    n_T(k)=1/(1+exp(A k/T)),   Omega=2 pi m T.          (9)

For nonzero denominator this follows by partial fractions of the two propagators and the Matsubara shift by the bosonic Omega. The product frequency sum is absolutely convergent at finite Lambda. Equivalently it is the imaginary-time free density correlator. At Omega=p=0 the static value is the continuous p limit. The distinct p=0, Omega!=0 values are zero, as charge conservation requires.

At T=0, n_0(k)=1_(k<0). For p>0 its difference in (9) is the indicator of (-p,0), and for p<0 it is minus the indicator of (0,-p). If |p|<kappa<Lambda, both cutoff factors are one on this interval. Therefore, as a distribution with the Hamiltonian equal-time convention,

    Pi_Lambda^0(Omega,p) = p/[2 pi(A p-i Omega)]
                                  on |p|<kappa.       (10)

This identity is uniform in all real Omega. It includes the static value 1/(2 pi A) and the dynamic vanishing at p=0, Omega!=0. The origin has the usual noncommuting static/dynamic limits; no continuity at that single point is claimed.

[proved] The high/high bubble is the full bubble minus the two cross bubbles with u, plus the u/u bubble. In each cross bubble one internal spatial momentum has |k|<=K. For |p|<kappa the other has magnitude at most K+kappa. Hence every regulator factor on these lines is one when Lambda>=Lambda0>K+kappa. Multiplication by chi(p)^2 makes the full transfer support lie in this interval. It follows from (10) that the smeared high/high coefficient K_Lambda^0 is exactly independent of Lambda for all such Lambda. This is stronger than pointwise convergence at zero temperature. It is special to the linear chiral bubble and is not an assertion about all higher coefficients.

### 4.2 Why its absolute spacetime moments are finite

[proved] The distribution fixed by (10) is

    B_full = -(i/(2 pi)) partial_x c,
    c(t,x)=1/[2 pi(A t-i x)].                          (11)

Indeed Fourier transformation of the right side gives p/[2 pi(A p-i Omega)]. Away from (0,0) it equals c(t,x)^2 by direct differentiation. Formula (11) specifies the short-distance extension, including its contact convention. An arbitrary additional delta function would change (10) and is not permitted.

Let eta(t,x) be smooth and compactly supported, equal to one near the origin. Since u is bounded and smooth there, the high bubble can be written

    B_high = eta[-(i/(2 pi)) partial_x c - 2 c u + u^2]
                                      +(1-eta)(c-u)^2. (12)

All products other than the displayed derivative are ordinary locally integrable functions. At large R the last term decays to arbitrary order by (4). Move the derivative through eta and then through the convolution in (8):

    K^0 = -(i/(2 pi)) phi' *_x (eta c)
          +(i/(2 pi)) phi *_x ((partial_x eta)c)
          +phi *_x [eta(-2 c u+u^2)]
          +phi *_x [(1-eta)(c-u)^2].                   (13)

Each function in brackets has finite weighted L1 norm with every fixed spacetime polynomial weight. Near the origin this uses integral_(R<1) R^(-1) dt dx < infinity. Away from it, use compact support or (4). The Schwartz functions phi and phi' have every fixed spatial absolute moment. The inequality

    1+|t|+|x+y| <= (1+|t|+|x|)(1+|y|)

then bounds the weighted norm of each convolution by the product of its two weighted norms. Thus for every fixed m,

    integral dt dx (1+|t|+|x|)^m |K_Lambda^0(t,x)|
                                              <= C_m, (14)

uniformly in Lambda>=Lambda0. Section 4.1 identifies K_Lambda^0 with (13), so no unjustified interchange of a nonabsolute bare loop integral is used. Two A derivatives satisfy the same estimate: derivatives of c remain homogeneous of degree minus one, derivatives of u are smooth locally, and derivatives of a obey (4).

[proved scope] Duplicating the two density-leg coordinates embeds (14) in the pinned four-field coefficient norm, with delta measures at the duplicated coordinates. Fermionic antisymmetrization is a finite sum of permutations and is bounded in the diameter norm. Consequently this actual generated quartic meets TP's fixed-scale fifth-moment requirement. This proves no bound on the rest of the prepared effective action.

## 5. A uniform thermal replacement for this coefficient

[proved] At finite Lambda, let a_Lambda,beta be the antiperiodization of a_Lambda. It is the high covariance sampled at fermionic Matsubara frequencies. For t in a period, suppressing x in the notation,

    a_Lambda,beta(t) = sum_n (-1)^n a_Lambda(t+n beta).

Using oddness under z -> -z, the high density bubble is its square. Subtract the periodic, bosonic periodization of the zero-temperature bubble. All equal-image terms cancel:

    B_Lambda,beta(t,x)-Per_beta B_Lambda(t,x)
      = sum_(n!=m) (-1)^(n+m)
                    a_Lambda(t+n beta,x)a_Lambda(t+m beta,x). (15)

At finite Lambda these are ordinary locally integrable products. After the subtraction, two singularities can never coincide. This observation permits a uniform absolute estimate that was impossible in (7).

[proved convolution estimate] For any prescribed moment order m>=0 and power P>1, (4) implies, for |R0|>=1,

    integral dt dx (1+|t|+|x|)^m
              |a_Lambda(t,x)a_Lambda(t+R0,x)|
                                      <= C_(m,P)|R0|^(-P). (16)

Here is an explicit bound. Put z=(t,x), e=(R0,0), and choose N sufficiently large in (4). On |z|<=|R0|/2, the second factor is at most C |R0|^(-N-1), and the weighted first factor has bounded integral. On the complement |z|>|R0|/2, bound the first factor by C |R0|^(-N-1). Use

    (1+|z|)^m <= (1+|R0|)^m(1+|z+e|)^m

and the finite weighted L1 norm of the second factor. The two contributions are bounded by C |R0|^(-N-1) and C |R0|^(m-N-1). Taking N>m+P proves (16). Norm equivalence between |z| and |t|+|x| changes only fixed constants. The argument is unchanged for a product of any two A derivatives through total order two.

[proved thermal estimate] Integrate (15) over t in [-beta/2,beta/2] with weight (1+d_beta(t,0)+|x|)^m. Set ell=m_image-n_image, with ell!=0, and unfold the n_image sum into the line. Circular time distance is at most |t+n_image beta|. The resulting bound is

    sum_(ell!=0) integral_R dt integral_R dx
          (1+|t|+|x|)^m |a_Lambda(t,x)a_Lambda(t+ell beta,x)|
                  <= C_(m,P) beta^(-P) sum_(ell!=0)|ell|^(-P).

The sum is finite. Convolution with phi in x is bounded in the same weighted norm. Therefore, for beta>=1, j=0,1,2,

    ||partial_A^j[K_Lambda^T-Per_beta K_Lambda^0]||_(m,beta)
                                            <= C_(m,P) T^P. (17)

All constants are independent of Lambda>=Lambda0. The norm is the pinned coefficient norm with circular time diameter. At a general fixed preparation scale q0 the dimensionless thermal factor is (T/q0)^P; constants and engineering powers are restored by rescaling. No factor q_h^theta at arbitrarily small infrared scale has been proved by this fixed-scale statement.

The estimate passes to regulator removal by the same uniform integrable cross-image bound. Equation (14) and boundedness of periodization also give uniform thermal fifth moments. Speed differentiation preserves (15), since the cutoff and periodization are temperature independent and the derivative majorants above are integrable. Varying the spatial interaction in its fixed Schwartz norms is handled by differentiating phi. Thus a finite difference between two speeds in the fixed interval gains their distance times T^P by the fundamental theorem of calculus. These are fixed-scale differentiated and marked bounds, not derivatives of a complete RG flow.

[proved checks] Taking P=3 and any 0<theta<1 satisfies P>2+theta for this coefficient. No loss of the thermal power is forced by its ultraviolet singularity. At T=0 the thermal difference vanishes. At g=0 the generated kernel vanishes. Its coefficient is proportional to g^2, while the logarithm in the absolute majorant has coefficient ||v0||_1^2/(4 pi A) before the fixed cumulant factor. The static response in (10) agrees with C65; the dynamic response vanishes for a conserved zero-momentum density. The inverse-speed dimensions in (6) and (10) agree.

## 6. What this does and does not establish for TP

[proved boundary] Smoothness repairs the sharp-edge cusp but does not turn the uncut continuum high covariance into a bounded-momentum covariance. Equations (6)-(8) isolate a failure of the direct absolute ultraviolet preparation estimate at a specified contraction and fixed order. Equations (10)-(17) repair that contraction by retaining its chiral cancellation before taking the coefficient norm. They do not repair every ultraviolet subgraph or prove a uniform analytic connected expansion.

[assumed, next ultraviolet lemma] Let U_Lambda^T be the exact effective action after integrating (2), including its scalar, in the C65 trace normalization. What is needed next is a construction of U_Lambda^T on a cutoff-independent small complex parameter domain such that

    sum_(l even>=2) H^l ||[U_Lambda^0]_l||_(M) <= C |g|,
    ||D_a^j(U_Lambda^T-Per_beta U_Lambda^0)||_(H',m)
                                         <= C |g| (T/q0)^P,
                 j<=2,   M fixed, H>H', P>2+theta,       (18)

with the scalar separately bounded per spacetime volume, appropriate coefficient dimensions understood at q0, and covariance/source derivatives included as needed by TP. If the free scalar is assigned separately, its known factor must be added back. At nonzero g, derivatives in a include speed and a fixed normalized parameterization of the interaction shape; the power |g| in (18) is not claimed for derivatives with respect to g itself. These bounds must survive removal of the time and spatial regulators in the specified trace order. They must permit the subsequent common slow-field matching. Equation (18) is an assumed target, not a result of this note.

An organization that simply estimates every ultraviolet two-line loop by the integral of its absolute product cannot prove (18), because of (7). One must retain cancellations within ultraviolet subgraphs, for example by integration before absolute values or by a suitable renormalized loop expansion. Such a reorganization must also handle open external fermion lines, arbitrary generated degrees, derivatives and vacuum roots. Exact Gaussianity of the full current partition function alone gives none of these off-shell kernel bounds. We have not proved that every possible proof of (18) must use one particular loop identity.

[assumed, remaining infrared lemma] Even after (18), the actual XXZ and reference histories must be matched in one complete normalized kernel space. At q_h=2^(-h), the missing bounds remain

    ||X_h^0-Y_h^0|| <= C q_h^theta,
    ||Y_h^T-Per_beta Y_h^0|| <= C(T/q_h)^P,
    ||(X_h^T-Y_h^T)-Per_beta(X_h^0-Y_h^0)||
                                  <= C q_h^theta(T/q_h)^P,          (19)

and the compatible complete scalar derivatives

    ||D E_(h,T)||+||D^2 E_(h,T)|| <= C q_h^2,
    ||D(E_(h,T)-E_(h,0))|| <= C q_h^2(T/q_h)^P.         (20)

Here periodization includes the specified finite-temperature localization identification. The stronger field and moment reserves cannot shrink once per scale. The common derivative domains must include the rectangles in C57's four-history identity. C62's critical-interpolation example rules out inferring them merely from control of two endpoint trajectories.

[proved logical qualification] Matching only the limiting reference after Lambda removal is not, by itself, finite-Lambda matching. For example, if a marginal coordinate of Y_(h,Lambda)^0 converges as h tends to infinity to g_*(Lambda), and the corresponding lattice coordinate converges to g_*, then a bound ||X_h^0-Y_(h,Lambda)^0||<=C q_h^theta at fixed Lambda forces g_*(Lambda)=g_*. This follows by taking h to infinity. Trace convergence at fixed L,T supplies no such equality. We do not assert a nonzero mismatch for C65's regulator, and the second-order bubble above supplies none. A proof must either match at each regulator with controlled limiting parameters, or first remove it in a uniform prepared-kernel topology. A scalar trace limit is not that topology.

[assumed, imported conditional implications] The collaborator's fixed thermal-mark allocation, C62's shared-short-memory scale sum, and C57/C59's four-history and scalar summations apply once their stated input estimates hold. They are not rederived here. Under (19)-(20), ultraviolet and terminal control, and the trace identification, they give

    [f_XXZ(T)-e_XXZ]-[f_ref(T)-e_ref] = O(T^(2+theta)).

C62/C59 evaluate the reference coefficient as -pi/(6v). No hypothesis newly established in this note completes that implication for interacting XXZ. No retarded-action or multichannel extension is claimed.

## 7. Verification and provenance

[numerical] The local checks.py controls the logarithmic slope in (6), the dynamic plateau identity (10), the chiral derivative identity behind (11), and a distinct-image thermal integral for an explicit smooth high covariance. The last is an illustrative covariance check, not a calculation of an interacting flow. Numerical controls supplement the proofs; they do not verify (18)-(20).

[proved work record] C57, C59, C62 and C65 were read, including scripts and saved checks. The collaborator's WINDING_ESTIMATE, THERMAL_LOCALIZATION, MOMENT_RESERVE and THERMAL_HISTORY are used only with their stated hypotheses. The lattice ultraviolet audit and the finite-cell and ensemble reductions are contextual inputs. No new microscopic hypothesis for the encoded chain is verified here. SOURCE_AUDIT.md records the limited primary-source check. No absence or novelty claim about thermal XXZ results is made.

[proved scope] All files are confined to work_C71. No subagents, shell network commands, large computation or contact with authors was used. The local STATE.md is the handoff under the write restriction. These arguments have not received independent or human review. The full uniform joint comparison remains unproved.
