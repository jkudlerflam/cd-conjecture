# Holomorphic length tilts of the massive matching map

Status: a proved extension of the abstract thermal matching theorem, conditional on the stated commuting length bounds and on its application hypotheses. This controls the **massive influence action and its scalar density**, with the logical fields retained. It is not a theorem about the full gapless pressure or its compressibility.

Provenance: Codex `infrared_theorem`, 2026-10-28; independently checked by the separate Codex agent `alternative_attack`, not human reviewed. The [real-tilt embedding and physical ensemble comparison](../../thermal_ensemble_20261028/TILTED_EMBEDDING.md) are developed separately by the coordinating chat. This note isolates the complex-tilt and derivative statements to avoid duplicating that work.

## 1. Relative length bound

Let the abstract auxiliary operator a have vacuum Omega and excited gap g>0, as in [the matching theorem](INFINITE_CELL_PROGRESS.md). Let q be the auxiliary length charge: q Omega=0, and on an excluded original cell of length ell,

    q=ell-b.

Assume a and q strongly commute. Suppose that outside finitely many exceptional length sectors,

    a >= gamma (1+|ell-b|),       gamma>0,

and on the exceptional excited sectors a>=g and |ell-b|<=b+1. Then

    |q| <= C_a a,
    C_a=max((b+1)/g, 1/gamma).                          (1)

This is an operator inequality by the direct sum over length sectors. The hypotheses match the legal/well-formed rooted length decomposition used in C44/C46; they are not asserted for arbitrary malformed strings without its invariant-sector qualifications.

Alternatively, the weaker hypotheses in [the real-tilt embedding note](../../thermal_ensemble_20261028/TILTED_EMBEDDING.md), namely

    a >= g Q,       |q| <= b Q+a/gamma,

give (1) with C_a=b/g+1/gamma, because Q<=a/g. All conclusions below then hold with this alternative C_a; the maximum formula in (1) uses the stronger sectorwise assumptions stated above. Either choice is fixed independently of the logical bandwidth.

Define in either case

    s_* = 1/(2 C_a)
    a_s = a+s q.                                       (2)

For the stronger sectorwise choice, s_*=min(g/[2(b+1)],gamma/2).

Here s can be complex. On the strip |Re s|<=s_*, strong commutation gives

    Re a_s >= a/2 >= (g/2) Q,
    ||exp(-t a_s)Q|| <= exp(-t g/2),                    (3)
    ||exp(-t a_s)||_1
          =Tr exp[-t(a+(Re s)q)]
          <=Tr exp(-t a/2).                            (4)

The trace statement is used only when the right side is finite. The vacuum is unchanged. No restriction on Im s is needed for these estimates, since the factor exp[-it(Im s)q] is unitary and commutes with the positive factor.

## 2. Holomorphy with an unbounded charge

[proved] Since |q|<=C_a a and a,q strongly commute, D(a_s)=D(a) throughout the strip. The upper graph-norm comparison follows from ||a_s psi||<=(1+|s|C_a)||a psi||. The lower comparison follows pointwise in the joint spectral representation from |a+s q|>=Re(a+s q)>=a/2. Thus a_s is a closed normal operator on the same domain, and s maps to a_s psi affinely for every psi in that domain.

For each t>0, exp(-t a_s) is operator-norm holomorphic in the interior of the strip. One may check this directly in the joint spectral representation: differentiating k times gives

    (-tq)^k exp(-t a_s).

On any smaller strip its norm is bounded by the supremum of a polynomial in ta times exp(-c ta), with c>0. This is finite independently of the upper auxiliary spectrum. The heat-kernel version is trace-norm holomorphic after reserving part of the positive exponential for a trace-class factor. These facts supply the required analytic dependence; q is not inserted as a bounded source vertex.

In the polymer proof, every Duhamel segment retains a_s inside its semigroup. Equations (3)-(4) make all block, passive-line and Schatten bounds uniform in s. To remove a spectral cutoff, use the projections onto the vacuum and auxiliary energies at most M; these commute with a and q and increase strongly to the identity. Compress each bounded vertex by these projections. The compressed vertices converge strongly, as do their adjoints, and remain uniformly bounded; operator-norm convergence of arbitrary bounded vertices is neither asserted nor needed.

For a fixed Duhamel order on a fixed finite spatial support, the resulting products converge strongly and their vacuum matrix elements converge. For a thermal coefficient, positive time gaps give Schatten-class semigroups whose Hoelder exponents can be chosen using the total retained heat time. Spectral semigroup tails tend to zero in these Schatten norms. Strong convergence of the uniformly bounded compressed vertices, on both sides, then gives convergence of the trace products; equivalently one first approximates the Schatten factors by finite-rank operators. Zero time gaps occur on a simplex set of measure zero and use the same uniform bound. The bounds (3)-(4) and the time-simplex majorant dominate the coefficient integrals, so the fixed-order vacuum and trace coefficients converge by dominated convergence. Apply this argument also after summing the absolutely summable input support and logical-monomial labels.

Finally the connected-polymer majorant controls the discarded large-order/support tail uniformly on every compact substrip. This gives convergence in the stated weighted kernel norm. The same argument, or Cauchy estimates with a slightly wider compact substrip, makes that convergence locally uniform in s. Finite-cutoff holomorphy therefore passes to the massive matching map as Banach-space holomorphy. This is a coefficientwise/cluster convergence argument, not an operator-norm approximation of the original Hamiltonian.

## 3. Uniform thermal constants

If Tr exp(-tau_0 a)<infinity, choose a single tau large enough that

    tau >= 4L/g,       Tr exp(-tau a/2)<=e.              (5)

Then the thermal matching theorem applies uniformly throughout |Re s|<=s_* with auxiliary gap g/2 and this fixed tau. Use its constants

    kappa=L/(8tau),      v_*=1/(8tau),
    B=100 exp(4 kappa tau)/tau,

after any harmless enlargement of B reserved for periodization. The input support exponent L and field weight are also fixed. These constants have no dependence on the logical bandwidth W or on s within the strip.

For C46's candidate heat bound Tr_Q exp(-t a)<=C0 exp(-tg), t>=tau_0, a sufficient explicit choice is

    tau >= max(4L/g, 2tau_0,
               (2/g) log[C0/(e-1)]).                   (6)

The heat trace in the Schatten estimate is uniform as well: its retained half-segments give Tr exp[-beta Re(a_s)/2], bounded by Tr exp(-beta a/4). Since beta>=2tau, (5) controls it by e per participating spatial site.

## 4. The static logical tilt is extracted exactly

On the valid embedding the physical length operator has the form

    L_cell=b+N_len,log+q_aux.                            (7)

The code lengths b and b+1 have N_len,log=0 and 1; excited original states are embedded with logical vacuum and q_aux=ell-b. Therefore a length tilt is represented by

    a -> a_s,
    H_log -> H_log+s sum_j(b+N_len,j).                  (8)

The scalar and the logical number term in (8) are retained exactly. They are not included in the small massive perturbation norm. This distinction prevents an artificial factor s/W in the massive expansion.

If the local SW unitary commutes with total physical length, its transformed identity is

    U H(s) U^dagger = D+R+s L_total.

Thus the same D, R and embedding may be used for every s: the residual does not acquire an uncontrolled length insertion. C45's generators have this property when their stated total-length conservation holds, because their projections and unperturbed reduced resolvents commute with that charge and all perturbing terms preserve it.

With the extracted reference (8), K_D continues to annihilate the auxiliary vacuum from both sides. This statement holds coefficientwise for every complex s. The lemma does not select a common Jordan-Wigner convention or prove the remaining payload/length Gaussian dictionary.

## 5. Cauchy bounds for massive-action derivatives

Let the coefficient norms of K_D and R be k,r and suppose k+r<=v_*/4 uniformly in the strip. Define

    E_s,beta=W_s,beta[K_D+R]-W_s,beta[K_D].

The thermal matching theorem gives

    ||E_s,beta|| <= C_match r,
    C_match=4B/v_*,                                    (9)

for beta>=2tau and every |Re s|<=s_*. The norm includes nonconstant connected kernels and the scalar density with their respective conventions.

For real s with |s|<=s_*/2, apply Cauchy's formula on circles of radius rho<s_*/2 and then let rho increase to s_*/2. This also covers centers whose maximal circle touches the boundary of the open holomorphy strip. For every integer m>=0,

    ||partial_s^m E_s,beta||
       <= m! (2/s_*)^m C_match r.                      (10)

All constants are fixed before choosing the final coupling lambda. The same argument applies to the zero-temperature endpoint and to the thermal-to-periodized-vacuum difference. In particular,

    ||partial_s^m [W_s,beta-Per_beta W_s,infinity]||
       <= m! (2/s_*)^m C B exp[-beta/(2tau)].            (11)

The baseline thermal term W_s,beta[K_D] obeys the analogous exponential estimate. These bounds concern massive matching; no logical propagator has been integrated.

For a dimensionless low-energy tilt y=s/W, chain rule gives partial_y^m=W^m partial_s^m. Thus, after normalizing an action coefficient by W, the corresponding massive error derivative is bounded by

    m! (2W/s_*)^m C_match r/W.                          (12)

This formula does not include differentiation of any separately chosen lambda-dependent retuning path. Along the direct affine physical-length tilt, D and R are fixed as explained in Section 4.


## 6. A fixed-coupling schedule for the C45 bounds

[proved conditional on the quoted local SW bounds and the invariant-space application] Choose the physical support exponent large enough to pay the logical coefficient conversion before computing C45's v0,G,C_SW. With that choice the conversion factor is absorbed into the support weight. Since K_D subtracts a norm-contracting code compression, one may use the conservative bounds

    k <= 4 lambda v0+4 lambda^2 G,
    r <= C_SW lambda^3.                                (13)

After the wrapper, tilt strip, heat time, and output kernel norm have been fixed, it suffices to add the restrictions

    lambda <= v_*/(64v0),
    lambda <= sqrt[v_*/(64G)],
    lambda <= [v_*/(16C_SW)]^(1/3).                     (14)

They give k+r<=3v_*/16<v_*/4. Retain all of C45's other local SW and sector restrictions. If W=c_n lambda^2, c_n>0 fixed, and the intended dimensionless tilt obeys |y|<=Y, add

    lambda <= sqrt[s_*/(2Y c_n)]                       (15)

for Y>0. Then s=Wy stays in the inner holomorphic strip. To make the normalized massive action error at most epsilon, also impose

    lambda <= epsilon c_n/(C_match C_SW).              (16)

Every right-hand side is a fixed positive number at fixed n. No quantity is chosen as a function of ring length or inverse temperature. Enlarging the support exponent changes the finite SW constants but creates no circular dependence, since lambda is chosen last. This is an explicit consistency check of the scale ordering, not a practical numerical parameter certificate.

## 7. What the derivative bound does not imply

The theorem does not prove a beta-uniform complex analyticity disk for the **full logical partition function**. Even one free fermion level has partition function 1+exp[-beta(e+s)], whose complex zeros approach the real Fermi energy at distance pi/beta. Uniform massive analyticity is possible because those gapless logical fields remain external.

Likewise, (10) is not by itself a C2 bound on the full interacting pressure. Differentiating that pressure also differentiates logical expectation values and introduces gapless connected correlations. Controlling them, their compressibility, the Legendre transform at fixed physical length, and the cyclic quotient requires the subsequent logical RG/ensemble argument.

The real-tilt invalid-sector comparison is separate: for real s the modified auxiliary gap is at least g/2, and the additional bounded logical number term must be included in the cut-bond comparison or matched in isolated valid and invalid cells. Complex tilts are used here for analytic estimates, not for positivity or ground-energy comparisons.

No central-charge conclusion follows from this lemma alone.
