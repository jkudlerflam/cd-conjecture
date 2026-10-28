# Adversarial audit of the vacancy infrared step

Status: **[physics inference, conditional on the frozen-chain and local
matching inputs]** The vacancy completion survives the targeted infrared
checks in this audit. I find no slow/fast instability that persists at
arbitrarily small positive J/t, no missing marginal locking operator in the
stated operator lattice, and no order-J charge dressing in the calculations
performed. These results support the conventional weak-coupling inference;
they do not independently verify the full microscopic c=24 conclusion.
The h23 phase and the uniform finite-J linked-sector matching remain inputs.
See [the response to the independent referee](REFEREE_RESPONSE.md) for the
distinction and correction to this report's initial stronger wording.

The dynamic smooth-density bound in
`contrib/vacancy_continuation/local_audit/README.md` is correct, with two important
scope restrictions: it is a **Euclidean Gaussian density** bound, not a real-time
bound and not a bound on the complete microscopic contact response. The
operator-valued interface argument can be strengthened beyond its scalar
version by explicitly contracting the connected internal correlation function,
as below. That is an actual infrared estimate, not an appeal to a missing
lattice-to-CFT theorem.

The remaining physical input is the usual local Wilsonian matching of a known
Luttinger liquid coupled weakly to finite-dimensional local degrees of freedom,
with the gapless Gaussian sector retained. There is no calculated universal
constant giving a numerically specified admissible J for the full 24-state bond.
No calculation here measures its central charge. The explicit second-order
and Gaussian estimates strengthen the local matching argument, but do not by
themselves establish its uniform finite-J completion. This is a substantive
remaining matching question even when a fully rigorous lattice-to-CFT theorem
is not being requested.

Audit by an OpenAI GPT-6-family Codex subagent, 2026-10-28. No human review.
Only this folder was written. The other chat's local_audit files were read and
not modified.

## Hypotheses actually used

1. The exact material representation is
   H=H_gap+J sum_j A_j h_j, with A_j=1[g_j=0] and unit-amplitude phase hopping.
2. The fixed h23 has a selected period-24 sector with a strictly positive
   material defect/interface margin. Its zero-defect low-energy theory has
   23 simple fermion bands, the stated quadratic and quartic momentum margins,
   and perturbations controllably small relative to its slowest band.
3. Occupied number and the internal U(1) are exact. A local neutral operator in
   the selected code has the usual even-fermion operator expansion. Marker
   defects are massive; they are not omitted gapless vertex operators.
4. Parameters g and lambda are fixed before choosing J/t. All constants below
   may depend, very badly, on that fixed h23 and on b=24.
5. The local matching has a remainder small in the normalized bulk couplings
   and per defect/interface, uniformly in system length and defect pattern.
   The calculations below support this hypothesis; they do not prove it at
   finite J for the full operator-valued microscopic coupling.

No argument here says that U(1) alone protects 23 bands. Without the certified
nonforward quadratic/quartic mismatches, pair transfer can be marginally
relevant and this construction would fail.

## Independent thermodynamic derivation

Let r=(L-N)/N and rho=1/(1+r). Free hard-core particles give the energy per
material gap

    e(r)=-(2t/pi)(1+r) sin(pi/(1+r)).

Direct differentiation and the physical-to-material wavevector conversion give

    e''(r)=2 pi t rho^3 sin(pi rho),
    v=2t rho sin(pi rho),
    K=pi v/e''=rho^-2.

Here the convention is Delta[exp(2i phi)]=K. Thus at rho=1/2:

    r=1, v=t, K=4, e''=pi t/4.

This agrees with the q=0 dressed charge and susceptibility in Bortz and Sergeev,
Section 4.1, equations (42)-(46), after multiplying their hopping 1/2 by 2t
and converting their K by pi. The primary source was opened independently:
[The q-deformed Bose gas](https://arxiv.org/html/cond-mat/0603093).

Wick's theorem for the physical charge sea yields

    alpha(r)=rho-sin(pi rho)^2/(pi^2 rho),
    alpha(1)=1/2-2/pi^2 = 0.297357632715...,
    a=d alpha/dr|_(r=1)=-(1/4+1/pi^2)=-0.351321183642....

The Gaussian action and density are

    S=(1/(2 pi K)) integral [ (d_tau phi)^2/v+v(d_x phi)^2 ],
    delta r=-d_x phi/pi.

Inverting the quadratic form gives

    G_rr(q,omega)=(K v/pi) q^2/(omega^2+v^2 q^2),
    0 <= G_rr <= K/(pi v)=4/(pi t).

The canonical q=0 density mode is absent. This bound is uniform in Euclidean
frequency, including omega of order J times any slow internal speed. It does
not have a hidden 1/J or 1/v_internal factor. The retarded continuation has
sound poles; it is false to bound its absolute real-time-frequency value by
4/(pi t). Keeping the coupled Gaussian fields, rather than treating those poles
as a perturbative instability, is essential.

For a static scalar source, with the convention
E^(2)=-J^2 sum_ij s_i chi_ij s_j, the smooth contribution is

    chi_ij=C(delta_ij-1/N),
    C=a^2/(2e'')=0.07857579749.../t.

The factor 1/2 comes from the one-sided energy denominator; G_rr is the
two-sided Euclidean response. The canonical pedestal contributes
J^2 C (sum s)^2/N. Changing k bounded sources changes this term by O(J^2 k/t),
even at k proportional to N. A translated period-b domain has the same mean
source; endpoint fragments change its sum only by O(b) per interface.

## What the new exact response calculation does and does not show

`response_check.py` independently builds the phase Hamiltonian using gap
compositions, obtains its ground state with sparse diagonalization, and solves
the spectral resolvent for both gap occupation and A_j. No earlier check code
is imported. At N material sites there are N holes, hence 2N physical sites.
Fourier operators are normalized by 1/sqrt(N). Units are t=1.

| N | dimension | G_rr(2pi/N,0) | G_AA(2pi/N,0) | G_AA(0,0) |
|---:|---:|---:|---:|---:|
| 4 | 35 | 0.96977208 | 0.19254833 | 0.04907381 |
| 6 | 462 | 1.09920890 | 0.19955340 | 0.04994507 |
| 8 | 6435 | 1.16098172 | 0.20275791 | 0.05024811 |
| 10 | 92378 | 1.19465244 | 0.20441890 | 0.05038805 |

The long-wavelength smooth-density value is 4/pi=1.27323954. The density
response at omega=q approaches half that value. The finite full contact
response is larger than a^2*4/pi=0.15715159, and its uniform canonical response
is nonzero even though the conserved density's is zero. This is expected:
the contact operator has UV and higher-dimensional contributions, in addition
to a delta r. Consequently the Gaussian bound must not be advertised as an
exact lattice bound for A. These extra contributions are finite local matching
coefficients, rather than evidence for a singularity. The saved data are not a
test of c=24 or of the full operator-valued h23 coupling.

## Slow/fast mixing and possible Wentzel-Bardeen instability

For one slow boson, use density stiffness

    B = [[v_c/K_c, g], [g,u/K_s]],
    A = diag(v_c K_c,u K_s)

for the conjugate current stiffness. Stability requires

    g^2 K_c K_s < v_c u.

For K_s=1 the mode speeds obey

    v_+^2+v_-^2=v_c^2+u^2,
    v_+^2 v_-^2=v_c^2 u^2-v_c K_c u g^2.

Thus there really is an instability if a finite g is kept while u tends to
zero. It is not enough to call the coupling marginal. In this construction,
however, u=J alpha v_* and g=J c with v_*>0 and c fixed by h23. The squared
relative coupling is O((J/t)c^2/v_*), which tends to zero. The raw mixed
coefficient is order J, while its normalized stiffness is order sqrt(J/t).
Both statements can hold simultaneously. The same conclusion follows from
the Schur complement of the full 24-mode stiffness matrix.

This is the familiar coupled-liquid instability mechanism, not a new one;
the primary discussion in
[Loss and Martin](https://arxiv.org/abs/cond-mat/9408086) was checked. The
equations above are independently derived from the quadratic matrices.

There is no order-J boundary energy hidden by using the charge UV scale t
as a cutoff for a static 1/u approximation. Internal dynamics occupy a
bandwidth O(J). At real charge frequencies well above that band the internal
response is not its static 1/u limit. Equivalently, the two-speed formula gives
v_+ + v_- - v_c - u = O(J^2/t) at fixed h23. Multiplying a static relative
correction by an unrestricted charge bandwidth would mix two incompatible
frequency regimes.

## Quantum sources and interfaces: an additional estimate

At second order around H_gap+J alpha H23, the charge and material correlation
functions enter as

    delta E_2 = -J^2 sum_ij integral_0^infinity d tau
        <delta A_i(tau) delta A_j(0)>_charge
        <h_i(tau) h_j(0)>_material.

Split the second factor into its static disconnected product and its connected
part. The disconnected part is exactly the scalar static problem: the
noncontact smooth contribution integrates to zero, its canonical pedestal is
handled separately, and dimension-two contact descendants have at worst a
|x|^-3 static tail. The first integer harmonic has dimension K=4 and gives
|x|^-7. These tails have a finite first spatial moment.

The connected quantum term must not be bounded using the charge kernel alone.
That would give an artificial 1/|x| envelope. At the free payload reference,
the lowest-dimension connected local even-fermion operators are bilinears of
dimension one. They include currents and nonforward bilinears. A local h_j
may contain the latter even though their uniform sum oscillates. Away from
the microscopic core, at this reference,

    |C_h(x,tau)| <= C_h/x^2,
    |G_rr(x,tau)| <= C_c/(x^2+t^2 tau^2),

up to fixed b-dependent constants and drift/velocity conventions. Hence

    integral d tau |G_rr C_h| <= pi C_c C_h/(t |x|^3).

The bound is uniform as a payload speed becomes small by reducing J: the
maximum over time of its correlator is its equal-time envelope. Interactions
can lower a nonforward bilinear's dimension below one. The actual condition
needed is only Delta_min>1/2. With relative Gaussian error eta<1/4,
Delta_min>=sqrt(3/5)>1/2, and the repaired bound is

    integral d tau |G_rr C_h|
        <= C/(t |x|^(1+2 Delta_min)).

Massive cell correlations and higher-dimension payload operators are less
dangerous. A boundary or interface permits at most O(|x|) linked pairs at
separation |x|, so the first moment converges because 2 Delta_min>1. This
supplies a connected quantum second-order boundary estimate of O(J^2/t),
not merely an extensive bulk norm. This refinement was identified by the
independent audit_vacancy worker; it removes an unnecessary claim that every
interacting local correlator retains its free |x|^-2 power.

For two Gaussian density currents of speeds v and u, the time integral can
even be evaluated explicitly (normalization constants omitted):

    integral_(-infty)^infty d tau
      [(x^2-v^2 tau^2)/(x^2+v^2 tau^2)^2]
      [(x^2-u^2 tau^2)/(x^2+u^2 tau^2)^2]
      = 2 pi u v / [(u+v)^3 |x|^3].

It vanishes linearly with u/v for u<<v, rather than diverging. The script
checks this integral for u/v=1,0.1,0.01,0.001. For operator insertions confined
to finitely many defect cores, the bounded Fourier kernel can instead be used
directly: the core support gives an O(number of cores) quadratic form. One
does not need to assume a critical bulk correlation bound for arbitrary
localized degeneracies.

This calculation is conditional on the usual Gaussian long-distance bounds.
It does not assert that all microscopic quantum sources are bounded by a
scalar static susceptibility matrix. Nor does a second-order calculation by
itself prove convergence of an arbitrary linked expansion. Its purpose is to
check the specific infrared failure mode that could invalidate ordinary local
matching. The only marginal bulk coupling involved is Gaussian and can be
resummed exactly; the other bulk charge operators have dimension at least two,
and the non-Gaussian uniform combined operators have dimension strictly above
two. There is no remaining marginally relevant bulk coupling with an unknown
sign that would spoil this estimate at an exponentially long scale.

## A safer route to the positive selector margin

At lambda=0, the vacuum-token sector is identity on the data bits. The g
conversions act only in excited-token sectors. Below the vacuum-walker gap,
the source in each selected cell is therefore a scalar periodic source,
with connected vacuum-walker correlations confined to that cell before
charge dressing. This makes selector stability simpler than the full critical
payload problem. Charge dressing of a fixed cell or interface starts at
O(J^2/t), with the preceding infrared controls. For any fixed finite list of
nearby competing cell lengths the unscaled walker gaps are positive.
Arbitrarily long competing cells require an extensive selector bound; their
constants cannot silently be treated as the selected-cell gap. This reduction
can be made concrete. Write delta_sel=2-mu_sel and e_b<0. Then

    E_l-l e_b >= |e_b| l-delta_sel.

For the b=24 midpoint,

    delta_sel=0.04610690275289...,
    e_b=-0.00126401272424....

Every l>=73 therefore has cost at least |e_b| l/2, or
0.00063200636212... per material site, before port corrections. An extensive
O(J^2 l/t) charge matching correction can be absorbed in that margin. Only
the finite list l<73 needs the per-cell comparison with fixed positive walker
gaps. This still uses local bulk subtraction for intervals of selected cells;
it does not turn an arbitrary O(J^2 L/t) error into a defect bound.

There is a useful quantum version of the extensive second-order estimate.
Include the finite UV/contact terms in a bound G_AA<=C_A/t for the full
Euclidean charge response; C_A is not the Gaussian coefficient 4a^2/pi.
Within the phase-model matching this is finite: the smooth part is bounded,
the remaining low-energy operators have integrable correlations, and their
UV matching coefficients are finite. The internal Euclidean spectral matrix
is positive and its frequency integral is its equal-time covariance. Therefore
the second-order spectral contraction on an l-site interval is bounded by

    |delta E_2| <= (C_A J^2/(2t)) sum_(j in interval) <h_j^2>
                <= C_A J^2 M^2 l/(2t),

with convention-dependent factors absorbed in C_A. This estimate uses no
internal excitation denominator at all and therefore no lower bound on the
walker gap as l tends to infinity. The assertion of a finite full-kernel
matching constant is a physics input, unlike the explicit exact Gaussian
coefficient; the finite ED data illustrate that the two must be distinguished.

Restoring ports costs at most J sqrt(3)|lambda| per port, because ||A_j||<=1.
Thus boundary self-energies from opening or closing finitely many payload
ports can be bounded coarsely by O(J|lambda|) without assuming anything about
their eventual critical boundary fixed point. This is smaller than the
J alpha gamma selector margin after lambda is chosen sufficiently small.
Relevant impurity flows or local degeneracies cannot defeat a finite
operator-norm boundary estimate. One must keep the gap chain intact while
cutting ports: cutting H_gap would introduce artificial order-t boundaries.
The marker pattern and grammar defects are conserved in material order, even
though their physical positions move with the particles. Consequently this
comparison is an energetic sector-selection question, rather than a tunneling
or quantum-proliferation problem for material domain walls.

This refinement addresses the concern that an order-one quantum boundary
self-energy in units of J has simply been assumed away. The exact port norm
bounds that contribution. The charge-induced part is controlled by a separate
small J/t ratio. A complete lattice inequality uniform in all defect patterns
is not derived here; the campaign's exact frozen-sector bounds plus local
matching supply the physics inference. In particular, proving only an
O(J^2 L/t) total energy error would still be insufficient.

## Why all 23 payload modes survive, subject to the band certificate

Physical marker harmonics must be included. If q is payload momentum per cell,
the neutral vertex dressed by an integer harmonic ell has

    k_phys=rho(q+2pi ell)/b,
    Delta=Delta_pay+[(q+2pi ell)/(2pi b)]^2.

At rho=1/2, physical uniformity imposes q+2pi ell=4pi b h and therefore

    Delta_uniform=Delta_pay+4h^2.

For h!=0 this is at least four. For h=0, the inherited quadratic and quartic
noncollision conditions leave only Abelian current bilinears and their
products, which are Gaussian. Higher neutral fermion vertices have free
dimension at least three. This includes pair transfer; checking only distinct
Fermi points would not suffice. The low-dimension detuning set is finite and
has the certified positive margin. A sufficiently small normalized stiffness
change keeps every allowed non-Gaussian uniform vertex above dimension two.
For example eta<1/4 gives the lower factor sqrt(3/5), so a reference dimension
three stays above 2.32. No exact conservation of each band number is required.

For the fixed full h23, let M bound local strengths, v_* the smallest unscaled
payload speed, Delta_* a relevant massive threshold, gamma_* the selector
margin, delta_k the finite momentum margin, and Omega_* a fixed nonforward
charge matching threshold. A sufficient schematic hierarchy is

    (J/t) C(h23) M^2 [1+|log(J/t)|]
       << alpha min(v_*,v_* delta_k,Delta_*,gamma_*),
    J M << Omega_*,
    ||G_ref^(-1/2) delta G G_ref^(-1/2)|| < 1/4.

All quantities on the right are fixed positive numbers before J is chosen.
The forward off-diagonal normalized norm tends to zero as sqrt(J/t); induced
internal diagonal terms tend to zero relative to J as J/t, with allowed logs.
There can also be order-J **charge** stiffness corrections. For example the
static background produces J<h> alpha''(r) (delta r)^2/2. Adiabatic response
of a machinery level with gap J Delta produces another coefficient of order
J/Delta, applicable only below that level's frequency. These are normalized
by the charge stiffness t and tend to zero as J/t at fixed h23. They are not
counterexamples to the O(J^2/t) statement for internal energy dressing; using
their adiabatic approximation all the way to charge frequencies of order t
would be invalid. This distinction was flagged by the audit_vacancy worker.
The earlier choices g then lambda remain essential. The hierarchy is very
restrictive but not contradictory. No finite-size charge gap is used.
Pre-existing payload interactions do not run across an artificial log(t/J)
interval: their entire bandwidth scales with J, and their matching interval
is the ratio of fixed h23 bandwidths and detuning scales. Charge-induced
corrections may contain log(t/J), but carry the additional vanishing J/t
factor. Thus taking J last does not undo the earlier small-g choice through
an unbounded payload-only RG time.

## Verdict and reproducibility

The specific infrared objections audited here do not refute the model. The
bounded Gaussian kernel, finite static first moment, connected quantum
first-moment estimate, exact port boundary norm, and explicit velocity
determinant give concrete support for the physical construction. The possible
reductions of the bare selector margin are assigned to independently reducible
O(J|lambda|) and O(J^2/t) terms at the orders controlled here. Completing this
assignment uniformly at finite J is part of the matching input. Subject to
the frozen band/sector hypotheses and that local matching input, the result
is 23 payload bosons plus the charge boson.

The precise unproved stronger statement would be a uniform microscopic
linked-sector inequality for the full operator-valued h_j, with a remainder
per interface tending to zero in units of J. Neither the Gaussian bound nor
the finite ED data supplies that statement. A rigorous theorem would be a
stronger validation; asking for a convincing physical derivation of this
uniformity is legitimate at the physics standard too. Concealing it behind an
extensive error bound would be a genuine gap in the argument.

Reproduce the finite checks with numpy/scipy:

    PYTHONPATH=/private/tmp/cd-campaign-python \
      /Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
      contrib/counterexample_audit_20261028/infrared/response_check.py

Default sizes are N=4,6,8,10. The run took about 1.2 seconds of measured model
construction/solve time on this host. Ground-state residuals were below
4e-13 and response solve residuals below 3e-10. The result file records all
finite-frequency responses and exact two-speed stability examples.
