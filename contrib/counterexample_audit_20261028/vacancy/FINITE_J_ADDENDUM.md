# Addendum: can matching constants secretly diverge as J tends to zero?

**Status: independent adversarial physics assessment, 2026-10-28.** This is an addendum to `REPORT.md`, not a replacement of its conditional theorem. It reviews the coordinating `../GAUSSIAN_BOUND.md` and the infrared worker's `../infrared/README.md` specifically for inverse-J factors and circular smallness assumptions. Corrections below were sent directly to that worker while its report was being revised.

**Conclusion:** the Gaussian matrix proof is correct. At fixed internal Hamiltonian h23, the proposed hierarchy has a conventional perturbative RG basin; the vanishing internal bandwidth does not itself force a matching constant proportional to 1/J. The concrete calculations support that physics inference, while a bound on one two-point kernel alone would not constitute an all-order microscopic matching theorem. No new refutation of the vacancy completion was found.

## 1. The Gaussian metric proof does what it claims

The chiral transformation and the identity GVG=Sigma V Sigma are correct. The geometric-mean integral has the correct normalization, including the scalar check sqrt(a b). Monotonicity proves the stated relative metric bound for arbitrary positive right and left velocities. Sylvester inertia supplies 24 modes of each chirality. This part neither assumes equal cone slopes nor hides a velocity denominator.

However, positivity alone is insufficient: it permits a charge umklapp dimension to fall below two. The quantitative relative neighborhood is essential. The proof explicitly states this, so it is not a circular proof of its own hypothesis.

## 2. Why the slow reference scale does not create a naked inverse J

Hold h23, its defect margin, its minimum payload speed v_*, and its massive threshold Delta_* fixed. Then the material reference Hamiltonian is J alpha H23. Equal-time material correlators and their amplitudes are independent of this overall J; time dependence is obtained by tau -> J alpha tau. That distinction matters when estimating a coefficient.

The smooth charge covariance satisfies

    G_rr(q,i omega)=(K t/pi) q²/(omega²+t²q²) <= K/(pi t).

Thus integrating its Gaussian field gives a quadratic influence kernel of order J²/t, uniform even at omega~J Delta_* or J v_* q. Rescaling material time tau'=J tau makes this correction order J/t relative to the material action. There is no substitution of a material gap for a charge denominator.

The full microscopic contact operator has finite UV contributions as well. Its response is not bounded by the smooth coefficient alone; the infrared worker's ED explicitly disproves that stronger numerical bound. Within the known phase-model infrared description, the additional contact operators have dimension at least two, while the first density harmonic has dimension four. Their low-frequency Fourier response is finite. Away from the soft momentum, the fixed charge threshold supplies the cutoff. The resulting constant may depend on the charge model and the fixed internal geometry, but the unperturbed charge model itself contains no J.

For a local massive internal transition, the charge's Ohmic spectral weight gives the schematic shift

    J² M² integral_0^O(t) [omega/t²]/[J Delta_*+omega] d omega
       = O(J² M²/t).

The small denominator occupies a correspondingly small spectral interval. Derivatives can produce logarithms, not an order-J shift. The associated orthogonality/spin-boson parameter scales as (JM/t)²; its logarithm of the transition frequency is harmless as J tends to zero.

In low-energy perturbation theory *after* charge matching, an induced internal vertex is O(J²/t). Repeated internal resolvents of size 1/(J Delta_*) then produce the pattern

    (J²/t)^k / (J Delta_*)^(k−1)
       = J (J/t)^k / Delta_*^(k−1).

At each fixed order this remains smaller than the leading material scale. Dependence on Delta_* is severe but fixed before J is chosen.

## 3. An O(J) charge stiffness correction is allowed and harmless

It would be incorrect to claim that every Gaussian coefficient changes only by O(J²/t). The disconnected material source already contributes

    J <h> alpha''(r) (delta r)²/2,

which shifts charge stiffness at order J. Likewise eliminating a massive transition adiabatically at frequencies below J Delta_* can yield

    (J delta r h_off)^2/(J Delta_*)
       = O(J/Delta_*) (delta r)².

These are charge coefficients and must be compared with its stiffness O(t); their relative size is O(J/t) with fixed h23. Mixed charge-payload current coefficients are O(J), normalized by sqrt(t J v_*), hence O(sqrt(J/t)/sqrt(v_*)). Pure induced material coefficients of size O(J²/t) are compared with J v_* and give O((J/t)/v_*).

Extending the adiabatic denominator 1/(J Delta_*) up to charge frequency t would falsely produce an O(J) material energy shift. The exact spectral integral above shows why that approximation is invalid. Retaining the low-frequency charge fields, and matching the high frequencies with their actual denominators, avoids the problem.

## 4. The apparent long RG interval does not multiply an unsuppressed payload coupling

The reference payload bandwidth and its low-energy matching cutoff both scale with J. Their ratio is fixed by h23. Therefore pre-existing weak material interactions run over a fixed, possibly enormous, interval determined by g and lambda; they do not acquire an additional log(t/J) merely because the charge bandwidth is t.

Between charge scale t and material scale J, the charge can produce matching logarithms. Their coefficients carry powers of J/t. Thus factors such as (J/t) log(t/J), or the smaller Ohmic (J/t)² log(t/J), vanish. A calculation that starts a linearized payload RG at charge bandwidth t would be using payload degrees of freedom outside their bandwidth and could manufacture the opposite conclusion.

The marginal intermode current interactions are retained and diagonalized as a Gaussian form. They are not assigned an unknown marginal beta-function sign. All remaining resonant non-Gaussian terms have a strict dimension margin; all dangerous nonforward terms have fixed momentum cutoffs. These are the specific features making the weak-coupling RG basin conventional. The smallness choice is a bootstrap around an established fixed reference and a strict neighborhood, not an inference that every marginal coupling is harmless.

## 5. Repair to the connected-interface correlation bound

An earlier infrared draft bounded every connected local even-fermion material correlator by C/x² because currents have dimension one. Local h_j can also contain oscillatory bilinears. In an interacting Luttinger liquid such a bilinear can have dimension below one even though its uniform Hamiltonian sum is momentum-forbidden. Thus x^−2 is exact at the free reference but is not the most general interacting bound.

The repair does not endanger interface control. Let Delta_min be a uniform lower bound on the dimensions of nonconstant neutral local material operators. The quantitative Gaussian neighborhood gives Delta_min>1/2 (at eta=1/4 it is at least sqrt(3/5)>0.774). Then

    |C_h(x,tau)| <= C_h/|x|^(2 Delta_min),
    integral d tau |G_rr(x,tau) C_h(x,tau)|
       <= C/(t |x|^(1+2 Delta_min)).

An interface has at most O(|x|) linked pairs at distance |x|, and the resulting first-moment sum converges because 2 Delta_min>1. Oscillatory phases can improve the estimate but are not needed. For exact current-current correlators the stronger x^−3 expression remains valid. Constants are uniform under the overall time rescaling J.

## 6. Arbitrarily long wrong cells need a separate bound

One cannot use the selected length-24 cell gap as a uniform gap for every competing cell length. A length-b' vacuum walker has a gap of order J/b'^2. Allowing a matching constant to contain its inverse and then declaring that constant fixed would be a genuine nonuniformity, since b' is unbounded in the thermodynamic problem.

The natural repair is the lambda=0 cell decomposition already identified by the infrared worker. Long wrong cells have a selector cost growing linearly in their length relative to the selected bulk. Their charge correction should be bounded by their volume, without dividing by their internal gap. At second order, a bounded Euclidean contact kernel together with the positive internal spectral covariance gives a bound proportional to

    (J²/t) sum_(j in cell) <h_j²>
       <= C J² M² b'/t.

The disconnected sources between cells are handled separately by the static kernel and its first moment. Restoring the ports costs a bounded J|lambda| per affected port. This supplies the appropriate power counting for long cells; it does not rely on their vanishing walker gaps. An all-pattern, all-order version is still the stronger microscopic linked-sector statement, and should not be claimed proved by this second-order estimate.

## Disposition

The concrete infrared estimates identify and remove the obvious inverse-J mechanisms. They support ordinary local Wilsonian matching at the requested physics standard. The statement that the full operator-valued microscopic remainder is uniformly small per defect remains stronger than the displayed two-point and Gaussian calculations. The correct assessment is neither that a 1/J divergence has been found nor that the Gaussian inequality alone proves the lattice counterexample. It is a controlled conditional existence argument with an explicit conventional matching input, now checked against the shrinking machinery gap, shrinking payload bandwidth, and long-cell nonuniformity.
