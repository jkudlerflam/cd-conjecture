# C92 result

[Status] The full finite-penalty mobile defect comparison remains open. A uniform contact bound and the projected thermal pin comparison are proved. The exact remaining estimates and their conditional implication for the full pressure are in [PROOF.md](PROOF.md).

[Proved] For the legally projected gap Hamiltonian with positive length tilt s, arbitrary bounded quantum material sources of norm at most M, and c=s_min-JM>0,

    <A_j> >= min(1/2,c^2/(256t^2)),
    0<T<=c/(4 log 2).

This holds uniformly in volume, cell lengths and pin pattern. The proof removes one vacancy from a chosen nonempty gap through a two-Kraus channel. Its entropy loss is at most log 2. It applies to the length-tilted ensemble, not directly to a fixed-vacancy-number state.

[Proved, d=13] In the source interval t/2<=s<=3t/4 and for J/t<=2^-13, a conservative C88 bound is <A_j> >= 1/4096. Combining this with the referee's pin replacement and the entropy loss log(24576) per pinned length-15 cell gives

    F_pinned-F_unpinned >= k[delta_pin-T log(24576)],
    delta_pin>=J/4096

at the scheduled port strength. Therefore the pressure contribution of all length-15 pin flags is at most

    (T/15) log[1+24576 exp(-delta_pin/T)] = o(T^2).

This theorem includes the legal auxiliary states and nonzero port, allows arbitrary other cell lengths, and assumes infinite grammar penalty. It does not yet establish the finite-K pin comparison.

[Proved] Entire legal marker-free and one-marker sectors have smaller tilted pressure at sufficiently low temperature. This does not compare an island inside selected bulk. The all-vacancy class, including finite-K vacancy color walls, is also excluded at the half-density source.

[Assumed, remaining] Let W be the total material length of wrong-length cells plus the number of pinned length-15 cells. The missing legal charge comparison is

    R_sigma(T,s)-R_good(T,s) >= -J epsilon(J/t) W-o(N),
    epsilon(x)->0,

where R subtracts the same free gap trace and the full contact-averaged material trace. PROOF equations (27)-(29) give a concrete sufficient one-sided bound on the difference of their imaginary-time response functions. A positive contact lower bound does not prove that response estimate.

The full finite-K theorem additionally needs the exact relative bound

    (F_sigma^K-F_sigma^infinity)
      -(F_good^K-F_good^infinity) >= -J eta_K W-o(N).

All illegal words remain in these full sector traces. A local remainder bound proportional to total volume cannot replace this relative estimate. The selected finite-K denominator retains its virtual grammar excursions and common scalar bulk.

[Proved conditional] If these two estimates fit below the inherited selector margin, a length-weighted trace sum gives

    0<=p_full(T)-p_good(T)<=C T exp(-delta/T)=o(T^2),

with the thermodynamic limit first. The proof covers wrong lengths, pins, illegal words and marker-free competitors, with the stated source-window conditions. Its structure applies to d=14 and d=24 with their own selector and penalty constants. The new contact lemma applies to all these legal vacancy models.

[Numerical] All 48 quantum-source contact-channel checks and 16 vacuum-clock pin trace checks pass. The largest matrix is 336 dimensional. These are algebra checks at small sizes, not a microscopic thermodynamic calculation. See [CHECKS.md](CHECKS.md).

[Record] No network commands, subagents or external edits were used. The project STATE.md was left unchanged under the folder restriction. No unconditional full-ring thermal or central-charge claim is made.
