# C65 result

[proved, regulator obstruction] The requested uniform estimate fails for C62's literal sharp fermion windows in the fifth-moment kernel norm of TP, before the windows are removed. At every finite physical cutoff Lambda and every T>0, the static chiral density bubble contains

    Pi_Lambda^T(p)
      = tanh(A Lambda/(2T))/(2 pi A)
        - n_Lambda(1-n_Lambda)|p|/(2 pi T) + O(p^2),
    n_Lambda=1/(1+exp(A Lambda/T)).

The coefficient of |p| is nonzero. The cusp survives smooth infrared subtraction and occurs in an actual connected order-g^2 four-field kernel. Fermionic antisymmetrization does not cancel it. A finite first spatial absolute moment would imply continuous differentiability in p, so even that weaker moment fails. Consequently the thermal replacement cannot satisfy a finite bound with (T/q_h)^p in TP's norm for these finite regulators, for any p. See [PROOF.md Sections 2 and 3](PROOF.md#2-the-finite-temperature-sharp-window-bubble).

[proved, scope] This refutes a uniform-in-volume-and-regulator moment estimate, not TP after removal of the regulator or the thermal formula. C62's order of trace limits, removing the fermion window at fixed L,T first, remains valid. The edge coefficient tends to zero as the cutoff is removed. Its exponential smallness at a finite cutoff cannot compensate an infinite absolute spatial moment.

[proved, smooth trace-compatible repair] [PROOF.md Section 5](PROOF.md#5-a-smooth-quartic-regulator-with-the-same-trace-limit) constructs a smooth regulator of the same Hamiltonian. A sea-preserving fermionic attenuation map multiplies each vertex leg by r(k/Lambda), keeps the interaction quartic, and preserves C62's relative form bound with the same constants. Spectral compactness, min-max and a dominated Gibbs sum prove

    Tr exp(-beta H_Lambda^sm) -> Tr exp(-beta H_ref)

at each fixed L,beta, including convergence of the ground energy. The exact finite-active-mode Duhamel identity retains all scalar factors. Thus this smooth regulator has C62's trace-defined thermal coefficient -pi/(6v).

[proved, what smoothing fixes] The smooth bubble has no cusp. The regulator-dependent part of the offending static coefficient now has finite spatial moments, uniformly in Lambda, with a thermal error smaller than every fixed power of T at a fixed preparation scale. Smooth infrared shells also have cutoff-independent Gram and weighted covariance bounds. The collaborator's one-step winding and four-node localization estimates can therefore use p=3 under their stated endpoint hypotheses. This is a repair of the identified estimate, not an all-depth comparison.

[proved, BM alternative boundary] BM's smooth fermion cutoff removes the same moving-boundary cusp. Smoothing only the interaction does not. BM's spacetime-smoothed reference still needs its own full scalar identification or an o(T^2) scalar comparison. No such identification is proved here. The new smooth Hamiltonian regulator avoids that additional scalar obligation.

[assumed, remaining task] The full joint normalized tree estimate, common critical domains, actual four histories, and compatible vacuum derivatives in [PROOF.md Section 7](PROOF.md#7-the-remaining-joint-estimate-after-this-repair) remain unproved for the smooth reference. Thus C65 provides an explicit failure and a trace-compatible repair, but does not close Lemma TP or establish the interacting XXZ thermal theorem.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. Controls cover the sharp bubble and cusp slope, the divergent first spatial moment, smooth-bubble differentiability, and attenuation identities and form bounds on a 16-dimensional CAR space with a 256-dimensional dilation. The attenuation residual is 2.23e-16. No interacting RG flow was computed.
