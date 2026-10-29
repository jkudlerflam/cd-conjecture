# C111 result

[Proved] The smooth-modulation lemma in specs/_smooth_lemma.md is false for the stated class. This remains true if the given bond is required to attain the gauge-optimized half-spread M. The counterexample does not violate the separate conjecture A <= C M f(d).

[Proved, construction] There is an explicit family at fixed d = 8 with M = 3/2. It is the tensor product of an XX chain and a four-state gapped spectator. For every fixed 0 < delta <= 1/20, the spectator has a unique product ground state and gap 2 delta. The uniform infrared theory is the usual XX sector, so A = 1 throughout the family. The XX infrared identification is the input supplied in the background; the spectator and energy estimates are proved directly in PROOF.md.

[Proved] With lambda = 1/10, a sector of the weighted spectator has the exact deficit

    D_delta(w) = sum_i [sqrt(delta^2 u_i^2 + lambda^2 v_i^2) - delta u_i],
    u_i = (w_i+w_(i-1))/2,       v_i = w_i-w_(i-1).

The full chain's deficit is at least D_delta. For a fixed nonconstant smooth positive periodic profile w_i = W(i/L),

    D_delta(w)/Q(w) -> lambda^2/(2 delta),
    liminf Def(w)/(M Q(w)) >= 1/(300 delta),
    Q(w) = sum_i (w_(i+1)-w_i)^2/w_i.

The length limit is taken first at fixed delta. Sending delta to zero afterwards disproves every finite C f(8), even with an o(Q) remainder.

[Proved] The same issue affects the modular geometry used to extract A. For w_j = j(N-j)/N on an open interval,

    Def(w) >= (1/(100 delta)) log N + O_delta(1),
    Q(w) = 2 log N + O(1).

Thus the unwanted contribution is at logarithmic order. It cannot be absorbed into a bounded boundary correction or a higher-gradient remainder. It comes from a gradient-induced field acting on a gapped degree of freedom, with a static response proportional to its inverse gap. Norm optimization does not eliminate it.

[Proved, proof details] The spectator has two local labels, A and B, each carrying a qubit. Its AA bond is -3/4 plus delta(XI+IX)/2 + lambda(ZI-IZ); its BB bond is the swap; its mixed-label bonds equal +I. All bond eigenvalues lie in [-1,1]. Equal-left/right-marginal states at eigenvalues -1 and +1 prove M = 1. Open-segment energy inequalities prove the uniform all-A product state is the unique ground state and has gap 2 delta. Adding XX gives M = 3/2 by the same marginal witnesses. PROOF.md gives the matrices, spectral bounds, gap proof, smooth expansion and logarithmic asymptotics in full.

[Proved, proposed routes] The three suggested approaches do not prove the stated lemma. A transformed trial state supplies the opposite variational direction. The exact static curvature in this example diverges as 1/delta despite bounded commutator norms. C105's interval obstruction already permits arbitrary lengths and overlaps, so choosing windows of length w/|grad w| does not evade it. Signed boundary corrections are outside that obstruction, but cannot certify the false universal inequality.

[Numerical] Small unrestricted checks passed on the counterexample, XX, Heisenberg, XXZ(1/2), SU(3), and four seeded random bonds. The full suite took 0.103 seconds, with largest matrix dimension 729 and largest eigenvector residual 3.08e-12. Scalar smooth-profile ratios at delta = .05, .01, .001, .0001 approach .1, .5, 5, 50, as proved. NUMERICS.md gives the values and distinguishes these tests from asymptotic results.

[Proved, uncompleted computation] Typhon returned "Connection timed out during banner exchange". The large Gaussian XX run is prepared in large_gaussian.py and cluster_job.sbatch but was not submitted or run. No heavy computation ran locally or on a login node.

[Assumed, remaining research question] A corrected route would need a physically specified energy-density representative and a controlled subtraction of nonuniversal local gradient response before identifying the deficit with A. No such general prescription or bound is proved here. The original velocity-weighted conjecture remains open for arbitrary infrared sectors.
