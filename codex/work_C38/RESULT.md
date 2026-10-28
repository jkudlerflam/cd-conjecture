# C38. All-n nonresonance from the scale hierarchy

[Proved, verdict] The proposed hierarchy argument does not establish all-n nonresonance for either the 23-letter frozen chain or C32's mobile chain. Its first missing implication is that a large bandwidth ratio makes cross-flavor effects small relative to the narrow band. That implication is false, including in the natural short-path Gram limit of these buses. This does not rule out an all-n parameter schedule for the actual families.

[Proved, main obstruction] In that limit the leading Gram matrix is

    Gamma=a^(-3) S [C+o(1)] S,   S_jj=a^(-2j),
    C_ij=(2i+2j+2)!/
         [(2i+1)!(2j+1)! 2^(2i+2j+3)].

For the first two sites,

    C=[[1/4,1/8],[1/8,5/32]].

The cross term changes the slow coefficient from 5/32 to the Schur pivot 3/32. Its relative contribution is 2/5, even as the energy ratio tends to zero. Adjacent squared correlations are (4j+4)(4j+3)/[(4j+6)(4j+5)], which approach one. The derivation covers the frozen large-mass stub and the positive mobile stub, after taking their existing large-express-energy Gram limit. See [PROOFS.md, Section 2](PROOFS.md).

[Proved, general mechanism] If h(k)=S F(k)S with invertible S, its Fermi zeros are those of F(k), independently of the scale separation in S. A two-flavor cross term sqrt(t_1 t_2)e produces a slow Schur correction t_2|e|^2/f_1. Dividing by the slow scale removes the bandwidth ratio. The fast flavor cannot be eliminated globally as a massive field because it also has Fermi points.

[Proved, loop-specific obstruction] In the real-hopping subfamily write

    h(k)=epsilon I-2S cos k+i[A+2K sin k],
    S=S^T>0,   A=-A^T,   K=-K^T.

Discarding K in the Hermitian Gram approximation makes h(k) even in k, including at nonzero A. The centers then coincide. For S=diag(s_a), the mixed small-A, small-K center response is

    theta_a=(2/epsilon) sum_(b!=a) A_ab K_ab/(s_a-s_b)
            + higher perturbative orders.

It depends on cross-scale matrix elements and their cancellations, not just on the separated eigenvalues s_a. A bound or sign theorem for these center differences and the remaining pair-sum functions is missing for general n. More general allowed phase choices require analyzing their full pencil instead of this restricted formula.

[Proved, schedule audit] C32's crossing-only choice ties t_d to g^2. The available finite-conversion error bound and an upper bound on loop-induced pair detuning can consequently both scale as g^2 in the normalized symbol. Reducing g alone does not establish that the former is smaller. One needs a symmetry-sensitive error estimate, or must fix a certified leading loop target before shrinking g. This is a failure of the current estimates, not a proof of resonance at finite g.

## What can be proved constructively

[Proved] A lacunary prescription does solve the complete arithmetic problem if imposed on all 2n endpoints. With r=1/16 and b=1/64, take

    p_(a,L)/(2pi)=-1/4+b r^(2a-2),
    p_(a,R)/(2pi)=+1/4+b r^(2a-1),  a=1,...,n.

Every distinct pair of endpoint sums, including repeated endpoints and modular umklapp, is separated by at least

    d_n=(7pi/256)16^(-(2n-1)).

The centers have ratio r^2. The proof cancels shared endpoints and bounds the largest surviving offset against the at most two opposing offsets. Modular wraparound is treated separately. This is stronger than checking distinct centers alone.

[Proved] The diagonal reference h_a(k)=2t_a[cos k_a-cos(k-theta_a)] realizes these points for any positive t_a, including arbitrary geometric energy scales. For a general Hermitian quadratic correction R(k), a sufficient condition is

    max_k {||t^(-1/2) R(k) t^(-1/2)||,
           ||t^(-1/2) R'(k) t^(-1/2)||}
       <= d_n^2/100000.

A Schur-complement proof retaining every off-diagonal entry then gives exactly 2n simple zeros, shifts at most twice the left-hand side, margins at least d_n/2, and velocities at least min(t_a)/2. In an ordinary norm it suffices to bound R and R' by min(t_a)d_n^2/100000. The accuracy is imposed by this budget; it does not follow from scale separation.

[Assumed, missing implementation] Neither candidate core has been shown to realize this endpoint prescription or satisfy this comparison bound for every n. C26 provides arbitrary target control in its abstract reflection cell, but importing that control does not supply a fixed-alphabet spatial realization. No new bond or enlarged fixed d0 with the required readout was constructed here.

## Frozen and mobile outcome

[Proved, frozen] C33 now reports a certified n=23 target, and C37 reports certified targets for every n=1,...,23. These are inherited finite results, not repeated here. C31 already proves that an energy hierarchy plus a common field can remove non-pair resonances of a centered reference. Its pair-transfer exception is the same missing step encountered here. No all-n frozen nonresonance result follows from C38.

[Proved, conditional mobile arithmetic] For the endpoint target above, every length density in [1/5,3/10] has mixed relevant/marginal momentum margin at least 59pi/160 in C20's cell-coordinate convention. One can choose the irrational density 1/4+sqrt(2)/100. Thus the length mode creates no further arithmetic obstruction for this designed target.

[Proved, mobile cross terms] C32's four residual length-pair bonds expand into constant, left-occupation, right-occupation, and occupation-product terms with norms at most e,2e,2e,4e when max||E_pq||<=e. This includes all scalar length interactions. A sufficient total bound is 9e. C32 can reduce e by choosing w/U first and s/delta second, using the actual small charge tunneling gap. The wrapper preserves a supplied target to any prescribed local tolerance; it does not create a missing payload nonresonance margin.

[Assumed, thermodynamic transfer] Turning those mobile operator bounds into dressed Fermi-point bounds needs local massive-sector elimination and continuum matching. Gaussian stiffness, oscillatory-flow control, and physical gluing are also inherited inputs. C40 is treating the RG budget separately. Neither c=n for the frozen family nor c=n+1 for d0=2736 is established for all n by this work.

## Checks and files

[Proved, exact finite checks] [checks.py](checks.py) uses rational arithmetic to check the endpoint construction at n=1,2,3,4,8,16,23,32,64,128. The largest check has 256 endpoints and 32896 unordered pair sums. Every margin passes; the smallest computed rational pair gap is 240/7 times the conservative bound. The Gram correlation identity is checked through j=63, and the scale-independent Schur correction at four scale ratios. The all-n claims have analytic proofs and do not rely on extrapolating these checks.

[Numerical, n=2] Ten 8-dimensional Lyapunov solves check the large-mass Gram asymptotics for the frozen and mobile stubs. Both approach the proved squared correlation 2/5. No microscopic many-body diagonalization, central-charge estimate, or all-n loop certificate was performed. [checks_results.json](checks_results.json) retains the exact ratios and numerical diagnostics. Run `python3 -B checks.py` from this folder.

[Proved, handoff] [PROOFS.md](PROOFS.md) contains all derivations, [PROGRESS.md](PROGRESS.md) records the work and neighboring status, and [STATE.md](STATE.md) is the local handoff. The unresolved step is a quantitative noncancellation theorem for the actual loop-dependent pencil, or a fixed-alphabet implementation of an explicit endpoint target. A stronger bandwidth hierarchy alone does not supply either.
