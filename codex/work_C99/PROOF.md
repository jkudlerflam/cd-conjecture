# C99 proofs and protocol audit

2026-10-29. Natural logarithms, lattice spacing and hbar equal to one. Claims labeled [Proved] are deductions from the hypotheses stated with them. A continuum formula is not asserted to hold at a microscopic lattice scale. All writes are confined to this folder.

## 1. Setting and the exact unitary resource

[Assumed setting] The Hamiltonian is H = sum_j h_(j,j+1), with one-site dimension d and the optimized bond scale

    M = inf_(o=o*, b real) ||h + o tensor I - I tensor o - b I||.

The infrared sectors have A = sum_a c_a v_a and vacuum-compatible Casimir term -pi A/(6L). We retain C69's convention for nonchiral c_a and its boundary-state calculation. No CFT identification of a lattice quench is implicit below.

[Published theorem] For a Hamiltonian g acting on the two d-dimensional sites adjacent to a cut, arbitrary ancillary systems and arbitrary initial pure states are allowed in

    |dS_left/dt| <= 18 ||g|| log d.                         (1)

This is the ancillary SIE theorem, [Van Acoleyen, Marien and Verstraete, arXiv:1304.5931](https://arxiv.org/pdf/1304.5931), equations (2) and (15). The paper obtains the constant 18; its cited stronger small-d bounds cover d=2. The theorem is used, not reproved.

[Proved] Write the finite-chain Hamiltonian as H_left + H_right + h_cut. Unitaries generated within either half preserve the reduced-state eigenvalues, so only h_cut contributes to the instantaneous entropy derivative. Adding a scalar or either one-site term to h_cut also contributes zero. Applying (1) to arbitrarily near-minimizing gauge representatives gives

    |dS_left/dt| <= K,       K = 18 M log d.                (2)

For mixed states, fix a purification and regard its reference as part of the complementary ancilla. For infinite systems all applications are limits of finite-volume entropy increments, not subtractions of infinite entropies.

Define the ancillary entangling capacity

    E(h) = sup_(finite ancillas, pure psi) |dS_left/dt|_psi.

[Proved] E(h) is unchanged by adding arbitrary local terms a tensor I + I tensor b, by local unitary conjugations, or by the bond gauge above. Thus

    E(h) <= 18 M log d.                                    (3)

At each instant, arbitrary control unitaries within either half merely choose another state in this supremum. Consequently any deterministic protocol using h across this cut, with arbitrarily strong local controls and local ancillas, satisfies

    |S(T)-S(0)| <= E(h) T.                                 (4)

If the crossing interaction is switched or rescaled by a real u(t), T is replaced by int |u(t)| dt. Sequential uses add this cost. Initially shared entanglement is allowed, but the left side is the net increment. Discarding or measurements require a different entanglement measure or an average-branch accounting; postselecting an entropy increase is not covered by (4).

[Proved reduction, not a no-go theorem] A protocol with Delta S >= alpha A T and a uniform alpha>0 would prove E(h) >= alpha A, hence A <= (18/alpha) M log d. Conversely E(h) >= kappa A makes near-optimizing states witnesses of such instantaneous rates, up to any fixed loss in kappa. Thus the unrestricted single-cut, unitary protocol route reduces to a microscopic-to-infrared lower bound on E(h). SIE proves only the opposite, upper estimate (3). Optimizing over states does not provide this missing lower bound.

This reduction leaves a real research problem. Neither a counterexample to E(h) >= kappa A nor a proof of that inequality is supplied here.

## 2. Weighted rates, integration, and many scales

[Proved conditional theorem] Consider entropy functionals

    F_j(t) = sum_k q_jk S(R_jk,t),
    n_j = sum_k |q_jk| number_of_boundary_bonds(R_jk).

First collect identical entropy terms. The elementary SIE estimate gives |r_j(t)| <= K n_j, where r_j = dF_j/dt. Suppose a continuum protocol at scale ell_j has a calibrated response

    r_j(t) = (A/ell_j) f_j(t),
    |f_j(t)| <= b n_j,        ell_j >= ell_* > 0,           (5)

with b independent of the scales and the chain. The theorem applies to exactly this class; an anomalously large shape response outside (5) must be analyzed separately. Corrections to matching can only be used after giving an error bound.

Take any signed integrable weights w_j and define

    Q = sum_j int w_j(t) r_j(t) dt = A beta,
    beta = sum_j int w_j(t) f_j(t)/ell_j dt,
    D = sum_j n_j int |w_j(t)| dt < infinity.

The triangle inequality and (5) give, step by step,

    |Q| <= K D,
    |beta| <= sum_j int |w_j| |f_j|/ell_j
           <= (b/ell_*) D.                                (6)

If beta is nonzero, the upper estimate delivered by these inequalities is

    A <= K D/|beta|,
    D/|beta| >= ell_*/b.                                   (7)

The second line concerns the size of the derived upper bound, not a lower bound on A. No choice of signs, replication, integration, or scale-dependent weights can make the SIE upper estimate in this class better than K ell_*/b. Discrete samples follow by replacing integrals with sums. Infinite families follow by absolute convergence. A divergent D supplies no finite bound from (2).

[Proved application] The boundary-state plateau in C69 has b=pi/12 per cut. Therefore all combinations of these plateaus at ell_j>=ell_* retain a factor at least 12 ell_*/pi, before errors. C69's lower matching factor alpha yields

    A <= 216 ell_* M log d/(pi alpha),                     (8)

when matching at ell_* is available. Negative coefficients can exploit correlations among actual rates only with a new joint inequality; independent SIE envelopes contain no such information.

[Proved reparametrization check] Under t=t(s), the generator is (dt/ds)H and

    dS/ds = (dt/ds) dS/dt,
    |dS/ds| <= K |dt/ds|.

Using a dimensionless time or a state-dependent clock multiplies both the response and the bound. Multiplication by a length is the same operation at the level of the inequality.

[Proved spatial parallelization check] A repeated arrangement with one active cut per ell sites has O(1/ell) cuts per unit length. Both the total CFT rate per unit length and its SIE bound acquire this cut density. Dividing by the number of cuts returns (7). A parallel construction does not make many cuts into one interaction channel.

## 3. An arbitrary state's correlation length cannot repair SIE

[Proved counterexample to an auxiliary inequality] There is one fixed translation-invariant nearest-neighbor chain with d=4, a conventional critical sector and a gapped spectator, for which

    sup_psi xi_Z(psi) |dS_left/dt| = infinity.              (9)

Here xi_Z is a genuine exponential correlation length of a local observable in the spectator state. Thus even a state-defined correlation length does not give a general bound xi |dot S| <= C M log d.

Each site consists of an active qubit a and spectator qubit s. Fix J,Delta>0 and take

    H = (J/2) sum_j (X^a_j X^a_(j+1)+Y^a_j Y^a_(j+1))
        + (Delta/2) sum_j (1-Z^s_j).                       (10)

Distribute each spectator field equally over its two neighboring bonds. The active bond has spectrum {-J,0,0,J}; the spectator bond has spectrum {0,Delta/2,Delta/2,Delta}. They commute. Hence this representative has half-spread J+Delta/2. For the lower bound, take the two active Bell eigenstates with energies +/-J and put both spectators in |0>. Their left and right one-site marginals coincide, so the expectation of every coboundary vanishes. Their energy difference remains 2J in every gauge representative. Therefore

    J <= M <= J+Delta/2,

uniformly for all the states below. The free active chain maps to hopping with dispersion 2J cos k, so its two Fermi points have speed 2J and its usual nonchiral charge is c=1. The spectator ground state is unique and has gap Delta. It contributes no infrared CFT sector. One can also extract A directly: on compatible rings L=4m, the occupied antiperiodic fermion momenta give a geometric cosine sum and

    E_0(L) = -2J/sin(pi/L)
           = -(2J/pi)L - pi J/(3L) + O(L^(-3)).

Thus A=2J is fixed, not growing in this example. With d=4 and M>=J, it obeys A<=pi M log_2 d.

Across the chosen cut put the adjacent active qubits in

    |phi_p> = sqrt(p)|01> + i sqrt(1-p)|10>,   p=1/4,

and put the other active sites in product states. The crossing Hamiltonian on {|01>,|10>} is J sigma_x. At time zero,

    p' = 2J sqrt(p(1-p)),
    r_0 = d[-p log p-(1-p)log(1-p)]/dt
        = 2J sqrt(p(1-p)) log((1-p)/p)
        = (sqrt(3)/2) J log 3 > 0.                        (11)

All other active terms are within one half and give zero entropy derivative, irrespective of this state's subsequent evolution.

For the spectators choose the pure state with positive amplitudes

    |chi_k> = Z_k^(-1/2) sum_(s_j=+/-1)
              exp[(k/2) sum_j s_j s_(j+1)] |{s_j}>.        (12)

Use finite open chains, then take their thermodynamic limit at fixed k. Squared amplitudes are the zero-field classical Ising measure. Its transfer matrix has eigenvalues 2 cosh k and 2 sinh k. Thus, in the infinite bulk,

    <Z^s_0>=0,    <Z^s_0 Z^s_r>=(tanh k)^r,
    xi_Z(k) = -1/log(tanh k) -> infinity as k->infinity.   (13)

The amplitude factorization in (12) gives Schmidt rank at most two across one open-chain cut, so this is not an infinite-entropy construction. The spectator Hamiltonian in (10) is a sum of on-site terms. Its evolution is a product of on-site unitaries, and its contribution to the cut entropy derivative is zero. It also leaves the Z correlations fixed.

Take psi=phi_p tensor chi_k, with the specified other active sites. Entropy adds between species. Equations (11)-(13) give xi_Z(k) dot S=xi_Z(k)r_0 -> infinity, while H,d,M,A stay fixed. This proves (9).

[Scope] The initial states need not be translation invariant or low energy; the Hamiltonian is. This refutes a supremum over arbitrary states with a full-state correlation length. It does not refute a quantity restricted to specially prepared low-energy states, or a length tied specifically to the degrees of freedom producing the entropy. Such a restricted length needs its own theorem. The example satisfies the A conjecture and is not a mechanism for violating it.

## 4. Local joining quenches and a positive integral for A

[Continuum input] For a regularized conformal join with spatial smoothing epsilon, the half-line entropy increment in sector a is

    Delta S_a(t) = (c_a/6) log[1+(v_a t/epsilon)^2].        (14)

This is the regulated continuum formula of [Calabrese and Cardy, arXiv:0708.3750](https://arxiv.org/pdf/0708.3750), equation (14), with its value at zero subtracted and physical velocity restored. Its use for a lattice near the regulator scale is an additional matching assumption.

[Proved from (14)] Direct differentiation gives

    r_a(t) = c_a v_a^2 t/[3(epsilon^2+v_a^2 t^2)]
           = (c_a v_a/(3epsilon)) x_a/(1+x_a^2),
    x_a=v_a t/epsilon.

Since 2x<=1+x^2,

    0<=r(t)<=A/(6epsilon),
    r(t) ~ c/(3t) at late times.                          (15)

With a common speed, the upper value is reached at t=epsilon/v. With unequal speeds the peaks do not generally coincide, so that equality cannot be used to infer A from a single peak. The late-time coefficient counts c, losing the velocity weighting.

[Proved continuum identity] There is nevertheless a positive weighted integral that measures A for arbitrary positive sector speeds:

    int_0^infinity r_a(t) dt/t
       = (c_a v_a^2/3) int_0^infinity dt/(epsilon^2+v_a^2 t^2)
       = pi c_a v_a/(6epsilon),

    A = (6epsilon/pi) int_0^infinity r(t) dt/t.            (16)

At finite endpoints t_-<t_+ the contribution of sector a is

    (c_a v_a/(3epsilon))
       [arctan(v_a t_+/epsilon)-arctan(v_a t_-/epsilon)].  (17)

This is a usable continuum observable family, but it does not solve the lattice estimate. SIE gives int K dt/t, divergent at both endpoints. On a finite window it gives

    (6epsilon/pi) int_(t_-)^(t_+) r(t)dt/t
        <= (6epsilon/pi) K log(t_+/t_-).                  (18)

To approximate (16) one also needs continuum control below the fastest peak and beyond the slowest peak. Short-time estimates specific to the lattice state might regularize the ultraviolet endpoint, but the large-time weighted bound and the epsilon prefactor remain uncontrolled by SIE. Using the CFT tail c/(3t) in the upper estimate reinserts the unknown central charge. Inferring a UV-sized epsilon from a one-bond join does not justify a uniform CFT approximation near t~epsilon/v.

[Scope of the literature input] The cited calculation relates epsilon to a nonuniversal boundary entropy constant through its equation (16). That relation supplies no uniform estimate over the present Hamiltonian class on the regulator or the onset of continuum scaling. It is not a theorem that all microscopic one-bond joins resolve all infrared sectors on one lattice spacing.

## 5. Inhomogeneous and Floquet evolution

[Proved] Suppose the crossing bond is g_cut(t). Linearity of the derivative and SIE give

    |Delta S| <= 18 log d int ||g_cut(t)||_local-quotient dt, (19)

where the norm can be minimized over a scalar and terms supported separately on the two sides. In particular a scalar envelope u(t)h_cut costs |u(t)|M. Inhomogeneity away from the cut and any within-half control leave this estimate valid. Repeating N periods multiplies the available crossing action by N; a rate per period is bounded by the average crossing norm.

[Proved obstruction to rescaling] Running a length-ell conformal process ell times faster by multiplying its physical generator by ell also multiplies its crossing norm by ell. If the norms remain bounded, the required physical time remains in (19). Sector-selective controls cannot be assigned a cost from the original M without a microscopic implementation.

[Assumed/open] A different bounded drive might realize a response >=kappa A of the original H. Neither SIE nor a Floquet effective Hamiltonian identity proves that lower response. Higher commutators can produce longer-range effective couplings; their support and norms must be charged to the actual drive in (19). A specially designed protocol remains possible if this missing comparison can be proved.

## 6. Distant mutual information and Lieb-Robinson bounds

[Proved SIE estimate] For separated intervals X,Y,

    I(X:Y)=S(X)+S(Y)-S(XY).

There are two, two, and four boundary bonds respectively. Hence

    |dI(X:Y)/dt| <= 8K = 144 M log d.                     (20)

This crude bound has no separation factor. For scale-covariant CFT interval geometries, derivatives at fixed scaled time carry 1/ell. Mutual information also depends on the full operator content, so its response need not be proportional to A without further assumptions.

[Proved diagnostic] A state-independent exponentially small MI-rate bound before the distance light cone is false. Let A be a qubit, E its neighbor and B a distant qubit. Initially put A in |0> and E,B in a Bell pair |Phi+>. Evolve by J X_A X_E. Then

    |psi(t)> = cos(Jt)|0>_A|Phi+>_EB
                -i sin(Jt)|1>_A X_E|Phi+>_EB.

The two Bell states are orthogonal, rho_B=rho_E=I/2, and the tripartite state is pure. Therefore

    I(A:B)=S(A)+S(B)-S(E)=H_2(sin^2 Jt).                  (21)

At Jt=pi/8 its derivative is positive and independent of the distance to B. If desired, this example uses the translation-invariant chain J sum X_j X_(j+1): put every intervening and exterior spin other than A,E,B in |+>. Those spins remain X eigenstates. All additional terms reduce to commuting local fields on A,E,B and do not change (21). Choose B far enough away to have |+> neighbors.

The correlation was already present between E and B. Local evolution transfers it to A without transmitting a new signal across the separation. This does not contradict LR. It is an exact lattice counterexample to the proposed all-state separation suppression, not to the Casimir conjecture.

The diagnostic can also be embedded in the target CFT class. Use the commuting spectator Hamiltonian -J sum X_j X_(j+1)-b sum X_j with b>0, which has a unique product ground state and a nonzero gap. The additional X fields are local, commute with the bond evolution, and leave (21) unchanged; reversing J also leaves sin^2(Jt) unchanged. Tensor this spectator with an independent critical XX chain in its ground state. The on-site dimension is four, the critical-sector MI is stationary, and the total MI derivative still has the distance-independent contribution (21).

[Conditional implication, inherited from C69] With dynamical visibility of each CFT sector in local observables, a conservative LR velocity yields v_a<=6eM and hence A<=6eM sum_a c_a. This does not bound the weighted number of channels. For initially uncorrelated preparations, useful LR bounds before the cone are available, but they do not control the amplitude of an inside-cone infrared entropy response. The ground state itself is stationary, with zero time derivative of every fixed-region MI.

## 7. Modular evolution costs a different generator

[Proved] For a finite regulated bipartition, the full modular unitary generated by K_X-K_Xc factorizes across X|Xc. Its action cannot change S(X) in any state. In particular the modular evolution of its own vacuum gives no entanglement-production rate across that cut. Moving the region or differentiating its modular Hamiltonian is a different operation.

[Continuum input] For a CFT interval (-R,R) of speed v, the vacuum modular generator is

    K_R = (2pi/v) int_(-R)^R [(R^2-x^2)/(2R)] H_density(x) dx
          + constant.                                   (22)

This follows from [Casini, Huerta and Myers, arXiv:1102.0440](https://arxiv.org/pdf/1102.0440), equation (2.20), with physical velocity restored.

[Proved accounting conditional on a local lattice implementation] A proposed lattice generator G=sum_j f_j h_j has crossing strength |f_cut|M in the instantaneous SIE bound. The weighted coboundary produces only extra on-site terms, since

    sum_j f_j(o_j-o_(j+1)) = sum_j (f_j-f_(j-1))o_j

up to endpoint terms. These terms do not contribute to the instantaneous cut entropy derivative. For a cut through the middle of (22), f_cut is of order R/v. Thus a response in modular time comes with an SIE bound of order (R/v)M log d. Dividing the generator by R/v divides both. At the interval boundary the generator is confined to one side and gives zero entropy production across that boundary.

The exact microscopic modular generator generally has no assumed local bounded decomposition resembling (22). Its spectrum -log rho has no uniform norm bound from d alone, because nonzero eigenvalues of rho can approach zero. Differences of modular generators, modular inclusions, and relative modular flows might have useful additional structure, but the physical M is not automatically their crossing norm. Establishing that structure and a response proportional to A would be a new theorem, not an application of (2).

## 8. A Euclidean ancillary SIE bound is false

[Proved counterexample] Let the left Hilbert space be A tensor a, with dim A=2, and let the right system have orthonormal vectors |0>,...,|N>. The right vectors and the left ancilla can be encoded in qubits. Consider

    |psi_p> = sqrt(p)|0>_A|0>_a|0>_R
       + sqrt((1-p)/N) sum_(j=1)^N |1>_A|j>_a|j>_R.      (23)

These are Schmidt terms. Let g=M Z_A act only on A, strictly within the left side, and evolve by the normalized imaginary-time map

    |psi(tau)> = exp(-tau g)|psi_p>/||exp(-tau g)|psi_p>||.

The total weight in the first Schmidt sector is

    p(tau) = p exp(-2M tau)/[p exp(-2M tau)+(1-p)exp(2M tau)],
    p'(tau) = -4M p(tau)(1-p(tau)).                       (24)

The cut entropy is

    S(tau) = H_2(p(tau))+(1-p(tau))log N.

At p=1/2,

    S'(0) = M log N.                                     (25)

The acted-on dimension is two and ||g||=M, while N is arbitrary. Consequently there is no constant C for a Euclidean analogue |S'|<=C||g||log d independent of ancillas. There is not even a boundary-only estimate, because g does not cross the cut. This is the normalized conditional filtering mechanism: a local filter reweights previously entangled Schmidt sectors of different multiplicity.

[Proved differential explanation] With rho the normalized pure-state density matrix,

    partial_tau rho = -{g,rho}+2<g>rho.

For full-rank reduced states, or on the fixed nonzero Schmidt support in (23),

    partial_tau S_left
       = <{g,log rho_left tensor I}> - 2<g><log rho_left tensor I>.

This is an anticommutator covariance. The real-time derivative is instead a commutator, for which terms supported on one side vanish. There is no such cancellation here. Analytic continuation of the real-time SIE inequality is invalid.

[Scope] The physical Casimir conjecture concerns unitary real-time Hamiltonians. Equation (25) does not violate it. Imaginary-time ground-state projection or replica identities may still give useful bounds with extra conditions on the state, correlations, or modular spectrum, but they cannot use an ancillary boundary-rate bound of the form (2). Integration cannot rescue an inequality already false pointwise, and a postselected physical realization needs branch-probability accounting.

## 9. Why none of these failures constructs a large-A chain

[Proved logical nonimplication] Let K=18 M log d. Fix a speed v_0 below a proposed LR cap and allow a formal arbitrary c, so A=c v_0. Choose a length Lambda large enough that

    Lambda >= pi A/(12K),
    Lambda >= pi c/(12 log d).

For every ell>=Lambda the formal boundary-state profiles

    Delta S_ell(t)=(c/6)log cosh(pi v_0 t/(2ell))

obey |dot S|<=K and the corresponding entropy density pi c/(12ell)<=log d. Their speed stays v_0. These profiles can satisfy all the independent rate and entropy-density inequalities for arbitrarily large A by delaying the regime of conformal validity. They are scalar response data, not a construction of local Hamiltonians, and impose none of the additional microscopic consistency conditions a counterexample needs. The weighted estimates of Section 2 do not change this observation.

[Heuristic candidate mechanism, not a counterexample] A violation would need many modes with an unsuppressed total sum c_a v_a, while the modes become visible only beyond a sufficiently long crossover scale. The given inequalities do not rule this out. They also do not realize it. The supplied background's fixed-alphabet constructions make additional modes slow and therefore do not provide this mechanism. A large crossover length by itself is not evidence that A/(M log d) grows.

[Open] The main conjecture A<=C M log d is unresolved here, as is the sharp C=pi/log 2. A universal impossibility theorem for all entanglement-rate arguments is also not proved. What is established is (i) a precise obstruction for weighted physical-rate bounds and all the protocols assessed above, (ii) exact counterexamples to three tempting auxiliary estimates, and (iii) the remaining microscopic comparison E(h)>=kappa A that would suffice to escape the obstruction.
