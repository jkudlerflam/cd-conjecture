# Handoff: cylinder c-theorem and the velocity-weighted bound (chat of 2026-11-07)

This file transfers a whole working session to a new chat. It records what was asked, what was established (with tags), what was refuted, the open questions, and the files. Tags follow the repository convention: proved, numerical, heuristic, assumed. Nothing here has been refereed by a human.

## 0. How the user wants to work

- Small steps. One conceptual move per response, 2 to 5 short paragraphs, then stop so the user can question or redirect. Do not jump over nontrivial steps. When there are several approaches, list them briefly and take them one at a time.
- Rigor, every step shown, and a consistency check whenever possible.
- Writing style: write like the user, not like an AI. Commas instead of dashes, few colons or semicolons, no italics for emphasis.
- Large computations go to the Typhon cluster (IAS, SLURM). Small ones can run locally.
- Only do computations that are on the critical path to a proof. The user explicitly declined exploratory tests that are not.
- Astra (gpt-6-astra) runs are launched by the user with `~/.claude/bin/gpt_solve.sh --model gpt-6-astra --effort high` on their machine. A cloud session cannot run Astra, so write a spec in `specs/` and push it. The user said not to wait for Astra results before continuing.

## 1. The original question (cylinder c-theorem)

Can the weak 2D c-theorem c_UV >= c_IR (unitary, Lorentz-invariant QFT flowing between CFTs) be proved from properties of the torus partition function Z(beta, L) or the finite-size ground energy E_0(R), without a monotone c-function at intermediate scales? The user's seven items: (1) is there a known proof purely from Z; (2) if so derive it and name the inequality; (3) if not, show why elementary properties (for example d^2 log Z/d beta^2 = Var H >= 0) fail, ideally with a counterexample having c_IR > c_UV; (4) keep monotonicity of c_eff(R) separate from the endpoint inequality; (5) whether the full finite torus with Z(beta, L) = Z(L, beta) and positivity in both channels gives more than L -> infinity; (6) the minimal extra input (is Z(beta, L) enough, or are metric or Weyl sources needed); (7) compare with entropic, thermal, spectral, modular and finite-size proofs.

### 1.1 Formulation (step 1)

- Z(beta, L) = Z~(M beta, M L), fixed up to a counterterm exp(-Lambda beta L). At fixed L it is a Laplace transform of a positive measure, so it is completely monotone in beta and determines the multiset {E_n(L)} without momentum labels. The two-variable function is the family of finite-volume spectra tied by Z~(x, y) = Z~(y, x). Var H >= 0 is implied by complete monotonicity (Cauchy-Schwarz), so it is not an independent input. [proved]
- C(R) = -(6R/pi)[E_0(R) - e_inf R], c_UV = C(0+), c_IR = C(inf). Strictly C measures c_eff = c - 12(h_min + hbar_min); unitarity in the untwisted sector gives c_eff = c. [standard]
- The two channels do different jobs. With R as the circle, E_0(R) is a ground energy. With R = 1/T and space L -> infinity, E_0(R)/R = f(T) is a free-energy density. Thermal-channel inequalities: s >= 0 is 2C - R C' >= 0 (C/R^2 nonincreasing); C_V >= 0 is (C/R)'' >= 0. [proved, given the swap]
- Class K: functions positive in both channels, symmetric, with the CFT asymptotics. A proof from Z alone means C(0+) >= C(inf) on all of K. Sufficient forms: C' <= 0; a sum rule with a positive integrand derived from double positivity; any functional inequality on K reducing to c_UV - c_IR. A negative answer is an element of K with C(0+) < C(inf).

### 1.2 Check: one channel is not enough [numerical, and analytic]

Free Majorana on a circle with periodic fermions (Ramond): E_0^P - e_inf R = +(m/pi) sum_k K_1(kmR)/k, exactly minus the periodic boson. C rises from -1 (h = hbar = 1/16) to 0. The ground-state channel Tr_R exp(-LH) is positive, but the crossed channel is Tr_NS (-1)^F exp(-RH), and indeed s < 0 and E_0 is convex there. So positivity of the ground-state channel plus existence of both limits does not imply C(0+) >= C(inf). Boson and antiperiodic Majorana satisfy all three conditions below; Ramond violates all three (`ceff_free_np.py`, `tilt_check.py`).

### 1.3 The real-tilt argument from a parallel session [checked here]

Claim from another chat: with P_+- = H +- P and F(b+, b-) = lim (1/L) log Tr exp(-b+ P+ - b- P-), boost covariance gives F = (b+ + b-) p(r), r^2 = 4 b+ b-, p(r) = -E_0(r)/r; convexity in b+ at fixed b- -> 0 gives Theta(q) = -L E_0(L) convex in q = L^2, hence c(L) nonincreasing and c_UV >= c_IR.

Checked in this chat:
- Hessian of F = beta^0 p(|beta|) at the point boosted by rapidity eta: T^2 [[cosh(eta) C_V, sinh(eta) s], [sinh(eta) s, cosh(eta) s]]. Positive iff C_V >= 0, s >= 0, C_V >= tanh^2(eta) s. At eta = 0 these are exactly the Euclidean thermal inequalities of 1.1. The boosted points add C_V >= s (sound speed <= 1). [proved, given F = beta^0 p]
- Theta''(q) = (C_V - s)/(4 L^3) and (L^2 d^2/dL^2 + L d/dL - 1) E_0 <= 0 confirmed. The free-gas formula in the other chat's note is only the k = 1 (Boltzmann) term; the full Bose sum is sum_k k m^4 K_1(kmL)/(4 pi m L), still positive. [proved]
- Convex in q plus a finite limit at infinity forces nonincreasing, so this route cannot separate item 4: the endpoint inequality comes with monotonicity.

### 1.4 F = b p(r) from the twisted torus, and the later proposition [checked here]

- Euclidean rotation of the torus with periods L and theta + i b gives (1/L) log Tr_L exp(-bH + i theta P) -> b p(sqrt(b^2 + theta^2)); continuing theta = i a gives b p(sqrt(b^2 - a^2)). [heuristic at first]
- Later update from the other chat, now checked: at fixed (L, r) write b = r cos(phi), theta = r sin(phi). Both sides of the real-angle rotation identity are analytic on the strip |Re phi| < pi/2, because with E - E_0 >= |p| each term is dominated by a heat trace at inverse temperature r cos(Re phi) exp(-|Im phi|). At phi = i eta this gives the exact finite-L duality Tr_L exp(-bH - aP) = Tr_r exp(-(L/r)(bH - aP)), r^2 = b^2 - a^2. The spectrum condition is needed on both circles. [proved, given the inputs]
- Independent sign check: holomorphic c = 24 theory with Z = J(tau). Left side J(i(b+a)/L); with -aP the right side is J(iL(b-a)/r^2) = J(-1/tau). With +aP it would be J(iL/(b-a)), which fails. Its ground state has P_0 = -2 pi/L, so this also exercises the P_0 version. [proved]
- Inputs of the continuum theorem: finite heat traces; Euclidean rotation invariance of the twisted torus; H - E_0 >= |P - P_0| on every circle; unique gapped ground state on every circle; endpoint CFTs with c_eff = c. Open gap (other chat): the spectrum condition, proved there only for states made by rotatable local insertions; closing it needs density of those states in every sector (a Reeh-Schlieder-type statement using only H >= 0 and local rotations; not checked).
- Answer to item 6 so far: beyond the rectangular Z(beta, L), the extra content is the spectrum condition, a statement about the joint (E, P) spectrum that the momentum-summed rectangular function cannot see. Local Weyl sources are not needed.

### 1.5 Not done on the original question

- Route 1: construct C(R) with C/R convex, C/R^2 nonincreasing, the correct asymptotics and C(inf) > C(0), to show the Euclidean thermal inequalities alone fail (the other chat also has no explicit counterexample). Route 2: whether finite-L double positivity of the rectangular torus rescues it.
- Literature cited by the other chat and not verified here: Delacretaz, Fitzpatrick, Katz, Walters arXiv:2105.02229 (Appendix B, c_s <= 1 implies s/T monotone); Gavassino et al. arXiv:2105.14621 (covariant stability implies c_s <= 1). Item 7 (comparison with Casini-Huerta, Cappelli-Friedan-Latorre and others) was not started.
- The other chat's files are outside this repository (Projects/Modular c-Theorem/, notes/cylinder_c_theorem/note/main.tex, referee_astra_B_2026-11-06.md, STATE.md).

## 2. Transfer to the c-d question (velocity-weighted bound A = sum_a c_a v_a <= f(d) M)

Target as relaxed by the user: any upper bound A <= f(d) M, any function f. Background in `../../STATUS.md` item 5 and `../../REPORT.md` sections 9.68 to 9.86.

### 2.1 The energy-current tilt [proved where stated]

- The lattice analog of P is the energy current J = sum_x i[h_{x-1,x}, h_{x,x+1}] (local, extensive; lattice momentum is not extensive). G(b, a) = lim (1/L) log Tr exp(-bH - aJ) is jointly convex even though [H, J] != 0 (log Tr e^X is convex in X). [proved]
- Low-temperature form: G = -b e + (pi/6) B T + (chi_0 T/2 + (pi/6) A T^3) a^2 + ..., B = sum c/v, with chi_0 the ground-state polarizability of J (zero when J is conserved). The a^2 coefficient carries A, as in C100. XX check: Var(J)/L over (pi/3) T^3 = 0.9998, 0.9993, 0.9972, 0.982 at T = 0.005, 0.01, 0.02, 0.05, deviation scaling like T^2 (`xx_tilt.py`). [numerical]

### 2.2 Convexity plus lattice bounds cannot bound A [proved, explicit model]

Model: G = -b e + sum_s (pi c_s T_s/(12 v_s)) [psi(T_s(b + a v_s)) + psi(T_s(b - a v_s))], psi(y) = (sqrt(y^2 + 4) - y)/2. It is convex, analytic and even in a; reduces to the multi-velocity CFT form when T_s(b - |a| v_s) >> 1; and all lattice constraints (G(0,0) = log d, |dG/db| <= ||h||, |dG/da| <= ||j_x||, tilted entropy in [0, log d]) reduce to three budgets with Dk_s = T_s/v_s: (pi/6) sum c_s Dk_s = log d; (pi/6) sum c_s v_s Dk_s^2 <= ||h|| - |e|; (pi/6) sum c_s v_s^2 Dk_s^2 <= ||j_x||. One species with fixed v (even v = v_LR) and Dk = 6 log d/(pi c) satisfies all of them with A = c v -> infinity.

### 2.3 Velocity-window inequality: proposed, then refuted [proved]

Proposed: Dk_s >= v_s/(kappa M) (equivalently T_s >= v_s^2/(kappa M)), which with the entropy budget would give A <= (6 kappa/pi) M log d. Refuted already for quasi-free chains at d = 4: a band -2t cos k crossing a flat level with on-site hybridization Delta has, with the Fermi level at the kink, one Dirac fermion with v = t and crossover of order Delta (curvature -v_c^2/(4 Delta)). Harmless for A because quasi-free chains have c <= m = log_2 d (det(h(k) - mu) has at most 2m roots) and A <= 2M (section 9.71). Bernstein controls single bands only.

### 2.4 Exact reformulation of "any f(d)" [proved]

A/M = integral_0^{u_LR} N(u) du with N(u) = sum_{v_s >= u M} c_s, u_LR = v_LR/M a universal constant (v_s <= 4eM, C109). With N_d(u) the sup over chains of local dimension d: necessary, N_d(u) < infinity for every u > 0; sufficient, integral of N_d finite. Thermal and Casimir data see only moments sum c v^k (k = -1 thermal, k = 1 Casimir and m_1, k = 3 next order), each valid only below each species' crossover, which can be pushed to zero for a fast species at fixed d (2.3).

### 2.5 Lattice monotonicity of c_A [numerical]

On the exact Bethe energies of the SU(N) Sym^2 chains (N = 2, 3, 4, 5, 7, 8; d up to 36; L = N to 120 in steps of N), c_A(L) = (6L/pi)(e_0 L - E_0(L)) is nonincreasing and -L E_0 is discretely convex in L^2 at every length, including L = N where E_0 = 0 (frustration free). c_A(N)/||h|| against A/||h||: 5.093 vs 3.142 (N = 2), 6.251 vs 5.236 (N = 4), 7.125 vs 6.597 (N = 8) (`lattice_cA.py`). An earlier objection (monotone from L = 2 would give A <= (48/pi)||h||, contradicted by the symmetric rank-N chains from N = 17) stands, but those chains are only vacuum-compatible at L = 0 mod N.

### 2.6 The ratio bound is the c-theorem in the continuum [proved in the continuum; lattice parts to verify]

- Continuum: by the swap, Z(beta, L) = Tr T^L with T = exp(-H_beta) >= 0, so R_L = log[Z(L)^2/Z(2L)] = S_2(T^L/Tr T^L) exactly, a crossed-channel Renyi-2 entropy. For beta >> L, R_L -> beta(E_{2L} - 2E_L) = (pi beta/(12L))(4c(L) - c(2L)) <= pi c_UV beta/(3L) by the c-theorem. The repo's classical proof R_L <= log d (section 9.80) is the same structure with a d-dimensional nonnegative transfer matrix.
- Lattice: Z(beta, L) = Tr T_beta^L with the quantum transfer matrix (Suzuki, Klumper), independent of L. For chains reflection positive across bond-centred planes (Dyson-Lieb-Simon class; believed to include XX and the Heisenberg antiferromagnet, not checked), T_beta can be chosen positive on even rings [assumed, to verify]. Then the Renyi identity is exact, and Renyi entropies of Gibbs states do not increase with inverse temperature (dS_2/dL = 2(<eps>_{2L} - <eps>_L) <= 0), so R_L is nonincreasing in L at every beta and D(L) = E_{2L} - 2E_L is nonincreasing in L [proved, given positivity of T_beta]. This does not bound A.
- Crossed-channel inputs on the lattice: T_beta >= 0; rotations of the Trotter circle are unitary and commute with T_beta, so Tr T^L exp(-a P~) is a positive log-convex trace; a crossed spectrum condition eps - eps_0 >= |p~|/v_LR from Lieb-Robinson via Matsubara decay (Hastings-Koma style) [heuristic, not checked in the literature]. Missing: rotation invariance mixing the Trotter circle and the lattice direction.

### 2.7 Velocity-mismatch obstruction [proved, under the stated generous assumption]

A species (c_s, v_s) becomes a crossed CFT with velocity 1/v_s; under the crossed tilt its chiral inverse temperatures are L/v_s +- a, null edge at a = L/v_s. The Lieb-Robinson condition guarantees only a < L/v_LR, so the reachable rapidity has tanh(eta_s) <= kappa_s = v_s/v_LR < 1. Even granting boost covariance up to there, the Hessian gives only C~_V >= kappa^2 s~, which integrates to s~(T~)/T~ <= [s~(T~_0)/T~_0](T~_0/T~)^{1 - kappa^2}, growing like L^{1 - kappa^2}: no bound unless kappa = 1. Heisenberg has kappa <= (pi/2)/(2e), about 0.29. Using the true scale-dependent maximal group velocity does not help: in the d = 4 avoided crossing the IR species moves at v_c/2 while its parent band moves at v_c above the arbitrarily small scale Delta. Physical reading: the continuum proof works because Lorentz invariance makes the IR light cone equal the UV light cone (kappa = 1); the c-d counterexamples live exactly in the gap kappa < 1.

### 2.8 Lemma VW: the velocity-weighted c-theorem must hold at every length [numerical and proved where stated]

- The natural shrinking-cone analog is monotone c_A(L) along vacuum-compatible lengths. Avoided-crossing check (Delta = 1e-3, crossover t/Delta = 1000): c_A/t = 2.036, 1.996, 1.861, 1.055, 1.000 at L = 8, 64, 512, 4096, 16384; c stays 1, the cone halves, c_A halves monotonically (`avoided_cA.py`). [numerical]
- Monotonicity alone cannot give f(d): it gives A <= c_A(L_0), useful only if the first vacuum-compatible length L_0 <= g(d). In route A a vacuum-compatible ring must hold a whole cell, which must hold a 2^n-dimensional code, so the cell length b >= n/log_2 d_0 and L_0 -> infinity at fixed d_0 (frozen and mobile cells alike). This refutes part (b) of the earlier program "(a) monotone c_A plus (b) L_0 <= g(d)". [proved]
- So the lemma must be an inequality at every vacuum-compatible L with constants depending only on d and M. In the continuum: L(E_{2L} - 2E_L) = (pi/12)(4c(L) - c(2L)) <= (pi/3) c_UV for every L. On the lattice this is Lemma VW: L(E_{2L} - 2E_L) <= C(d) M, with finite-temperature form R_L(beta) <= K(d) + kappa(d) M beta/L. That is exactly the repository's reconnection/ratio bound (sections 9.80, 9.86). Monotonicity is neither needed nor sufficient.
- What the continuum adds: for reflection-positive chains, e^{R_L} is the effective dimension of the positive spatial transfer operator T_beta^L, so Lemma VW says one lattice step in x carries at most K(d) plus kappa M/L per unit imaginary time at crossed temperature 1/L, at every L. The classical proof is the same statement with transfer dimension exactly d. This is where log d must enter: the capacity of a single site's imaginary-time history measured at crossed temperature 1/L.

## 3. Where the session stopped

Proposed next step (not started): bound the crossed effective dimension e^{R_L} = (Tr T_beta^L)^2/Tr T_beta^{2L} for reflection-positive chains directly from the structure of the transfer-matrix column (one site's imaginary-time history), uniformly in L. Do not repeat the routes listed as obstructed in `../../specs/_reconnection_lemma.md` and REPORT sections 9.85 and 9.86 (operator-order comparison, Araki constants, seam pinching with G <= 0, dimension-bounded divergences, local energy functionals, gradient norms, fixed-support bootstraps).

## 4. Artifacts

- `../../specs/C125_crossed_c_theorem.md`: Astra spec (crossed-channel c-theorem for the quantum transfer matrix), not yet run. It deliberately omits the conclusion of 2.7 so Astra's answer is independent.
- Draft PR https://github.com/jkudlerflam/cd-conjecture/pull/6 on branch `claude/loving-franklin-kj8ay2` (the spec and this folder).
- Scripts in this folder (numpy and scipy, each under a minute): `ceff_free_np.py` (C(R) for periodic boson, NS and Ramond Majorana), `tilt_check.py` (the three conditions 2C - RC' >= 0, (C/R)'' >= 0, C convex in R^2 on the same fields), `xx_tilt.py` (energy-current variance on XX), `lattice_cA.py` (monotonicity and L^2 convexity on the Bethe data; default path `../casimir_survey/astra_bethe/bethe_exact_energies.json`), `avoided_cA.py` (c_A(L) for the d = 4 avoided crossing).
