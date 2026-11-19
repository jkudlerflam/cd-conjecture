# Is the doubling constant C(d) finite? Prove it or disprove it.

You are working alone on one research question. Read only the files in this directory and write only inside it. Do not run network commands. There is no time limit. Work until you have a proof, a disproof, or a precise account of where each route you tried stops.

## Setting

A one-site translation-invariant nearest-neighbour chain H_L = sum_x h_{x,x+1} on a ring of L sites with local dimension d, h a Hermitian d^2 x d^2 matrix. M is the gauge-optimized half-spread of the bond,

    M = min over alpha, o of (lambda_max(h') - lambda_min(h'))/2,   h' = h + alpha 1 + o (x) 1 - 1 (x) o,

which leaves every ring Hamiltonian unchanged up to a constant. E_L is the ring ground energy. When the infrared is a finite set of CFT sectors with central charges c_a and velocities v_a, the Casimir coefficient is A = sum_a c_a v_a, with E_L = e_0 L - pi A/(6L) + o(1/L) and L (E_{2L} - 2 E_L) -> pi A/4.

The target is

    C(d) = sup over bonds h with M(h) > 0 and over vacuum-compatible L of L (E_{2L} - 2 E_L)/M.

## What counts as an answer

(P) A proof that C(d) is finite for every d, with any explicit function of d, however fast it grows. The two variants in `_Cd_target.md` also count, namely a Casimir lower bound E_L >= e_0 L - K(d) M/L on vacuum-compatible lengths, or the finite-temperature ratio bound R_L(beta) <= K(d) + kappa(d) M beta/L with R_L(beta) = log[Z(L,beta)^2/Z(2L,beta)]. The proof must be rigorous and must give "vacuum-compatible" a precise meaning.

(D) A disproof, namely a family of nearest-neighbour bonds at one fixed d along which A/M (or L(E_{2L} - 2E_L)/M on vacuum-compatible lengths) is unbounded. For a disproof the standard is that of the known c > d - 1 constructions, namely an explicit Hamiltonian, controlled perturbation theory and renormalization group, and numerics calibrated against exactly solvable references. A harmonic velocity ladder c_k = 1, v_k of order M/k for k = 1, ..., N already gives A of order M log N.

A partial result is worth reporting only if it lies on the path to (P) or (D). That means a lemma the proof needs, an obstruction showing a route cannot work, or a construction whose failure point is located precisely.

## Read first, in this order

1. `_bound_background.md`, the conventions and the early rounds.
2. `_reconnection_lemma.md`, the seam-pinching reduction.
3. `_Cd_target.md`, the target and the five routes that are already ruled out. Do not repeat them.
4. `REPORT_9.68-9.109.md`, every round on this bound so far.

Paths such as `../work_C121/` mentioned in those files are not available to you. The REPORT sections summarize what they contain.

In brief, the bound is proved for quasi-free chains (A <= 2M), for classical and commuting chains (C = 0), for independent XX stacks, for WZW sectors tied to on-site symmetries (conditionally on the WZW dictionary), and for the d = 2 U(1) class (XXZ plus DM along z plus field) with explicit constants away from a thin strip next to the ferromagnetic point. No chain checked so far violates A <= pi M log_2 d.

## The mechanism a disproof has to beat

Chains with c far above d - 1 at fixed d exist at the physics standard. In the smallest ones (d = 9 and 10), a soliton lattice of kinks forms cells of about n sites, each cell carries an exactly degenerate 2^n-dimensional ground space (flat bands), and neighbouring cells talk through a single bond by second-order superexchange, through virtual creation of a pair of fractionally charged kinks. After Jordan–Wigner the effective theory is exactly quadratic, H_eff = sum_m [c_m^dag (mu + f) c_m + c_m^dag T c_{m+1} + h.c.], with an n x n hopping matrix T. All entries of T are of the same perturbative order, yet its eigenvalues spread geometrically, so the n Luttinger sectors have velocities v_a proportional to r^a (r near 10^-2.74 in one d = 9 example). Thus c = n grows while A converges. In an older family with much larger d, all n payload velocities carry a common suppression factor D^(-2n-4).

The heuristic reason (REPORT §9.93, §9.98) is a bond bottleneck. Everything passing between two cells passes through one bond, whose operator space has dimension d^2. At second order at most d^2 - 1 independent inter-cell channels open, and further sectors appear only at higher order, hence exponentially slower. Nothing proves that this bottleneck is unavoidable, and nothing excludes a sector that emerges late but fast.

## Two observations not yet in the files

1. Partition functions alone cannot give the bound. The momentum-twisted traces Tr(T^j e^{-beta H_L}), with T the lattice translation, are traces of the spatial transfer matrix with a rational rotation of the Trotter circle. With g = gcd(j, L), q = L/g and b j = g mod L, the space-time torus spanned by (L, 0) and (j, beta) is also spanned by (g, b beta) and (0, q beta). So they lie inside the C125 no-go of §9.88. More generally, every infrared consistency condition (modular invariance at any genus, crossing, positivity) is unchanged when each sector's velocity and its moduli are rescaled together, so none of them can produce the bond scale M on the right-hand side. M has to enter through an operator inequality for h at the scale of the ring. The dimension d enters partition functions only through state counting (log d per site), which bounds sum_a c_a k_a* <= (3/pi) log d, with k_a* the momentum window over which sector a is linear, and not A.

2. Every known obstruction is a long-crossover effect (lesson 3 of `_Cd_target.md`), meaning sectors whose conformal regime begins only at lengths much larger than one lattice spacing. A proof has to bound the doubling defect at length L using only data at length L.

## Directions you may take (suggestions only)

- For (P), make the bond bottleneck rigorous. One form is a velocity-resolved count N(u) = sum over sectors with v_a >= u of c_a <= B(d) (M/u)^p with p < 1, which gives A <= B(d) (4e)^(1-p) M/(1-p) (§9.93). A finite-size form that avoids the infrared assumption would be better. Keep in mind that unfiltered counting of low ring levels fails (C115).
- For (P), compare scales directly, bounding D(2L) in terms of D(L), with D(L) = E_{2L} - 2E_L, up to errors controlled by data at length L, iterated over dyadic lengths. Exact convexity of -L E_L in L^2 fails near the ferromagnetic point of XXZ (lesson 5), so any such comparison must survive a z = 2 to z = 1 crossover.
- For (D), look for a fast tower, meaning sectors produced at bounded perturbative order whose velocities decay only polynomially. Examples are payload carried across cells by first-order resonant transport, or a hierarchy of cell sizes. Then decide whether the bond bottleneck forces geometric decay.

## Deliverables

- `PROGRESS.md`, kept current while you work, so that a reader can follow the work if it is interrupted.
- `RESULT.md`. Give the outcome first, (P), (D) or neither. Then list each claim tagged [proved], [numerical], [heuristic] or [assumed].
- `PROOF.md`, the complete proof or construction, with every step checkable by hand.
- Python check scripts (python3, numpy, scipy) in this folder for every numerical claim, with one line in `RESULT.md` saying how to run each.

If the target is false or ill-posed as stated, for example through the meaning of "vacuum-compatible" or the choice of M, say so and give the corrected statement. No probabilities, odds or percent-confidence anywhere. Plain prose, no em-dashes.
