# The velocity budget for cell constructions at all orders: prove it or break it

You are working alone on one research question. Read only the files in this directory and write only inside it. Do not run network commands. There is no time limit. This is a long task. Do not stop after a first partial result. If the full statement is out of reach, push each remaining hypothesis as far as it goes and record exactly what is left.

## Background in one paragraph

For a nearest-neighbour chain H = sum_x h_{x,x+1} with local dimension d, whose infrared is a set of CFT sectors with central charges c_a and velocities v_a, the Casimir coefficient is A = sum_a c_a v_a. The central charge itself is unbounded at fixed d. Known constructions reach c = n at d = 9 by building cells of about n sites, each with an exactly degenerate ground space, coupled through one bond per interface at second order in perturbation theory, so that n Luttinger sectors appear with velocities spread geometrically. The conjectured replacement for c <= d - 1 is A <= C(d) M, with M the gauge-optimized half-spread of the bond (see `_bound_background.md` and `_Cd_target.md`). A previous run (`Cd62_PROOF.md`, second half, "A velocity budget for quadratic exchange through isolated cells") proved that bound for the second-order effective model. Your task is to remove the restriction to second order.

## The class

A chain in the class, in the sector of interest, is H = sum_cells H_c + sum_i V_i.

- Each H_c acts on a cell of n consecutive sites with nearest-neighbour terms of norm at most J_*, has an exactly degenerate ground space of some dimension D (possibly growing with n), and has gap Delta above it.
- V_i is a single bond across interface i, with norm g < Delta.
- Cells may be defined by conserved quantities (for example frozen kink positions, as in the d = 9 chain) rather than by an explicit decoupling of the Hamiltonian. Say how such sectors fit.

Everything below must be uniform in n, in D (equivalently in the number of encoded modes), and in the ring length L.

## What counts as an answer

(P) A theorem that every chain in the class with g/Delta below an explicit threshold has A <= C(d) J_*, or equivalently a velocity-resolved count N(u) = sum over sectors with v_a >= u of c_a <= B(d) (J_*/u)^p with p < 1, which gives A <= B(d)(4e)^(1-p) J_*/(1-p) (REPORT §9.93). The derivation of the low-energy theory may be at the physics standard, meaning Schrieffer–Wolff to all orders, perturbative renormalization group and bosonization. State which steps are rigorous and which are physics-level. Rigorous lemmas are welcome wherever they can be had, for example on the structure of the order-k effective Hamiltonian.

(D) A family in the class at fixed d along which A/J_* is unbounded, at the physics standard of the known constructions. Natural places to look are the window Delta/d^2 <~ g <~ Delta, intracell terms that couple the two ends of a cell, and interacting payloads.

## Points the answer has to address

1. Order-k structure. The order-k effective Hamiltonian on the cell ground spaces has on-cell terms, one-interface terms and terms spanning several interfaces. Bound each. The heuristic count says an order-k one-interface term has operator rank at most about d^(2k) across the interface and size about Delta (g/Delta)^k, which gives N(u) <~ (J_*/u)^p with p near 2 log d / log(Delta/g), summable when g < Delta/d^2. Either make this precise or show the true count is smaller (for example subexponential in k), which would remove the threshold. The second-order argument already replaces the energy denominator by a short sum of products, which is the mechanism to extend.

2. Long cells. The second-order proof needs neighbouring interface terms to nearly commute, which a Lieb–Robinson estimate inside gapped cells provides once n is much larger than J_*/Delta. At higher order, intracell terms that couple the two ends of a cell and terms spanning a whole cell appear. Show they are suppressed in n, or find where they are not.

3. Interactions. Beyond second order the effective theory is not quadratic. The infrared is then an interacting multicomponent Luttinger liquid (or partly gapped). Interactions of size epsilon per cell can raise slow velocities substantially (for one band, from v_F to about the square root of v_F epsilon), and physical speeds carry a factor n relative to cell units. Bound the sum of renormalized velocities, or show how interactions could make it large.

4. From the effective model to the chain. The previous run noted that a residual of norm epsilon per cell changes the scaled doubling defect by up to 4 epsilon L^2 / n when energies are compared in operator norm. At the physics standard one bounds the velocities of the infrared sectors instead. Say clearly which route you use.

5. Test case. Check whether the d = 9 frozen chain (`d9_construction_notes.tex`, `d9_README.md`, `d9_bond9.py`) belongs to the class, identify its n, D, Delta, g and J_*, and compare its velocity ladder (v_a proportional to r^a, r near 10^-2.74 at n = 9) with your bound.

## Files

- `Cd62_PROOF.md`, `Cd62_RESULT.md`: the previous run. Its second-order cell theorem is your starting point. Its first half (an XXZ strip estimate) is not needed.
- `REPORT_excerpts.md`: REPORT sections 9.93 (the fast-tower analysis and the p < 1 criterion), 9.98 (a physics-level reduction, emergent velocities and the second-order operator-rank budget), 9.207 (recent theorems, including a rank lemma for superlattices) and 9.210 (summary of the previous run).
- `_bound_background.md`, `_Cd_target.md`: conventions and the routes already ruled out for the general problem.
- `d9_construction_notes.tex`, `d9_README.md`, `d9_bond9.py`: the d = 9 construction. Paths mentioned inside these files are not available to you.

## Deliverables

- `PROGRESS.md`, kept current while you work.
- `RESULT.md`. Give the outcome first, (P), (D) or neither. Then list each claim tagged [proved], [physics-level], [numerical], [heuristic] or [assumed].
- `PROOF.md`, the complete argument or construction, with every step checkable by hand.
- Python check scripts (python3, standard library or numpy/scipy) for every numerical claim, with one line in `RESULT.md` saying how to run each.

No probabilities, odds or percent-confidence anywhere. Plain prose, no em-dashes.
