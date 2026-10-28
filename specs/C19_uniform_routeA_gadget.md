# C19: A uniform route-A cell, so unbounded c needs no universal machine

You are a long-running research agent. Your working directory is your own folder. If `PROGRESS.md` exists there, read it first and continue from where it stops. Update `PROGRESS.md` after every substantive step, and write `RESULT.md` at the end.

## Background

The repository is two levels above your folder (`../../`). Read `../../STATUS.md` first (the c–d conjecture and the cell construction), then `../../contrib/minimal_d0/README.md` in full (this run attacks its open direction (ii)), `../../contrib/unbounded_c/G1/G1_result.md` (the second-order exchange cell), `../work_C9/RESULT.md` (any logical bond through two qubit ports at second order, cell dimension 4⌈(6N³ + 13N)/4⌉), `../../contrib/programmable_cell/README.md`, `../../contrib/ir_control/README.md` (the stable target: n logical fermion flavours with distinct, incommensurate Fermi momenta and an on-site Z₂ acting as logical fermion parity, c = n), and `../work_C13/RESULT.md` (mobile cells). Astra C8 (`../work_C8/`) and C18 (`../work_C18/`) work on related parts of the computation track; read their PROGRESS.md and do not duplicate them.

## The problem

Route A gives unbounded c by second-order exchange through small (qubit) ports. As written, its cell gadget depends on the code size N without structure, so the cell program needs a universal history machine (Bausch–Cubitt–Ozols class, about 42 letters), giving d0 ≈ 130 with C7's walker tensored on. minimal_d0 speculates that a uniform gadget, whose program depends on n only through loop lengths, would remove the universal machine and bring unbounded c down to alphabets of the size of non-universal history constructions, perhaps d0 ≈ 20 to 30.

## Tasks

1. Choose the logical target to make the program uniform: the n-flavour incommensurate fermion chain with c = n is the right one (not SU(N)₁, which is unstable). Find a second-order cell whose code holds n flavours as n logical qubits or modes laid out in space, and whose exchange through the ports produces, at second order, independent hopping and chemical potential for each flavour with distinct values, using a program that is the same loop body repeated n times (or repeated a number of times derived from the cell length).
2. Show the second-order effective bond is the intended logical bond up to terms that only renormalize couplings, and that the on-site Z₂ (logical fermion parity) is exact. Compute it for small n by exact diagonalization of two and three cells.
3. Write the cell program as translation-invariant nearest-neighbour transition rules with the port idle, count the alphabet, and give d0 (with C7's walker tensored on, and with the merged selector if C18 provides one).
4. If the uniform gadget is impossible, prove the obstruction precisely, and fall back to writing out the Bausch–Cubitt–Ozols adaptation to marker-bounded cells with an idle port (the few extra letters minimal_d0 assumes), with the exact count.

## Rules

The standard is theoretical physics: standard RG, controlled perturbation theory and calibrated numerics; prove the finite-dimensional statements. Tag every claim as proved, numerical (with sizes), heuristic, or assumed. A precisely located obstruction counts as progress. Do not put probabilities, odds or percent-confidence estimates in any file. Python (numpy and scipy) only, inside your folder; no network; never modify files outside your folder; keep each computation to minutes. Write prose plainly, without em-dashes. Deliverables: `PROGRESS.md`, `RESULT.md` (construction, effective bond, alphabet and d0, checks with reproduction instructions, ranked next steps), code and data.
