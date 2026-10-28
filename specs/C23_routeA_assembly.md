# C23: Assemble route A end to end

Read `../../specs/_routeA_background.md` first. If `PROGRESS.md` exists in your folder, continue from it.

**Your part.** Be the integrator. Write the complete route-A chain as one document, with every interface between the pieces made explicit, and find every place where the pieces do not yet fit.

Tasks.
1. The parts list. For the target with n flavours, list each component (selector or clock-selected length, mobile wrapper and equalizer, history clock, cell core and port interface, the logical target's couplings), its alphabet, its bond terms, and the source file for each claim with its status tag. Mark explicitly which interfaces between components have been checked and which are assumed.
2. The scale hierarchy. Write every energy scale (selector margin, clock gap, cell gap, port visibility, second-order logical bandwidth, marker hopping and length-fluid bandwidth, payload incommensurability margins of order n^{-4} from C20, third-order corrections) as explicit functions of n and the cell length, and check that a consistent ordering exists for every n. Give the numbers for n = 2 and n = 3.
3. Corrections. Compute or bound the leading corrections the pieces generate on each other: third-order non-Gaussian terms from C19 at the target's margins, the O(1/J_c) clock errors, the residual length–payload coupling after equalization, leakage. Decide at the physics standard which renormalize couplings only and which could destabilize the target.
4. The accumulation point. As n grows, identify the limiting normalized bond and the nature of its gapless spectrum, consistent with C13 §9.
5. Output: RESULT.md as a self-contained description of the route-A counterexample as it stands, with a list of the remaining gaps ranked by how much they matter, and for each a concrete next calculation. If you find an incompatibility that kills route A, say so plainly.
