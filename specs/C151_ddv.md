# C151: The ferromagnetic strip by the Destri–de Vega equation

Your folder is `../work_C151/`. Read `specs/_strip_context.md` first; it states the target and what is known.

## The idea

The Destri–de Vega (Klümper–Batchelor–Pearce) nonlinear integral equation is an exact finite-L reformulation of the ground-state Bethe equations. It expresses root sums as contour integrals on lines Im x = ±ε inside the analyticity strip of the kernel, with correction terms that are functionals of log(1 + e^(±iLZ(x ± iε))). These are exponentially small, of order e^(−2πLσε), wherever the counting density is not small. Near the ferromagnet the strip width is limited by η, so the decay rate is of order λ = ηL r, and the shell is precisely where these terms are of order one. The hope is that the exact representation shows the shell's contribution with a sign, or with a cancellation, that the per-cell quadrature bound cannot see.

## Tasks

1. Derive the DdV equation for the XXZ ring at flux Lπ, in every magnetization sector (finite Fermi boundary as well as half filling), for −1 < Δ ≤ 0, with exact statements of the analyticity strip, the contour choice, and the hypotheses (real ground-state roots, monotone counting function from C136). Check the literature (Destri–de Vega 1992 and 1995, Klümper–Batchelor–Pearce 1991, Kozlowski) for rigorous versions.
2. Write F_L(n) − L e_bulk exactly in terms of the nonlinear terms, and bound them uniformly in η, L and n, with explicit constants. Pay special attention to the region where LZ has small imaginary part (the shell and tails, λ of order 1 to log(1/η)), and to the Fermi edge at finite filling.
3. If a uniform sector bound is out of reach, bound the matched defect L·D directly in the DdV language, where the bulk terms cancel between 2L and L.
4. Assemble the strip theorem with explicit η₀ and C, or state exactly which estimate is left and its numerical margin. Check numerically against the Bethe solvers, down to 1 + Δ = 1e-6.
