# C21: Compile C19's uniform cell into a constant-alphabet nearest-neighbour chain

Read `../../specs/_routeA_background.md` first. If `PROGRESS.md` exists in your folder, continue from it.

**Your part.** Turn C19's uniform second-order cell into local, translation-invariant, nearest-neighbour terms on a fixed alphabet, overcoming the frustration-free interface obstruction (bottleneck 2).

Tasks.
1. Design the port interface. Since positive checks that annihilate the code give scalar resolvents, find an interface that produces logical-operator-valued second-order exchange: for example energy terms that do not annihilate the code (non-frustration-free), positional energy ladders written into the cell by its letters, or port couplings to excited "source" states whose energies and weights are set by local letters. Prove that the second-order source resolvent reproduces C19's filters (its n + 1 poles and Jacobi weights), or a variant that still gives independent flavour hopping, phases and chemical potentials with exact parity.
2. Realize the n + 1 poles with a fixed alphabet. The pole-count bound says longer cells must carry more poles; show how a cell of length growing with n carries n + 1 distinct accessible excitation energies using only local repeated terms (for example a ladder of sites whose local energies follow from position-independent rules), and how C19's coefficients become local amplitudes or follow from the geometry.
3. Write the full rule table: alphabet, every bond term with its coefficient, initialization, and the exclusion of illegal sectors. Count d0.
4. Verify by exact diagonalization on two and three cells at n = 1, 2 (and 3 if feasible) that the second-order effective bond matches the intended logical bond, with residual scaling as λ³, and that the code is exactly degenerate with exact parity.
