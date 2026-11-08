# C154: c ≤ 1 for qubit bonds without symmetry

Your folder is `../work_C154/`. Read `specs/_qubit_c_context.md` first.

## Setting

A generic point of B_2 has no symmetry (no on-site unitary or antiunitary symmetry, beyond translation). Without symmetry, every gapped phase in one dimension is trivial, so a generic one-parameter path need not close the gap, and gapless points occur only on fine-tuned loci. The typical generic critical point is an Ising critical endpoint of a first-order surface (codimension 2, like the liquid–gas critical point: transverse-field Ising plus a longitudinal field). A CFT fixed point with k relevant translation-invariant scalar operators (other than the identity) needs k tunings, so c > 1 points would have to be multicritical loci of codimension at most 8 inside B_2.

## Tasks

1. Develop the codimension and RG argument into a derivation. For each candidate c > 1 CFT reachable from a qubit chain (Ising × boson c = 3/2, Ising² × boson, SU(2)₂, tricritical Ising × Ising, two bosons, and so on), count the relevant translation-invariant scalar operators that a generic qubit bond generates, using lattice momentum and the one-qubit-per-site constraint (LSM-type anomalies apply only with symmetry, so say exactly what survives without it). Then decide whether such a fixed point can be realized inside the 8-dimensional B_2 at all, and where its multicritical locus would have to sit (for example at the closure of a symmetric stratum, where emergent symmetry lowers the codimension).
2. Show, or refute, that every gapless point of a symmetry-free qubit bond has an emergent symmetry at its fixed point and is continuously connected, along its critical locus, to a point of a symmetric stratum (U(1), parity, staggered charge, or a singular stratum), so that its c is fixed by the symmetric analysis (C152, C153, the U(1) class). Treat critical endpoints and their emergent Z2.
3. Analyse the most dangerous regions directly: the YX and ZY directions at the XZ point (C144 found a slow sector with c ≈ 1 at L = 64; collect the L = 128 to 256 outputs of Typhon jobs 17580192 and 17580432 from /scratch/lustre/jkudlerflam/c144/ if they have finished, and decide c there), and intersections of an Ising critical surface with a c = 1 region where a 3/2 coexistence could occur.
4. Calibrate numerically on the critical path (DMRG, ED). State the outcome as a conditional derivation of c ≤ 1 for all of B_2, with every remaining assumption tagged.
