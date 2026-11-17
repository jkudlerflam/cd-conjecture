# Primary sources and provenance audit

These sources were opened during this contribution on 2026-11-17.
Local proofs and the finite tests are separate from literature inputs.

* [Lukyanov, Low energy effective Hamiltonian for the XXZ spin chain,
  Nucl. Phys. B 522 (1998) 533](https://arxiv.org/html/cond-mat/9712314).
  Read for the Gaussian low-energy description and symmetry-allowed
  corrections. The finite-a stability argument and its hypothesis list
  are given in ZERO_FIELD; the paper is not cited as a theorem for
  the generic asymmetric Ashkin–Teller lattice model.
* [Dugave, Göhmann and Kozlowski, Functions characterizing the ground
  state of the XXZ spin-1/2 chain in the thermodynamic limit,
  SIGMA 10 (2014) 043](https://sigma-journal.com/2014/043/sigma14-043.pdf).
  Read the stated integral-equation framework and its scope. It gives
  mathematical control of thermodynamic functions in the disordered
  regime, not a stand-alone universal CFT theorem for every Hamiltonian
  possessing U(1). The table labels the full XXZ infrared classification
  as established exact-solution physics.
* [Kohmoto, den Nijs and Kadanoff, Hamiltonian studies of the d=2
  Ashkin–Teller model, Phys. Rev. B 24 (1981) 5229](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.24.5229).
  The accessible primary abstract describes dualities, the XXZ map and
  a critical fan. Full text was access restricted. No blanket theorem
  about all unequal-colour couplings is inferred from the abstract.
* [O'Brien, Bartlett, Doherty and Flammia, Symmetry-respecting real-space
  renormalization for the quantum Ashkin–Teller model,
  Phys. Rev. E 92 (2015) 042163](https://arxiv.org/html/1507.00038v3).
  Read for the standard symmetric quantum AT model and its C=1 critical
  line. Its approximate real-space RG is not an exact solution of
  our generic unequal-velocity model.
* [Peotta et al., The XYZ chain with Dzyaloshinsky–Moriya interactions,
  J. Stat. Mech. (2014) P09005](https://arxiv.org/html/1403.4568).
  The paper studies a broader interacting parity family, with numerical
  evidence for liquid and ordered phases and possible split Ising
  transitions. It expressly leaves a meeting region unresolved.
  It is not used to close the whole h=0 face.
* [Zhang and Song, Exact phase diagram of the XXZ Heisenberg chain
  in a staggered magnetic field, arXiv:2609.29728v1](https://arxiv.org/html/2609.29728v1).
  The exact boundary and surrounding DMRG analysis are distinguished.
  The exact-gap proof used here is independently available in the
  repository's staggered_slice.md, including factors in Pauli units.

Prior project results read: `d2_gap_reduction/additive_charge.md`,
`staggered_slice.md`, C152 sections 8 and 11, C153 RESULT, C144 RESULT,
C157 RESULT, and the 2026-11-10 comparable-coupling search. The
Ashkin–Teller mapping and its exchange duality are prior results, not
claimed as discoveries here. The contribution is the auditable face
table, the explicit onsite implementation and seam of that duality,
rotation-invariant simplicity audit, and restricted stability statement
with its full leading operator selection rule.
