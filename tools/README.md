# Tools

Digests of techniques that might bear on the c–d conjecture, written in our own words with every claim tagged [read] (source opened, with arXiv number and section), [recalled], or [own] (the digest author's computation, unchecked unless stated). They were produced by Claude agents in October 2026 and have not been refereed.

## entanglement_bootstrap/

- `gapless_bootstrap.md` covers the 1+1D and CFT side of the entanglement bootstrap: Lin–McGreevy's vector fixed-point equation and reconstructed Hamiltonian (2303.05444), conformal geometry from entanglement and the stationary c_tot (2404.03725), the chiral Virasoro algebra from a single wavefunction (2403.18410), the Li–Lin–McGreevy search that tests c–d directly (2509.04596), modular commutators and the finite-dimension no-go for chirality (2206.00027, 2408.10306, 2609.09019), the Markov gap, and lattice Kac–Moody constructions (2206.01656).
- `gapped_bootstrap.md` covers the gapped program of Shi, Kato, Kim, McGreevy and collaborators: the axioms, the merging lemma, information convex sets and the simplex theorem, Markov block structure, and the chiral central charge from a single wavefunction.
- `mcgreevy_notes.md` covers John McGreevy's lecture notes and reviews: entanglement c-theorems, Lieb–Robinson bounds, anomalies and Lieb–Schultz–Mattis constraints, non-invertible symmetries, Luttinger's theorem as an anomaly, and tensor-network bounds.

## The common conclusion

Every quantity in this program that equals c is an entropy combination, a commutator controlled by modular-Hamiltonian variances, or a conditional mutual information. At a fixed region size each is bounded by a power of log d. Used alone, these tools therefore give c ≲ (onset length of conformal behaviour) × log d, and the SU(d)₁ chains show that the onset length must grow at least like d/log d.

What the program does supply is (i) velocity-free normalizations, (ii) exact multi-scale consistency relations, (iii) generators written as sums of local terms on the same (ℂ^d)^⊗L, and (iv) a template for finite-dimension no-go theorems. The most promising unexplored lever is building Kac–Moody currents from the wavefunction, where a count tied to the d² − 1 traceless on-site operators could enter through the Sugawara construction.
