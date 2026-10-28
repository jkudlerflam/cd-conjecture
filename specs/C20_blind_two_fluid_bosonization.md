# C20: Stability of a multi-component Luttinger liquid with a "length" mode (blind check)

Closed book. Work only inside your current working directory. Do not read or list any file outside it. Write `PROGRESS.md` as you go and `RESULT.md` at the end.

## Setting

A one-dimensional lattice system has the following low-energy degrees of freedom.

- A "length" mode: a conserved U(1) charge (call it r) with density ν per unit cell, ν irrational and in (0, 1), described as a Luttinger liquid with parameters (v_ℓ, K_ℓ) and fields (φ_ℓ, θ_ℓ). Microscopically it is a free-fermion-like gas of "long cells" on a lattice of cells; its Fermi momenta are ±πν in cell units.
- A "payload" of n complex fermion flavours α = 1, …, n living on the cells. Flavour α has Fermi points p_{α,R} = θ_α + k_α and p_{α,L} = θ_α − k_α (cell units), with the centres θ_α and half-widths k_α chosen freely and generically. The payload is described by n Luttinger liquids with parameters (v_α, K_α), coupled by general quadratic (marginal, density–density and current–current) terms.
- Exact symmetries: one-site translation of the underlying lattice (which acts on cell variables as translation together with a shift tied to the length fluid, through an Ogata–Shiba-type gluing), conservation of r, and conservation of the total payload fermion parity. Do NOT assume conservation of individual flavour numbers.
- Couplings between length and payload are generic, local, small, and respect the symmetries above.

## Questions

1. Classify every local operator allowed by the symmetries that is relevant or marginal at the Gaussian point and does not oscillate (zero total momentum) for generic irrational ν and generic (θ_α, k_α). Include pairing terms, flavour-changing terms, pair tunnelling between flavours, terms that combine length-fluid vertex operators e^{2imφ_ℓ} with payload operators, and umklapp of every kind. Give their scaling dimensions in terms of the Luttinger parameters and their momenta.
2. State exact conditions on (ν, θ_α, k_α) under which no non-Gaussian relevant or marginal operator survives, and show that these conditions hold on an open set (or find an unavoidable resonance). Quantify the margins needed as a function of n.
3. Give the condition on the full quadratic stiffness matrix for the (n + 1)-boson theory to stay positive and gapless, and say what fails if it does not.
4. Conclude: is the low-energy theory an (n + 1)-component Gaussian theory with c = n + 1 on an open set of parameters? Describe the finite-size spectrum on a ring, including how zero modes behave at irrational ν.

Tag every claim as proved, derived with standard bosonization input, heuristic, or assumed. Do not put probabilities, odds or percent-confidence estimates in any file. Python (numpy and scipy) only for any check. Write prose plainly, without em-dashes.
