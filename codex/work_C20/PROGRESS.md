# C20 progress

[Assumed] This is a closed-book check using only the setting in the prompt. No files outside the current directory have been read or listed. The directory initially contained no visible research files.

[Derived with standard bosonization input] The useful organizing object is the integer lattice of local vertex operators, with length-charge conservation and total payload parity imposed before any momentum test. A finite list of named scattering processes would miss higher harmonics when stiffnesses are unrestricted.

[Assumed] The action of underlying translation on the compact bosons, including the stated gluing, has not been specified explicitly. I will distinguish a cell-coordinate momentum calculation from any conclusion that requires a particular physical translation action.

[Proved] An open stability claim must concern only the finitely many vertices of dimension at most two in a region where the Gaussian scaling metric stays uniformly positive. Excluding resonances of every order is a stronger condition and need not define an open set.

[Derived with standard bosonization input] With vertex coefficients `(2m,a_1,...,a_n;b_1,...,b_n)`, length-charge conservation sets the length dual coefficient to zero. Payload locality requires `a_alpha = b_alpha mod 2`, and total parity requires `sum b_alpha = 0 mod 2`. In a cell-coordinate convention the momentum is `2 pi m nu + sum(a_alpha k_alpha - b_alpha theta_alpha)`, up to an overall sign. The diagonal Gaussian dimension is `m^2 K_length + sum(a_alpha^2 K_alpha + b_alpha^2/K_alpha)/4`.

[Proved] At the decoupled free metric, squared norms of allowed nonzero vertices are positive multiples of four. If the full scaling metric is bounded below by `(1-epsilon) I` with `epsilon < 1/3`, all possibly relevant or marginal vertices have squared norm four or eight. Including one representative of each conjugate pair gives

`M(n) = 1 + 8 n + 32 binomial(n,2) + 96 binomial(n,3) + 128 binomial(n,4)`.

[Proved] A finite union of resonance hyperplanes cannot fill the momentum parameter domain. A quantitative detuning bound of order `1/M(n)` is sufficient for a nonempty open complement. Requiring irrationality itself restricts that open complement to the irrational-density subspace; irrational densities do not constitute an open subset of the real line.

[Proved] The independent integer enumeration in `check.py` agrees with the shell-count formula for `n = 1,...,8`. The resulting candidate counts are `9, 49, 217, 737, 1961, 4369, 8569, 15297`.

[Derived with standard bosonization input] The general scaling metric is `G = H^(-1/2) |H^(1/2) J H^(1/2)| H^(-1/2)`. Its matrix identity `G J G = J` and its decoupled limit were checked numerically, with maximum residuals `4.98e-15` and `6.66e-16`. These numerical checks supplement the derivation in `RESULT.md`.

[Proved] Quadratic energetic stability requires strict positivity of the full stiffness matrix, equivalently the block Schur-complement conditions given in `RESULT.md`. An entrywise perturbation bound of order `1/(n+1)` suffices to preserve a fixed positive stiffness margin. This is independent of the momentum-detuning bound of order `1/M(n)` near the free scaling metric.

[Derived with standard bosonization input] The ring spectrum has oscillator towers and a quadratic winding-sector energy. At irrational thermodynamic length density, finite rings use integer charges with offset `r - nu L_c`. These offsets can contribute at order `1/L_c` and must be retained when interpreting the Casimir term.

[Assumed] The exact physical translation character, the winding-sector selection rules, and the conversion between cell circumference and physical circumference still require an explicit gluing map. The prompt does not supply one. `RESULT.md` gives the full cell-coordinate answer, the general replacement test for a supplied gluing, and a conditional Gaussian stability conclusion.

[Proved] Execution status: the initial `python` command failed with `zsh:1: command not found: python`. The full deterministic check subsequently passed with `python3 check.py` and NumPy. A text audit initially flagged the en dash in the name Ogata-Shiba; that typography was replaced and the final text audit passed.

[Assumed] Work is complete within the supplied information. Deliverables are `RESULT.md`, this progress record, and the reproducible check `check.py`. No project memory outside the current directory was accessed or modified.
