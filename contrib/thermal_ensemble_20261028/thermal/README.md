# Thermal comparison after massive matching

The remaining scalar RG problem is to compare the microscopic free energy
with the **dressed** multichannel Luttinger reference to accuracy o(T^2).
Smallness of a fixed interaction norm, or matching normalized current
responses, does not provide that accuracy. No microscopic central charge
or c-d violation is proved in this folder.

## Results

- [WINDING_ESTIMATE.md](WINDING_ESTIMATE.md) proves a single-scale thermal
  comparison with retarded endpoints and determinant contractions. Tree
  unfolding, a mixed Gram cofactor bound, and one temporal moment give
  an arbitrary fixed thermal power p>=3. Cauchy differentiation extracts
  a reference defect linearly. This closes the fixed-history winding step;
  multiscale localization and finite-temperature histories still need work.
- [THERMAL_LOCALIZATION.md](THERMAL_LOCALIZATION.md) constructs a projection
  using four allowed Matsubara frequencies. Its mismatch with line Taylor
  localization has power three for the temporal kinetic coefficient and
  power four for relevant/quartic coefficients. It supplies the other
  one-step thermal forcing term, with an explicit finite moment cost.
- [THERMAL_CRITERION.md](THERMAL_CRITERION.md) proves an abstract scalar
  comparison theorem. A residual scale modulus omega(s) tending to zero,
  together with a thermal comparison power p>2, suffices for o(T^2).
  Logarithmic convergence is allowed. It also proves a causal thermal
  history bound with the zero-temperature counterterm kept fixed, and
  identifies which vacuum-root, terminal-scale, and winding estimates
  the actual RG expansion must still establish.
- [QUADRATIC_MEMORY.md](QUADRATIC_MEMORY.md) proves exact thermal matching
  for a finite-orbital quadratic model with a positive massive block.
  The frequency residue in its retarded action changes the velocities.
  The eliminated scalar trace is exponentially small, so the correctly
  normalized action retains the original low-temperature coefficient.
- [checks.py](checks.py) and [checks.json](checks.json) check a two-band
  example, the geometric thermal modulus estimate, and a relevant causal
  recursion. They are controls of the stated lemmas, not tests of the
  full candidate construction.

The most relevant published blueprint is the scalar calculation in
[Giuliani–Mastropietro, Section 3.4](https://arxiv.org/pdf/1211.3030).
Its Ising theorem does not cover the present retarded multichannel action.
The source audit and the specific differences are recorded in the criterion.

## Reproduction

The script requires Python and numpy, and runs in about a second:

    python3 -B contrib/thermal_ensemble_20261028/thermal/checks.py

On the author's desktop, plain `python3` lacked numpy. The bundled runtime
was used; no packages were installed. The script requires no scipy.

For Delta=3, eta=0.8, and mu=0.2, the exact dressed speed is about 1.827.
The thermal coefficient converges to the dressed prediction as T decreases;
the saved data retain the precise values and errors. Discarding the
frequency residue gives a different, incorrect coefficient.

## Review and scope

The thermal worker independently checked the parent's length-tilted
embedding and the direct-attack theorem's new thermal periodization
comparison. The latter provides a beta-independent initial action up to
exponentially small terms under its polymer hypotheses. It does not supply
the later gapless scalar RG estimates in THERMAL_CRITERION.md.

Authored by an OpenAI GPT-6-family Codex agent on 2026-10-28. All claims
tagged proved are accompanied by derivations. Numerical checks are
identified separately. No human review or formal verification is claimed.
