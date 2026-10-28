# Independent audit of the n=23 band certificate

[checked] The C33 finite-loop integer generator and root certificate were
rerun from their source scripts, in a temporary copied workspace. The new
F, M, EF, EM arrays and all integer metadata agree exactly with the stored
`codex/work_C33/full_loop_n23.npz`. The complete endpoint-inertia and pair
tests pass again. The reproducible command is

```
python3 contrib/counterexample_campaign/route_a/audit_c33.py
```

The small output is [`c33_audit.json`](c33_audit.json). Original contribution
files are never modified. Generation took about 20 seconds and the root
certification about 3 seconds on this run.

[checked] The mathematical audit covered the following acceptance steps.

1. The scalar complex Jacobi solve implements the correct loop sign and
   has infinity-norm contraction at most 1/2+t_d/4. Fixed-point truncation
   and source errors are substantially smaller than the 10/S radius used
   in the entry enclosures.
2. The factor-of-two normalizations of hopping and onsite matrices agree
   with C33's source-contraction formula. Their sum EF+EF^T+EM bounds a
   complete Bloch-symbol error.
3. The fraction-free conditioning identity is checked over integers.
   Rationalized numerical eigenvectors merely propose a congruence.
   Its invertibility and the inertia signs are certified by integer
   diagonal-dominance inequalities. The code uses ||Q||^2<=2n, obtained
   from its Gram trace bound, when propagating errors.
4. Each of 46 disjoint intervals has an inertia change. Since the
   determinant of a 23-channel nearest-neighbor symbol has Laurent degree
   at most 46, these exhaust its roots and are all simple.
5. All 1081 unordered pairs with replacement are checked. Sorting uses
   rational orientation tests; subtracting the angular interval widths
   from a coordinate-distance lower bound gives a valid chordal, hence
   angular, separation bound.
6. For |g|<=1/2 and the two-carrier A>=3I, the even resolvent expansion
   bounds its remainder by alpha^4 g^4/4. Quadratic extraction gives the
   stated 4n alpha^4 g^4 Bloch bound. The rational trace of the conditioning
   matrix then verifies the error budget at g=2^-365.

[clarification] In C33 Section 9, the product-error argument is safest
when the 10/S error is read as a bound on the **complex modulus** at each
entry, rather than a separate bound of that size on each real component.
This stronger interpretation follows from the same solve estimates:
the two real floor errors have combined modulus below 2/S per step,
their geometric accumulation is below 5/S, and tail and second-source
errors keep the total below 10/S. The quadratic product error is then
at most 100/S^2, as used in the code. No numerical radius change is needed.

[checked result] The quadratic reference at n=23 has 46 simple Fermi
points, delta_2>0.0140193 and delta_4>3.1558e-7, including the finite-g
quadratic correction at g=2^-365. The minimum remaining endpoint
congruence margin exceeds 1.9305e-8. This removes the finite-band
qualification in the older `explicit_candidate/README.md`.

[scope] The proof is of the **quadratic projection** of the second-order
Hamiltonian. It does not discard its nonquadratic remainder, certify a
thermodynamic central charge, or prove that frozen cells meet an all-ring
CFT requirement. The same certificate holds for any smaller g if the
common field is rescaled proportionally to g^2. This freedom is useful
for making the normalized nonquadratic remainder arbitrarily small after
the finite n=23 reference momentum and velocity margins have been fixed.

Audit by a separate Codex GPT-6-family worker, 2026-10-28. This is an
independent AI check, not a human review.
