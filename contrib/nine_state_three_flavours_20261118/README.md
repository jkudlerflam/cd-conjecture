# Nine states and three flavours: six roots, but a coupled instability

## Claim and status

**[Proved finite band algebra; numerical physical matching; conditional
one-loop instability.]** The same literal nine-state, Hermitian,
one-site-translation-invariant nearest-neighbour alphabet admits six
simple Fermi roots using the five proposed physical controls. All five
full-code projections are verified; no arbitrary symmetric logical
onsite matrix is inserted. An exact rational calculation certifies the
roots and a small box of independently varied physical controls.

At the selected rational point, the complete generated leading quartic
has three nonzero resonant pair transfers. Its coupled weak-coupling
flow lies in an invariant growth cone. This is an instability of the
three-mode Gaussian theory, **not a stable c=3 result**. No ultimate
microscopic Virasoro central charge is assigned. Neither this point nor
the bounded search excludes stable points elsewhere in the five-control
family. There is no arbitrary-c claim and no c>d-1 example.

The physical dimension is nine, including all three frozen marker
states and both charged carriers. The selected cell has five physical
sites and eight logical states. There is no vacancy or mobile lift.
Finite K, dressed sector selection and uniform infrared matching remain
unproved, exactly as for the preceding seed.

## Result at the rational point

In the quadratic energy unit `t0=lambda² alpha⁴ g²`, the physical
controls `(mu,eta,boundary,kappa_L,kappa_R)` are

```
(-3385, 9387, -7486, -2560, 11999) / 500000.
```

The cosine roots are approximately `(-.671732,.206478,.469500)`.
The three positive velocities are `(7.42972,3.93152,5.71505) 10^-7 t0`.
An exact conservative lower bound is `1.20214 10^-8 t0`.
The root box has independent control halfwidth `10^-9`; it is narrow.
This point is rational but is not a simple or broadly robust Hamiltonian.

With dimensionless transfer coefficients `c_ij` and forward dimensions
`2+D_ij`, the leading coefficients in edge order 12,13,23 are

```
D/g² = (-195.766578, -7.451236, -362.961693),
c/g² = (  53.842721,-70.814527,  189.004243).
```

Their sign product is negative. A field rephasing makes all three
transfers negative; then the coupled flow increases all their magnitudes
and makes every D more negative. The smallest initial negative-D margin
is 7.45123 g². The proof includes transfer through the third channel,
not three independent applications of the earlier two-flavour test.

## Files and checks

- [BOND.md](BOND.md): all local rules, physical controls and exact projections.
- [BANDS.md](BANDS.md): rational determinant, six-root certificate and velocities.
- [RG.md](RG.md): complete generated vertex and coupled instability argument.
- [LIMITS.md](LIMITS.md): error scales, unresolved hypotheses and all-n limits.
- [PROTOCOL.md](PROTOCOL.md): bounded scope and initial benchmark.

The literal and scalar-resolvent implementations agree through fourth
conversion order within 5.6e-17. Full-code projection errors are at most
4.5e-16. All 21 root-pair sums, including repeated roots, are audited;
only the three opposite-momentum pairs collide at the selected point.
The complete 15 by 15 quartic matrix and all 21 allowed ordered quartic
vertices are saved in `audit.json`.

The optional exploration uses 108 local solves with prescribed cosine
roots. All produced six-root symbols. Sixty generic-momentum points
reached a twentyfold growth threshold in the coupled flow; 48 points
with extra or nearly resonant umklapp were excluded from that generic
flow test. This is a bounded diagnostic, not a no-go theorem.

Author and AI assistance: OpenAI Codex, GPT-6, for derivation, code and
literature audit. No separate-agent or independent human review is
claimed. The independent implementations are an internal cross-check.

## Reproduction

From the repository root, with numpy and scipy and one BLAS thread:

```sh
OPENBLAS_NUM_THREADS=1 python contrib/nine_state_three_flavours_20261118/audit.py
OPENBLAS_NUM_THREADS=1 python contrib/nine_state_three_flavours_20261118/explore.py
```

The audit takes about .2 seconds here, and the optional bounded search
about one second. `exact.py` can also be run separately. The scripts
reuse the unchanged carrier and literal-word implementations from
[the two-flavour audit](../ten_state_resonance_audit_20261117/README.md).
No many-body size scan, DMRG or extra worker is used.

This contribution resolves the proposed three-flavour band-synthesis
test positively and supplies a located interacting instability. The
coordinator owns any changes to the repository's top-level summaries.
