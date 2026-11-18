# Quantitative control of finite markerless mobile sectors

2026-11-18. This contribution concerns the corrected physical d=10
nearest-neighbour, Hermitian, one-site translation-invariant model in
[`markerless_d9_d10`](../markerless_d9_d10/README.md). Coefficients are
fixed as chain length grows. It does not establish a microscopic
unbounded-c family or a unitary compact Virasoro CFT.

## Results and scope

- **[Proved, integer certificate]** At n=9, canonical N=20 occupied
  particles and V=2 vacancies (rho=10/11), J/t=1/25 and the original
  selector midpoint, the (10,10) word beats (9,11), (20), and (6,7,7).
  Their certified lower gaps are respectively **0.0001200227 t,
  0.0032431549 t, and 0.0039278859 t**. The same bounds hold throughout
  1.77<=nu<=1.7703. A small interval around J/t=1/25 is also enclosed.
  These are specified finite competitors, not an all-word or
  thermodynamic ground-sector theorem.
- **[Proved]** At zero field and port, charged heads cannot lower any
  legal word below its vacuum energy when
  Delta/2>=max(1,|g|,|t_d|) and D>=0. This holds for the certified n=9
  masses at every density and J/t. It is a domination inequality, not
  a charged excitation gap.
- **[Proved, finite extension]** At J/t=1/25 the same three word
  comparisons survive arbitrary charged/bit states and finite ports
  with sqrt(2)|lambda|+f_ext/2<=1/1000, for the external field below.
  The smallest remaining certified gap is 0.0000400227 t. This bound
  is deliberately not extended to infinitely many spectator cells.
- **[Proved formula; numerical coefficients]** The actual markerless
  clock has an explicit contact form factor and a bulk-subtracted
  second-order response. Balanced defects and the weak change to more
  length-9 cells have unfavorable corrections; the change to fewer
  length-11 cells has a favorable correction in the evaluated baths.
  The earlier biased-clock favorable-sign theorem does not transfer.
- **[Proved algebra]** Moving only the conserved charge field outside
  the contact embedding removes its bit-position error exactly at
  every J/t. To match the leading frozen band target at fixed density,
  set f_ext=alpha_infinity(rho) f9. Other material energies stay contact
  weighted. This leaves the zero-port vacuum selector unchanged but
  does not control the other dressed band errors.

The detailed statements, integer proof, and the remaining uniform
relative-curvature inequality are in [PROOF.md](PROOF.md). The field
normalization, charge algebra, term counts, and coefficient hierarchy
are in [FIELD.md](FIELD.md).

## What the calculation adds

The previous check compared one cell with a headless ring. Here the weak
balanced defect is tested in the literal moving-vacancy Hamiltonian,
with both head-number directions included. The larger, genuinely weak
head-number comparisons are also calculated to second order:
9 length-10 cells versus 10 length-9 cells at N=90,V=1, and 11
length-10 cells versus 10 length-11 cells at N=110,V=1. Their single
vacancy responses reduce exactly to a cycle Green function. These two
densities are stated separately and are not half-occupancy evidence.

No inverse gap uniform in volume is assumed. No extensive error from
unchanged cells is charged to one defect. The finite full-matrix
certificate and its finite-port norm bound are explicitly confined to
their finite ring. At half occupancy the precise remaining requirement
is a uniform bound on the curvature of the energy difference for a
same-length local replacement, with the interacting vacancy fluid left
intact. This is the most decisive next calculation.

The selector gap near n=9 is about 1.82e-4 per material site. The logical
band's normalized 2^-35 allowance is a different budget. There is still
no certified J/t allowance for the dressed mobile band or a complete
infrared stability theorem. The field modification also does not make
the model simple: it keeps 50 Hermitian off-diagonal pairs, adds 12
diagonal bond entries, and the existing n=9 band schedule has a nonzero
coefficient spread exceeding 10^407 at half occupancy.

## Reproduction and verification

Author and AI assistance: OpenAI Codex, GPT-6 Astra / Max research
assignment. No subagents, background jobs, DMRG, or large scans were used.
No independent human or second-agent review has been performed.

With Python and `numpy`, `scipy`, `mpmath`:

```sh
python3 contrib/markerless_mobile_control_20261118/checks.py
python3 contrib/markerless_mobile_control_20261118/checks.py --verify-only
```

The first command regenerates four positive integer witnesses and the
small [checks.json](checks.json) report. The second checks the saved
[witnesses.npz](witnesses.npz) without solving the four many-clock
eigenproblems. The script forces one numerical-library worker. The
recorded local run took under one second; allow several seconds on
other machines. Largest matrix: 61,740 states, 581,560 nonzero entries.
All files written are in this folder.

Floating eigensolvers only find positive trial vectors. Each energy
interval is subsequently certified by integer arithmetic and the
stoquastic minimum/maximum local-energy inequalities. The script also
checks the literal material sources, the charge commutator, the exact
clock form-factor identity, the one-vacancy Green function, and a
rational enclosure of the original selector value. The finite-volume
response table is floating-point evidence, not an interval certificate.

## Relation to the project

This addresses assignment 4 of the focused 2026-11-18 round. It closes
some finite microscopic checks missing from the corrected markerless
model and identifies the actual relative response that a thermodynamic
selection proof must bound. Arbitrary-n band algebra is owned by the
parallel contribution; it is neither assumed solved nor edited here.
Shared root summaries are left to the coordinator; proposed text is in
[SUMMARY_FOR_COORDINATOR.md](SUMMARY_FOR_COORDINATOR.md).
