# Exact-diagonalization screen

OpenAI Codex, GPT-6 family, 2026-11-10. These are full microscopic spin-chain
results, with no blockade or particle-number projection. No point establishes
c>1, and the finite-ring null readings do not exclude it. The most striking
preferred-ratio flags at L=10,12 weaken at L=14,16. The best rational attractive
point gives readings near c=1 over these sizes.

The frozen [protocol](PROTOCOL.md) produced 145 points at L=10 and 12, with
both exact parity sectors and two low eigenvectors per sector: 290 completed
point-size records in [screen.jsonl](screen.jsonl). The separately registered
[follow-up](FOLLOWUP.md) produced six selected points at L=14,16 and four
controls at L=16,20, using four low eigenvectors per sector: 20 records in
[followup.jsonl](followup.jsonl). There are no unfinished ED jobs. Optional
L=18 calculations were not run; open-chain follow-up is recorded separately.

## Calibration and detector limitations

For periodic rings the size diagnostic is
`c_size=3*(S_half(L2)-S_half(L1))/log(L2/L1)`. The profile diagnostic is three
times the fitted slope of S against `log[(L/pi)sin(pi*l/L)]`, using cuts in
`[L/4,3L/4]`. These are finite-size diagnostics, not measurements of an
established conformal limit.

| Control | c_size, 10→12 | c_size, 16→20 |
|---|---:|---:|
| Nearest-neighbour XX, c=1 | 0.9804 | 0.9927 |
| Critical transverse-field Ising, c=1/2 | 0.4976 | 0.4991 |
| Gapped Ising | -0.0042 | -0.00009 |
| Two interleaved XX chains, c=2 | 3.0742 | 1.9402 |

The two-XX control has range-two interactions and is **not** a target
counterexample. Its L=10,12 legs have incompatible odd/even lengths; the
initial pair fails c=2 calibration. At L=16,20 both legs have even lengths,
and the control separates c=1 from c=2. Its profile estimates still vary,
from 2.073 to 2.155, illustrating cut and shell sensitivity.

The known r=40 reference is almost invisible to this small-ring screen:
its mean density is 0.002189→0.002242 and c_size is 0.0531 at 10→12. This is a
dressed-vacuum finite-size state, not the gas seen in the much larger
open-chain work. Consequently low ED entropy is not a veto on this mechanism.

## Selected microscopic points

All bridge points use
`Jx=1/2, Jy=-1/2, Jz=r/4, Dz=u/(4r), hz=-r/2+k/(4r)`.
Each nonzero-field bridge has four grouped terms and six expanded Pauli
strings. The attractive points have five strings and four grouped terms
in the `{XX+YY, XX-YY, ZZ, XY-YX, Z}` basis.

| Short label | Grid label / parameters | (Jz,Dz,hz) | Coefficient spread R | Exchange singular values |
|---|---|---|---:|---|
| A | r=3/2,u=3/2,k=3 | (3/8,1/4,-1/4) | 2 | (3/4,3/8,1/4) |
| B | r=2,u=3/2,k=4 | (1/2,3/16,-1/2) | 8/3 | (11/16,1/2,5/16) |
| C | r=3/2,u=2,k=1 | (3/8,1/3,-7/12) | 7/4 | (5/6,3/8,1/6) |
| D | r=4,u=1,k=4 | (1,1/16,-7/4) | 28 | (1,9/16,7/16) |
| E | Attractive rational | Jx=1/2,Jy=-1/4,Jz=-1/2,Dz=1/2,hz=0 | 2 | (0.89039,0.5,0.14039) |
| F | Existing exploratory reference | phi=1,Delta=0.7,V=-2.4,t=1,hz=0 | 7.7666 | (0.85,0.6,0.15) |

A–C and E meet the preferred literal coefficient spread. F is in the
exploratory band. D is a continuation bridge outside it. Singular values
alone conceal D's weak antisymmetric interaction: its antisymmetric/symmetric
Frobenius ratio is 0.07217, compared with 0.44173,0.30619,0.58897 for A,B,C.
E's singular-value spread is 6.34 despite literal coefficient spread 2.
All these invariants and the field norm are saved per record.

There was no precision field fit in this stage. The bridge ladder spacing
was 1/(4r); the selected fields are the exact fractions shown. The bare
combinations `4Jz+2hz` are 1,1,1/3,1/2 for A–D. Their exact `Jx+Jy=0` is a
chosen symmetry condition, not subtraction of fitted decimal coefficients.
F is a previously supplied non-rational reference, not a simplification.

| Point | Raw c_size 10→12 | Momentum-resolved c_size 14→16 | c_profile 14→16 | Density 14→16 |
|---|---:|---:|---|---|
| A | 1.9364 | 0.9974 | 1.0701→0.8628 | 0.43227→0.43498 |
| B | -1.8025 | -0.0548 | 1.4983→1.3328 | 0.43122→0.44003 |
| C | 2.2317 | 0.9167 | 1.1604→1.2943 | 0.37603→0.38668 |
| D | 1.3100 | -0.5349 | 2.2494→1.8195 | 0.30591→0.28758 |
| E | 1.0248 | 1.0272 | 1.0231→1.0218 | 0.50000→0.50000 |
| F | 0.9916 | 0.9927 | 1.0073→1.0047 | 0.50000→0.50000 |

A was the cleanest initial flag, with momentum eigenstates at both initial
sizes; its excess does not persist in this larger-ring comparison. B was
retained despite a negative initial size slope, as registered. C changes
ground parity and shell occupancy. D changes density strongly: its L=12
density was 0.15923, so these sizes do not compare a settled dilute gas.
E's full seven-point attractive grid has its largest initial c_size at this
rational point, about 1.025; none supplies a clean larger-than-one signal.

## Momentum degeneracy and numerical checks

Single-vector Lanczos can return a non-translation eigenstate in a degenerate
ground space, even when its requested spectrum misses the degenerate partner.
For example C at L=10 has translation modulus 0.892 and a nonuniform density
profile in the original output. Large entropy from such a vector is not
interpreted as a conformal signal.

The follow-up retains the raw vector's entropy and resolves translation in
the span of all translations of the returned low ground vectors. This also
recovers missed momentum partners. It reports every retained branch, selects
the smallest momentum index without optimizing entropy, and checks closure,
energy residuals and repeated momentum. Raw and resolved observables remain
separate. The original Stage A physics outputs were not replaced.

All 20 selected ground states pass momentum resolution. Four **non-ground**
odd-sector XX/two-XX control states retain a repeated-momentum internal
degeneracy and are explicitly marked unresolved; their arbitrary entropies
are not used for calibration. Ground-sector selection always uses energy.
The returned low-spectrum gaps are not certificates of multiplet completeness.

A review found that the original `ground_gap` field was only the spacing of
the returned Lanczos eigenvalues, even when translation completion recovered
a missing ground partner. That original field is retained unchanged as legacy
data, with explicit semantics and a duplicate named `raw_lanczos_ground_gap`.
Added derived fields count the resolved ground vectors across both parity
sectors. B and C at L=14, and D at L=14 and16, each have at least two ground
vectors within the recorded energy tolerance, so their
`gap_counting_resolved_multiplicity` is zero. In particular D's raw spacings
6.49e-5 and1.21e-4 are **not** gaps counting ground multiplicity.
`gap_above_manifold` is explicitly null for every record because completeness
of the remaining low spectrum is not certified. A separately labeled raw
candidate spacing above the observed manifold remains available. The
[metadata audit](gap_metadata_audit.json) verifies that this correction added
derived metadata without changing any original field or running eigensolves.

Independent L=6 Kronecker-product assembly agrees exactly with the bit-built
matrix, for both parities and both open/periodic boundaries. Maximum Stage A
eigen residual is 1.20e-9. Follow-up maxima are 2.16e-10 for raw eigenvectors,
1.23e-9 for resolved momentum vectors, and 3.54e-13 for translation-subspace
closure. These checks test the finite Hamiltonian and state calculation;
they do not certify infrared conformality, exclude phase separation, or
rule out an ordered cat state. Density structure factors and profiles are
saved, but a uniform translation eigenstate can still be a cat state.

## Reproduction and provenance

From the repository root, with NumPy and SciPy installed:

```sh
python3 contrib/d2_order_one_20261110/ed_screen.py --validate
python3 contrib/d2_order_one_20261110/ed_screen.py --config contrib/d2_order_one_20261110/grid.json --output /tmp/order_one_screen.jsonl
python3 contrib/d2_order_one_20261110/ed_followup.py --config contrib/d2_order_one_20261110/followup_grid.json --output /tmp/order_one_followup.jsonl
```

The actual runtime was `/private/tmp/cd-conjecture-venv/bin/python`, NumPy 2.5.4,
SciPy 1.18.1, with one BLAS thread. Completed JSONL rows are resumable only when
couplings, geometry, sectors, eigenvector count and tolerance match. An
interrupted final fragment is retained separately before truncation.

After Stage A, the grouped-term count was corrected to the explicit paired
operator basis above, and cancellation metadata were added. A record-level
amendment identifies that metadata-only change. Later edits hardened resume
handling. None changed the assembly, eigensolver, entropy or observables.
Exact executed source versions, including the pre-correction follow-up wrapper,
are archived under `executed_sources/` and
match the source hashes recorded in the raw files; [ed_provenance.json](ed_provenance.json)
audits these hashes and the unchanged physics functions. Use the current
scripts for new runs because they include the safe resume fixes.

Checkpoints for L>=16 are stored only in
`/private/tmp/d2-order-one-checkpoints`; no vectors or checkpoints are in Git.
The ED raw outputs total about 2 MB. Every screened point, including null
results, is retained.
