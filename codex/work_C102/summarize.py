"""Render the C102 numerical record as a short Markdown report."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
r = json.loads((ROOT / "numerics.json").read_text())
out = ["# C102 numerical checks", "", "2026-10-29", "",
       "[Numerical] Reproduce with `python3 checks.py`, then `python3 summarize.py`, from this folder. "
       "The first script writes numerics.json and asserts the consistency checks. The second writes this report. "
       "Python, NumPy, and SciPy are the only computational tools. No network commands or cluster jobs were used.", "",
       "[Numerical] The sparse ring matrix is assembled from the two-site bond in a base-d product basis. "
       "XXZ uses total Sz=0 and its known zero-field ground-state sector. Spin-1 SU(2) multiplets all contain Sz=0, "
       "so the spin-1 sector calculation includes a global ground state. Generic rotations and gauges use the full space. "
       "Independent small full-space comparisons check both sector implementations.", "",
       f"[Numerical] Recorded {r['number_recorded_eigenpairs']} eigenpairs, with largest matrix dimension "
       f"{r['max_matrix_dimension']}. The largest eigenpair residual was {r['checks']['max_eigenpair_residual']:.3e}. "
       "These are modest local checks. No large-volume computation was performed.", "",
       "## Estimators and Bethe input", "",
       "[Proved] For lengths L<K, the bulk-canceling two-size estimator is", "",
       "    A_pair=(6/pi)[E(K)/K-E(L)/L]/[L^(-2)-K^(-2)].", "",
       "[Numerical] The second estimator fits E(L)/L=e+a/L^2+b/L^4 and sets A_fit=-6a/pi. "
       "Where the exact bulk energy is known, the third estimator is A_L=(6/pi)L[L e-E(L)]. "
       "The polynomial fit is a diagnostic and is not a controlled treatment of logarithmic or noninteger-power corrections.", "",
       "[Assumed, Bethe input] XXZ uses A=pi sin(gamma)/(2 gamma) and delta=cos(gamma). "
       "Its bulk energy is evaluated independently by quadrature of", "",
       "    e(delta)=delta/4-sin(gamma) integral_0^infinity",
       "              sinh((pi-gamma)x)/[sinh(pi x) cosh(gamma x)] dx.", "",
       "[Numerical] The integrand is implemented with decaying exponentials to avoid overflow. "
       "Bethe testing here means this bulk-energy quadrature and the closed velocity formulas, "
       "not a finite-size Bethe-root solver. TB and ULS use the inherited thermodynamic data specified in PROOF.md.", "",
       "## XXZ", "",
       "[Numerical] Every row uses L=6,8,10,12,14,16, in spin units J=1. All optimized M values are 1/2.", "",
       "| delta | Bethe A | Pair 14,16 | Polynomial fit | Bulk subtraction L=16 |",
       "| :-- | --: | --: | --: | --: |"]
for m in r["models"]:
    if m["name"].startswith("XXZ"):
        out.append(f"| {m['name'][4:]} | {m['exact_A']:.9f} | {m['diagnostic']['pairs'][-1]['A']:.9f} | "
                   f"{m['diagnostic']['quadratic_in_inverse_L_squared']:.9f} | {m['rows'][-1]['A_bulk_subtraction']:.9f} |")
out += ["", "[Numerical] The decisive delta=-1/2,0,1/2 midpoint defect from bulk subtraction is", "",
        "| L | Midpoint defect |", "| --: | --: |"]
for x in r["XXZ_Jensen_identity"]:
    out.append(f"| {x['L']} | {x['A_midpoint_defect']:.12f} |")
out += ["", "[Proved conditional] The Bethe limit is 1-9 sqrt(3)/16=0.0257214207425. "
        "All sampled signs agree. The finite-size defect is not monotone in L. "
        "The bulk Jensen defect is 0.00644964023537; at L=16 the ring Jensen defect is 0.102357202698. "
        "Their difference after the extensive subtraction has the positive sign in Eq. (1) of PROOF.md.", "",
        "## XYZ", "",
        "[Numerical] The convention is h=sum_a J_a S_a tensor S_a. Full spaces at L=6,8,10,12 are used. "
        "The first two bonds are rotated XX chains. The last two are massive anisotropic free XY chains. "
        "The middle bond is a rotated easy-axis XXZ chain, with a massive thermodynamic sector.", "",
        "| (Jx,Jy,Jz) | Thermodynamic A input | Pair 10,12 | Polynomial fit |",
        "| :-- | --: | --: | --: |"]
for m in r["models"]:
    if m["name"].startswith("XYZ"):
        out.append(f"| {m['name'][4:]} | {m['exact_A']:.6f} | {m['diagnostic']['pairs'][-1]['A']:.9f} | "
                   f"{m['diagnostic']['quadratic_in_inverse_L_squared']:.9f} |")
out += ["", "[Numerical] The gapped (1,1/2,1/2) model still gives an apparent A near 0.87 at these lengths. "
        "It is not evidence for a nonzero asymptotic Casimir coefficient. "
        "The two anisotropic XY endpoints have A=0, while their midpoint is (1,1,0), with A=1. "
        "This additional convexity failure uses the gapped extension of A, unlike the principal XXZ example.", "",
        "## Spin-1 bilinear-biquadratic", "",
        "[Numerical] Here h_b=X1+b X1^2, X1=S_i.S_{i+1}. For b<1 use L=4,6,8,10. "
        "For b>=1 use L=3,6,9, compatible with the period-three ULS vacuum. "
        "The last column is a fit diagnostic, not an asymptotic determination.", "",
        "| b | Thermodynamic A input | Last pair | Polynomial fit |",
        "| :-- | --: | --: | --: |"]
for m in r["models"]:
    if m["name"].startswith("BB"):
        a="unresolved here" if m['exact_A'] is None else f"{m['exact_A']:.9f}"
        out.append(f"| {m['name'][5:]} | {a} | {m['diagnostic']['pairs'][-1]['A']:.9f} | "
                   f"{m['diagnostic']['quadratic_in_inverse_L_squared']:.9f} |")
s = r["spin1_critical_midpoint_defects"]
out += ["", f"[Numerical] For b=1,1.5,2, the midpoint defect is {s['pair_L6_L9']:.9f} "
        f"from the last pair and {s['quadratic_in_inverse_L_squared']:.9f} from the fit. "
        "Both support the XXZ direction, but three small lengths do not establish the limiting sign. "
        "The b=1/2 fit is particularly unstable, and the nonzero Haldane-phase values are finite-size artifacts. "
        "AKLT at b=1/3 supplies an exact zero check.", "",
        "## Seeded critical bonds", "",
        "[Numerical] Four independent triples were generated from XXZ with random J and anisotropy, "
        "a common on-site SU(2) rotation for each triple, and independently chosen scalar and Hermitian coboundary gauges. "
        "Each midpoint is the literal average of the two full matrices. Each entire segment is known to be critical. "
        "All computations use full spaces at L=6,8,10. Seed=102 fixes the run. "
        "This is a constrained critical-family sample, not a claim that unconstrained random Hermitian bonds are critical.", "",
        "| Triple | J endpoints | delta endpoints | Exact midpoint defect | ED defect L=10 |",
        "| --: | :-- | :-- | --: | --: |"]
for i,p in enumerate(r["random_pairs"]):
    js=", ".join(f"{x:.6f}" for x in p["J"])
    ds=", ".join(f"{x:.6f}" for x in p["delta"])
    out.append(f"| {i} | {js} | {ds} | {p['exact_midpoint_defect']:.9f} | {p['ED_midpoint_defect_L10']:.9f} |")
out += ["", "[Proved] Gauge-optimized M for each transformed bond is inherited from its XXZ source, J/2. "
        "The JSON field raw_half_width is the unoptimized width and is not substituted for M. "
        "[Numerical] Rotations, scalar shifts, and coboundaries preserve the checked ring energies to "
        f"{r['checks']['max_unitary_gauge_scalar_identity_error']:.3e} after accounting for the scalar.", "",
        "[Numerical] Two additional randomly rotated and gauged critical transverse-field Ising bonds use "
        "h=-J XX-J(ZI+IZ)/2, with Pauli matrices. Their exact A=J comes from c=1/2 and v=2J. "
        "These are an independent free critical class, with on-site field terms retained in the bond.", "",
        "| Bond | Exact A | Pair 8,10 | Bulk subtraction L=10 |", "| :-- | --: | --: | --: |"]
for m in r["models"]:
    if m["name"].startswith("random_Ising"):
        out.append(f"| {m['name']} | {m['exact_A']:.9f} | {m['diagnostic']['pairs'][-1]['A']:.9f} | {m['rows'][-1]['A_bulk_subtraction']:.9f} |")
out += ["", "## Opposite sign and extreme-bond checks", "",
        "[Proved] Hopping amplitudes t=exp(+/-i pi/3)/2 have A=1 each, whereas their midpoint has A=1/2. "
        "The nonzero straight segment stays critical. L=6,12 remove the endpoint phase twists. "
        "[Numerical] All six energies obey E=-A csc(pi/L) within "
        f"{r['checks']['free_hopping_ring_formula_max_error']:.3e}.", "",
        "[Proved] The exposed extreme bonds XX and YY with Pauli matrices have E=-L and A=0. "
        "Their midpoint has M=1 and A=2. [Numerical] Rings L=4,6,8 reproduce the classical endpoint "
        f"energies within {r['checks']['classical_ring_energy_max_error']:.3e}.", "",
        "## Validation", "",
        f"[Numerical] Script validation status: {r['validation_status']}. All error checks require <1e-9; "
        "the eigenpair residual check requires <1e-8. The analytic proof, not the grid, establishes the curvature sign.", "",
        "| Check | Value |", "| :-- | --: |"]
for name,value in r["checks"].items():
    out.append(f"| {name} | {value:.12g} |")
out += ["", "[Scope] No broad classification of random critical bonds, finite-volume Bethe-root solution, "
        "large spin-1 extrapolation, or proof for arbitrary projector bonds was completed or claimed. "
        "None is needed for the exact XXZ obstruction. Only files in work_C102 were modified.", ""]
(ROOT / "NUMERICS.md").write_text("\n".join(out))
print("Wrote NUMERICS.md")
