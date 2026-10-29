"""Small-ring C102 checks. All generated files stay beside this script.

Run with python3 checks.py. Only numpy, scipy, and the Python standard library.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import json
import time
import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parent
RNG = np.random.default_rng(102)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])
I = np.eye(2)
PAULI = [X, Y, Z]
XYZ = [np.kron(s, s).real / 4 for s in PAULI]


def xxz(delta, J=1):
    return J * (XYZ[0] + XYZ[1] + delta * XYZ[2])


def ba_A(delta, J=1):
    gamma = np.arccos(delta)
    return J * np.pi / 2 * np.sinc(gamma / np.pi)


def ba_e(delta):
    gamma = np.arccos(delta)
    def integrand(x):
        if x == 0:
            return (np.pi - gamma) / np.pi
        return (2 * np.exp(-2 * gamma * x)
                * (-np.expm1(-2 * (np.pi - gamma) * x))
                / ((-np.expm1(-2 * np.pi * x))
                   * (1 + np.exp(-2 * gamma * x))))
    value, error = quad(integrand, 0, np.inf, epsabs=2e-13, epsrel=2e-13)
    return delta / 4 - np.sin(gamma) * value, error


def ring(h, d, L, sector=False):
    """Sparse periodic ring in full space or total Sz=0.

    SU(2) spin-1 multiplets always contain Sz=0, so that sector contains
    a global ground state. XXZ uses its known zero-field ground sector.
    """
    all_states = np.arange(d ** L, dtype=np.int64)
    powers = d ** np.arange(L, dtype=np.int64)
    digits = (all_states[:, None] // powers) % d
    if sector:
        keep = (2 * digits.sum(axis=1) == L * (d - 1))
        states, digits = all_states[keep], digits[keep]
    else:
        states = all_states
    lookup = np.full(d ** L, -1, dtype=np.int64)
    lookup[states] = np.arange(len(states))
    rows, cols, vals = [], [], []
    for i in range(L):
        j = (i + 1) % L
        local = d * digits[:, i] + digits[:, j]
        for a, b in zip(*np.nonzero(np.abs(h) > 1e-13)):
            src = np.flatnonzero(local == b)
            dst = states[src] + (a // d - b // d) * powers[i] + (a % d - b % d) * powers[j]
            target = lookup[dst]
            assert np.all(target >= 0), "bond does not preserve selected sector"
            rows.append(target)
            cols.append(src)
            vals.append(np.full(len(src), h[a, b]))
    mat = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                     shape=(len(states), len(states))).tocsr()
    mat.eliminate_zeros()
    if np.max(np.abs(mat.data.imag)) < 1e-13:
        mat = mat.real
    return mat


def ground(h, d, L, sector=False):
    mat = ring(h, d, L, sector)
    v0 = RNG.normal(size=mat.shape[0])
    if np.iscomplexobj(mat.data):
        v0 = v0 + 1j * RNG.normal(size=len(v0))
    if len(v0) < 20:
        vals, vecs = np.linalg.eigh(mat.toarray())
        E, vec = vals[0], vecs[:, 0]
    else:
        vals, vecs = eigsh(mat, k=1, which="SA", tol=2e-12, v0=v0, maxiter=5000)
        E, vec = vals[0], vecs[:, 0]
    return {"L": L, "E": float(E), "dimension": mat.shape[0],
            "residual": float(np.linalg.norm(mat @ vec - E * vec))}


def diagnostic(rows):
    L = np.array([r["L"] for r in rows], float)
    e = np.array([r["E"] for r in rows]) / L
    pair = []
    for i in range(len(L) - 1):
        value = 6 / np.pi * (e[i + 1] - e[i]) / (L[i] ** -2 - L[i + 1] ** -2)
        pair.append({"lengths": [int(L[i]), int(L[i+1])], "A": float(value)})
    fit = np.polynomial.polynomial.polyfit(L ** -2, e, 2)
    return {"pairs": pair, "quadratic_in_inverse_L_squared": float(-6 / np.pi * fit[1])}


def measure(name, h, d, lengths, sector=False, exact_A=None, exact_e=None):
    rows = [ground(h, d, L, sector) for L in lengths]
    ev = np.linalg.eigvalsh(h)
    out = {"name": name, "d": d, "raw_half_width": float((ev[-1]-ev[0])/2),
           "rows": rows, "diagnostic": diagnostic(rows), "exact_A": exact_A}
    if exact_e is not None:
        out["exact_e"] = exact_e
        for r in rows:
            r["A_bulk_subtraction"] = float(6 / np.pi * r["L"] * (r["L"] * exact_e - r["E"]))
    print(name, "last pair", out["diagnostic"]["pairs"][-1]["A"], "exact", exact_A, flush=True)
    return out


def transformed(h, unitary, scalar, o):
    u2 = np.kron(unitary, unitary)
    return u2 @ h @ u2.conj().T + scalar * np.eye(4) + np.kron(o, I) - np.kron(I, o)


def main():
    t0 = time.time()
    result = {"seed": 102, "models": [], "random_pairs": [], "checks": {}}
    for delta in [-.75, -.5, 0., .5, .75]:
        e, error = ba_e(delta)
        row = measure(f"XXZ_{delta}", xxz(delta), 2, [6, 8, 10, 12, 14, 16], True,
                      float(ba_A(delta)), e)
        row["quadrature_error"] = error
        result["models"].append(row)
    triples = [(1, 1, 0), (1, 0, 1), (1, .5, .5), (1.5, .5, 0), (.5, 1.5, 0)]
    for js in triples:
        h = sum(j * b for j, b in zip(js, XYZ))
        exact = 1. if js in triples[:2] else 0.
        row = measure("XYZ_" + str(js), h, 2, [6, 8, 10, 12], exact_A=exact)
        result["models"].append(row)
    # Spin-1 bilinear-biquadratic family in literal h=a X+b X^2 units.
    sz = np.diag([1., 0., -1.])
    sp = np.array([[0., np.sqrt(2), 0.], [0., 0., np.sqrt(2)], [0., 0., 0.]])
    sx, sy = (sp + sp.T)/2, (sp-sp.T)/(2j)
    dot = sum(np.kron(s, s) for s in [sx, sy, sz]).real
    for b in [-1., -0.5, 0., 1/3, .5]:
        h = dot + b * (dot @ dot)
        exact_A = float(3*np.pi) if b == -1 else 0.
        exact_e = -4. if b == -1 else (-2/3 if b == 1/3 else None)
        result["models"].append(measure(f"BB_b={b}", h, 3, [4, 6, 8, 10], True,
                                        exact_A, exact_e))
    for b in [1., 1.5, 2.]:
        h = dot + b * (dot @ dot)
        exact_A = float(4*np.pi/3) if b == 1 else None
        exact_e = float(2-np.log(3)-np.pi/(3*np.sqrt(3))) if b == 1 else None
        result["models"].append(measure(f"BB_b={b}", h, 3, [3, 6, 9], True,
                                        exact_A, exact_e))
    # Same random rotation for each triple keeps its full segment in a critical cone.
    for p in range(4):
        axis = RNG.normal(size=3)
        U = expm(1j * sum(a*s for a,s in zip(axis, PAULI)))
        js = RNG.uniform(.6, 1.5, size=2)
        ds = np.array([RNG.uniform(-.85, -.2), RNG.uniform(.2, .85)])
        hs, scalars = [], []
        for J, delta in zip(js, ds):
            scalar = float(RNG.uniform(-.3, .3))
            o = sum(a*s for a,s in zip(RNG.normal(size=3)*.2, PAULI))
            hs.append(transformed(xxz(delta, J), U, scalar, o))
            scalars.append(scalar)
        hs.append((hs[0]+hs[1])/2)
        jmid = js.mean()
        dmid = np.dot(js, ds) / js.sum()
        ex = [float(ba_A(ds[k], js[k])) for k in range(2)] + [float(ba_A(dmid, jmid))]
        ms = [float(js[0]/2), float(js[1]/2), float(jmid/2)]
        runs = []
        for k,h in enumerate(hs):
            J = js[k] if k<2 else jmid
            delta = ds[k] if k<2 else dmid
            scalar = scalars[k] if k<2 else np.mean(scalars)
            e = J*ba_e(delta)[0] + scalar
            runs.append(measure(f"random_XXZ_{p}_{k}",h,2,[6,8,10],False,ex[k],e))
        defect = ex[2]-(ex[0]+ex[1])/2
        ed_defect = runs[2]["rows"][-1]["A_bulk_subtraction"] - sum(r["rows"][-1]["A_bulk_subtraction"] for r in runs[:2])/2
        invariance_errors = [abs(runs[k]["rows"][0]["E"] -
                            (js[k]*ground(xxz(ds[k]),2,6,True)["E"] + 6*scalars[k]))
                             for k in range(2)]
        result["random_pairs"].append({"J":js.tolist(), "delta":ds.tolist(), "M":ms,
                 "runs": runs, "exact_midpoint_defect":defect, "ED_midpoint_defect_L10":ed_defect,
                 "unitary_gauge_scalar_identity_errors_L6":invariance_errors})
    # Independent free critical class, with a field in the bond and complex rotation.
    for p in range(2):
        J = float(RNG.uniform(.5,1.5))
        h = -J*np.kron(X,X)-J/2*(np.kron(Z,I)+np.kron(I,Z))
        U = expm(1j*sum(a*s for a,s in zip(RNG.normal(size=3),PAULI)))
        o = sum(a*s for a,s in zip(RNG.normal(size=3)*.2,PAULI))
        h = transformed(h,U,0.,o)
        result["models"].append(measure(f"random_Ising_{p}",h,2,[6,8,10],False,J,-4*J/np.pi))
    # Independent exact identities and geometry checks.
    checks = result["checks"]
    checks["XX_bulk_quadrature_error"] = abs(ba_e(0)[0]+1/np.pi)
    xx = next(r for r in result["models"] if r["name"] == "XXZ_0.0")
    checks["XX_ring_formula_max_error"] = max(abs(r["E"]+1/np.sin(np.pi/r["L"])) for r in xx["rows"])
    aklt = next(r for r in result["models"] if r["name"] == "BB_b=0.3333333333333333")
    checks["AKLT_exact_energy_max_error"] = max(abs(r["E"]+2*r["L"]/3) for r in aklt["rows"])
    # On a small ring, compare full and reduced bases.
    checks["XXZ_sector_full_error_L6"] = abs(ground(xxz(.37),2,6,True)["E"]-ground(xxz(.37),2,6,False)["E"])
    checks["spin1_sector_full_error_L4"] = abs(ground(dot- dot@dot,3,4,True)["E"]-ground(dot-dot@dot,3,4,False)["E"])
    gamma = np.linspace(.01,np.pi-.01,2001)
    curv = (2-gamma/np.tan(gamma)-gamma**2/np.sin(gamma)**2)/(gamma**3*np.sin(gamma))
    checks["largest_sampled_XXZ_f_second_derivative"] = float(curv.max())
    # Finite-volume Jensen identity: D_A=(6/pi)L(L D_e-D_E).
    m = {r["name"]:r for r in result["models"]}
    left,mid,right = [m[f"XXZ_{d}"] for d in [-.5,0.,.5]]
    defects=[]
    for a,b,c in zip(left["rows"],mid["rows"],right["rows"]):
        L=b["L"]
        de=mid["exact_e"]-(left["exact_e"]+right["exact_e"])/2
        dE=b["E"]-(a["E"]+c["E"])/2
        defects.append({"L":L,"bulk_Jensen_defect":de,"ring_Jensen_defect":dE,
                        "A_midpoint_defect":6/np.pi*L*(L*de-dE)})
    result["XXZ_Jensen_identity"] = defects
    # Noncollinear free hopping gives the opposite Jensen sign, also at d=2.
    # The phase is removable on both even rings because exp(i L pi/3)=1.
    hopping = []
    hp = np.zeros((4,4),complex)
    hp[1,2] = np.exp(1j*np.pi/3)/2
    hp[2,1] = hp[1,2].conjugate()
    hm = hp.conjugate()
    for name,h in [("plus",hp),("minus",hm),("mid",(hp+hm)/2)]:
        ts = abs(h[1,2])
        rows = [ground(h,2,L,True) for L in [6,12]]
        for r in rows:
            r["A_bulk_subtraction"] = 6/np.pi*r["L"]*(-2*ts/np.pi*r["L"]-r["E"])
        hopping.append({"name":name,"A":float(2*ts),"M":float(ts),"rows":rows})
    result["free_hopping_opposite_Jensen"] = hopping
    checks["free_hopping_ring_formula_max_error"] = max(abs(r["E"]+m["A"]/np.sin(np.pi/r["L"])) for m in hopping for r in m["rows"])
    # Both half-width-one commuting bonds are exposed quotient-ball points.
    # Their midpoint is critical XX with A=2 and M=1.
    classical=[]
    for name,h in [("XX_pauli",np.kron(X,X)),("YY_pauli",np.kron(Y,Y)),
                   ("midpoint",(np.kron(X,X)+np.kron(Y,Y))/2)]:
        classical.append({"name":name,"rows":[ground(h,2,L) for L in [4,6,8]]})
    result["extreme_bond_screen"] = classical
    checks["classical_ring_energy_max_error"] = max(abs(r["E"]+r["L"]) for m in classical[:2] for r in m["rows"])
    result["spin1_critical_midpoint_defects"] = {}
    b1,bm,b2=[m[f"BB_b={b}"] for b in [1.,1.5,2.]]
    for key in ["quadratic_in_inverse_L_squared"]:
        result["spin1_critical_midpoint_defects"][key]=bm["diagnostic"][key]-(b1["diagnostic"][key]+b2["diagnostic"][key])/2
    result["spin1_critical_midpoint_defects"]["pair_L6_L9"] = bm["diagnostic"]["pairs"][-1]["A"]-(b1["diagnostic"]["pairs"][-1]["A"]+b2["diagnostic"]["pairs"][-1]["A"])/2
    all_runs = result["models"] + [run for pair in result["random_pairs"] for run in pair["runs"]] + hopping + classical
    checks["max_eigenpair_residual"] = max(r["residual"] for run in all_runs for r in run["rows"])
    checks["max_unitary_gauge_scalar_identity_error"] = max(e for pair in result["random_pairs"] for e in pair["unitary_gauge_scalar_identity_errors_L6"])
    result["number_recorded_eigenpairs"] = sum(len(run["rows"]) for run in all_runs)
    result["max_matrix_dimension"] = max(r["dimension"] for run in all_runs for r in run["rows"])
    for name,value in checks.items():
        if "error" in name:
            assert value < 1e-9, (name,value)
    assert checks["max_eigenpair_residual"] < 1e-8
    assert checks["largest_sampled_XXZ_f_second_derivative"] < 0
    assert all(x["exact_midpoint_defect"]>0 and x["ED_midpoint_defect_L10"]>0 for x in result["random_pairs"])
    result["validation_status"] = "passed"
    result["elapsed_seconds"] = time.time()-t0
    (ROOT/"numerics.json").write_text(json.dumps(result,indent=2)+"\n")
    print("checks",checks,flush=True)
    print("elapsed_seconds",result["elapsed_seconds"],flush=True)


if __name__ == "__main__":
    main()
