"""Independent small-state entropy and analytic-calibration checks."""
from pathlib import Path
import importlib.util
import json
from math import factorial
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
NUM = HERE.parent / "numerics"
spec = importlib.util.spec_from_file_location("vacancy_control", NUM/"xx_vacancy.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

L, up, down, J = 8, 2, 2, .1
codes, mats = module.build(L, up, down)
dim = factorial(L)//(factorial(up)*factorial(down)*factorial(L-up-down))
assert dim == len(codes) == len(set(codes))
assert mats[0].nnz == dim*2*(up+down)*(L-up-down)//(L-1)
assert mats[1].nnz == dim*2*up*down//(L-1)
assert all(np.all(m.data == 1) for m in mats)
energy, states, residual = module.ground(-mats[0]-J*mats[1])
state = states[:, 0]
full = np.zeros(3**L)
full[codes] = state
saved = json.loads((NUM/f"result_L{L}.json").read_text())
record = next(r for r in saved["runs"] if r["J"] == J)
errors = []
for ell in range(1, L):
    # A literal tensor reshape and complete Schmidt decomposition,
    # independent of the worker's block-by-block sector tracing.
    schmidt = np.linalg.svd(full.reshape(3**(L-ell), 3**ell), compute_uv=False)
    p = schmidt**2
    p = p[p > 1e-28]
    entropy = -p@np.log(p)
    block_entropy = module.entropy(codes, state, L, ell, up, down)
    errors.append(abs(entropy-block_entropy))
    if ell <= L//2:
        assert abs(entropy-record["entropy"][str(ell)]) < 1e-10
assert max(errors) < 1e-12
assert abs(record["entropy"]["1"]-1.5*np.log(2)) < 1e-12

size_fits = []
for J in [.2, .1, .05]:
    values = []
    for size in [8, 12, 16]:
        data = json.loads((NUM/f"result_L{size}.json").read_text())
        row = next(r for r in data["runs"] if r["J"] == J)
        values.append(row["entropy"][str(size//2)])
    fit = np.polyfit(np.log([8, 12, 16]), values, 1)
    c0816 = 3*(values[-1]-values[0])/np.log(2)
    c1216 = 3*(values[-1]-values[1])/np.log(16/12)
    size_fits.append({"J": J, "c_three_size": 3*fit[0],
                      "c_L8_L16": c0816, "c_L12_L16": c1216})

out = {"small_sector_dimension": dim,
       "energy_error_vs_saved": abs(float(energy[0])-record["energies"][0]),
       "max_full_tensor_entropy_error": float(max(errors)),
       "small_sector_ground_residual": residual[0],
       "half_chain_fits": size_fits,
       "casimir_slope_first_order": float(2-12/np.pi**2),
       "casimir_scope": "v_charge+v_spin = 2 + J*(2-12/pi^2) + O(J^2), t=1"}
(HERE/"numerics_audit.json").write_text(json.dumps(out, indent=2)+"\n")
print(json.dumps(out, indent=2))
