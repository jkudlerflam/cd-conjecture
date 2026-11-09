"""Full-qubit finite-ring onset diagnostics; not a CFT fit."""
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse.linalg import eigsh

from compare_full import choices, bond, est, Y, I


def scan(name, L, fields):
    table = est.table(L)
    c = choices()[name]
    ybond = (np.kron(Y, I)+np.kron(I, Y))/2
    ymat = {m: table.sector_matrix(ybond, m) for m in range(L)}
    particles = {m: np.array([int(s).bit_count() for s in table.sector(m)["R"]]) for m in range(L)}
    previous = {}
    output = []
    for hz in fields:
        h = bond(dict(c, hz=hz))
        levels = []
        best = None
        for m in range(L):
            mat = table.sector_matrix(h, m)
            e, v = eigsh(mat, k=2, which="SA", tol=2e-11, v0=previous.get(m), ncv=28)
            order = np.argsort(e)
            e, v = e[order], v[:, order]
            previous[m] = v[:, 0]
            levels.extend(e.tolist())
            if best is None or e[0] < best[0]:
                best = (e[0], m, v[:, 0], float(np.linalg.norm(mat@v[:, 0]-e[0]*v[:, 0])))
        energy, m, v, residual = best
        n = particles[m]
        prob = np.abs(v)**2
        density = float(prob@n/L)
        y = float(np.real(np.vdot(v, ymat[m]@v))/L)
        full = table.full_vector(v, m)
        singular = np.linalg.svd(full.reshape(1<<(L//2), -1), compute_uv=False)
        weights = singular**2
        weights = weights[weights > 1e-15]
        entropy = float(-weights@np.log(weights))
        levels.sort()
        row = dict(name=name, L=L, hz=hz, energy=float(energy), momentum_index=m,
                   density=density, Y=y, half_ring_entropy=entropy,
                   gap=float(levels[1]-levels[0]), residual=residual,
                   particle_variance=float(prob@(n*n)-(prob@n)**2))
        output.append(row)
        print(json.dumps(row), flush=True)
    return output


if __name__ == "__main__":
    start = time.monotonic()
    rows = []
    for name in ("C183", "binary"):
        rows += scan(name, 14, (-7.948, -7.9475, -7.94725, -7.947, -7.94675, -7.9465, -7.9463))
    rows += scan("binary", 16, (-7.9475, -7.947, -7.9465))
    result = dict(scope="finite full-qubit rings; no central-charge estimate", rows=rows, seconds=time.monotonic()-start)
    target = Path(__file__).with_name("onset_scan.json")
    target.write_text(json.dumps(result, indent=2)+"\n")
    print(target)
