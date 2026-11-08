"""Central charge on the critical line of the Kitaev spin chain in a [110] field, uniform gauge, PBC ED."""
import sys, numpy as np, scipy.sparse.linalg as spla
from ed import build, ent
K = 1.0
Ju = np.zeros((3, 3)); Ju[0, 1] = K
for h in [float(x) for x in sys.argv[1].split(',')]:
    hv = -h * np.array([1, 1, 0]) / np.sqrt(2)
    res = []
    for L in [int(x) for x in sys.argv[2].split(',')]:
        H = build(lambda j: Ju, lambda j: hv, L)
        w, v = spla.eigsh(H, k=4, which='SA', tol=1e-10)
        o = np.argsort(w); w = w[o]; v = v[:, o]
        psi = v[:, 0]
        Ss = [ent(psi, L, l) for l in range(1, L // 2 + 1)]
        res.append((L, w[0] / L, w[1] - w[0], w[2] - w[0], Ss[-1], Ss))
        # chord fit S = c/3 log(L/pi sin(pi l/L)) + a for l=2..L/2
        ls = np.arange(2, L // 2 + 1); X = np.log(L / np.pi * np.sin(np.pi * ls / L))
        c = 3 * np.polyfit(X, np.array(Ss[1:]), 1)[0]
        print(f"h={h} L={L} e0={w[0]/L:.8f} gaps={w[1]-w[0]:.5f},{w[2]-w[0]:.5f},{w[3]-w[0]:.5f} S(L/2)={Ss[-1]:.5f} c_chord={c:.3f}", flush=True)
    for a, b in zip(res[:-1], res[1:]):
        print(f"   c from S(L/2) pair {a[0]}->{b[0]}: {3*(b[4]-a[4])/np.log(b[0]/a[0]):.3f}")
