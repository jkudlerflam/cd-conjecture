"""Check what I P K forbids in the full transverse-DM qubit chain.

Independent literal Pauli-bond construction and periodic-ring ED. Small
rings corroborate symmetry algebra only; they do not identify a CFT.
"""
from pathlib import Path
import json
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

I = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1])


def bond(p):
    kron = np.kron
    return (
        p["Jx"] * kron(X, X) + p["Jy"] * kron(Y, Y)
        + p["D"] * (kron(X, Y) - kron(Y, X))
        + p["Jz"] * kron(Z, Z)
        + p["hz"] * (kron(Z, I) + kron(I, Z)) / 2
        + p["Gamma"] * (kron(X, Z) - kron(Z, X))
    )


def ring(L, p):
    h = bond(p)
    rr, cc, vv = [], [], []
    for s in range(1 << L):
        for j in range(L):
            k = (j + 1) % L
            a, b = (s >> j) & 1, (s >> k) & 1
            col = 2 * a + b
            for row in np.flatnonzero(abs(h[:, col]) > 0):
                ar, br = divmod(int(row), 2)
                f = s ^ ((a ^ ar) << j) ^ ((b ^ br) << k)
                rr.append(f); cc.append(s); vv.append(h[row, col])
    return coo_matrix((vv, (rr, cc)), shape=(1 << L, 1 << L)).tocsr()


def measure(psi, L):
    states = np.arange(1 << L)
    x = y = 0j
    for j in range(L):
        flipped = states ^ (1 << j)
        x += np.vdot(psi[flipped], psi)
        y += np.vdot(psi[flipped], 1j * (1 - 2 * ((states >> j) & 1)) * psi)
    n = sum(s.bit_count() * abs(psi[s]) ** 2 for s in states)
    return dict(X=float(x.real / L), Y=float(y.real / L), density=float(n / L))


def main():
    p = dict(Jx=0.4979983, Jy=-0.5020017, D=0.0152122,
             Jz=4.0, hz=-7.9463, Gamma=0.13693064)
    rows = []
    for L in (6, 8, 10):
        h = ring(L, p)
        states = np.arange(1 << L)
        parity = np.array([(-1) ** int(s).bit_count() for s in states])
        inv = np.array([sum(((int(s) >> j) & 1) << ((-j) % L) for j in range(L)) for s in states])
        transformed = h.conjugate()[inv][:, inv].multiply(parity[:, None] * parity[None, :])
        residual = float(abs(h - transformed).max())
        assert residual < 1e-12
        e, vec = eigsh(h, k=2, which="SA", tol=1e-12, maxiter=20000,
                      v0=np.random.default_rng(L).normal(size=1 << L))
        order = np.argsort(e); e = e[order]; vec = vec[:, order]
        r = dict(L=L, energy=float(e[0]), gap=float(e[1]-e[0]),
                 symmetry_residual=residual, **measure(vec[:, 0], L))
        assert abs(r["X"]) < 1e-8
        rows.append(r)
    out = dict(
        scope="full-chain small-ring symmetry check; no central-charge claim",
        parameters=p,
        exact_symmetry_action={"X": "odd", "Y": "even", "Z": "even"},
        rows=rows,
    )
    text=json.dumps(out,indent=2)+"\n"
    Path(__file__).with_name("transverse_symmetry.json").write_text(text)
    print(text,end="")


if __name__ == "__main__":
    main()
