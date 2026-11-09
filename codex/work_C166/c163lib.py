"""C163 library (Claude Opus agent, 2026-11-08): bond constructions for the attractive-window path scans and the
twisted curve J = -R(phi).  Conventions of est.py / tk.py: H = sum_j [s_j^T J s_{j+1} + f . s_j], Pauli matrices,
v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz) after a common rotation that diagonalizes the symmetric part of J."""
import numpy as np


def Rz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])


def R_axis(m, t):
    m = np.asarray(m, float); m = m / np.linalg.norm(m)
    K = np.array([[0, -m[2], m[1]], [m[2], 0, -m[0]], [-m[1], m[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


def J_to_v(J, f=np.zeros(3)):
    """General real J (3x3) and field f -> v (common rotation diagonalizing (J + J^T)/2)."""
    J = np.asarray(J, float); f = np.asarray(f, float)
    w, R = np.linalg.eigh((J + J.T) / 2)
    Jp = R.T @ J @ R; fp = R.T @ f
    Dz = (Jp[0, 1] - Jp[1, 0]) / 2; Dx = (Jp[1, 2] - Jp[2, 1]) / 2; Dy = (Jp[2, 0] - Jp[0, 2]) / 2
    return np.array([w[0], w[1], w[2], Dx, Dy, Dz, fp[0], fp[1], fp[2]])


# non-axial symmetric directions about z (U(1)_z charge 2: xx-yy, xy; charge 1: xz, yz)
def _E(a, b):
    M = np.zeros((3, 3)); M[a, b] += 1; M[b, a] += 1; return M / np.linalg.norm(M)


HARM = [(np.diag([1., -1, 0]) / np.sqrt(2)), _E(0, 1), _E(0, 2), _E(1, 2)]


def path_bond(spec, lam):
    """spec: dict(mode, Delta, sperp, theta, H (4 weights), hm, b (3 field weights)).
    T : J = (S0 + lam |S0| H) Rz(theta), f = 0.
    KI: J = (S0 + lam |S0| mu (xx - yy)/sqrt2) Rz(theta), f = lam |S0| b_inplane  (parent at zero field, sheet K_I).
    N : J = (S0 + lam |S0| H) Rz(theta), f = hm e_z + lam |S0| b  (parent in a field along m = z)."""
    sp, De = spec['sperp'], spec['Delta']
    S0 = np.diag([sp, sp, De * abs(sp)])
    n0 = np.linalg.norm(S0)
    H = sum(w * M for w, M in zip(spec['H'], HARM))
    f = np.zeros(3)
    if spec['mode'] == 'KI':
        H = spec['H'][0] * HARM[0]
        f = lam * n0 * np.array([spec['b'][0], spec['b'][1], 0.0])
    elif spec['mode'] == 'N':
        f = spec['hm'] * np.array([0, 0, 1.0]) + lam * n0 * np.asarray(spec['b'])
    J = (S0 + lam * n0 * H) @ Rz(spec['theta'])
    return J, f


def make_specs(seed=163):
    rng = np.random.default_rng(seed)
    specs = []

    def theta():
        while True:
            t = rng.uniform(0.25, np.pi - 0.25)
            if abs(t - np.pi / 2) > 0.15 and abs(t - np.pi / 3) > 0.1 and abs(t - 2 * np.pi / 3) > 0.1:
                return t
    for De in (-0.95, -0.85, -0.75, -0.67):
        for k in range(4):
            Hw = rng.normal(size=4); Hw /= np.linalg.norm(Hw)
            specs.append(dict(mode='T', Delta=De, sperp=float(rng.choice([-1, 1])), theta=theta(), H=Hw.tolist(), hm=0.0, b=[0, 0, 0]))
    for De in (-0.95, -0.85, -0.75, -0.67):
        for k in range(3):
            w = rng.normal(size=3); w /= np.linalg.norm(w)   # (mu, bx, by)
            specs.append(dict(mode='KI', Delta=De, sperp=float(rng.choice([-1, 1])), theta=theta(), H=[w[0], 0, 0, 0], hm=0.0, b=[w[1], w[2], 0]))
    for De in (-0.85, -0.7, -0.5, -0.3):
        for k in range(3):
            Hw = rng.normal(size=4); b = rng.normal(size=3)
            n = np.sqrt(np.sum(Hw ** 2) + np.sum(b ** 2)); Hw /= n; b /= n
            hsat = 2 * (1 + De)
            specs.append(dict(mode='N', Delta=De, sperp=1.0, theta=theta(), H=Hw.tolist(), hm=float(rng.uniform(0.15, 0.85) * hsat), b=b.tolist()))
    for i, s in enumerate(specs):
        s['tag'] = f"{s['mode']}{i:02d}"
    return specs


def twisted_specs(seed=1631):
    """J = -R_z(phi) + eps X, f = eps b (b = 0 for type T)."""
    rng = np.random.default_rng(seed)
    out = []
    for phi in (1.0, 2.0, 2 * np.pi / 3, np.pi / 2):
        for eps in (0.15, 0.3):
            for typ in ('T', 'N'):
                X = rng.normal(size=(3, 3)); X *= np.sqrt(3) / np.linalg.norm(X)
                b = rng.normal(size=3) if typ == 'N' else np.zeros(3)
                if typ == 'N':
                    b *= 0.5 / np.linalg.norm(b)
                out.append(dict(phi=phi, eps=eps, typ=typ, X=X.tolist(), b=b.tolist(), tag=f'G_{phi:.3f}_{eps:.2f}_{typ}'))
    return out


def twisted_bond(s):
    J = -R_axis([0, 0, 1], s['phi']) + s['eps'] * np.array(s['X'])
    f = s['eps'] * np.array(s['b'])
    return J, f


if __name__ == '__main__':
    import est as E
    rng = np.random.default_rng(0)
    # J_to_v round trip: spectra of 8-site rings agree
    for _ in range(3):
        J = rng.normal(size=(3, 3)); f = rng.normal(size=3)
        v = J_to_v(J, f)
        J2, f2 = E.J_f(v)
        h1 = E.bond_from_JB(J, f); h2 = E.bond_from_JB(J2, f2)
        w1 = np.sort(np.concatenate([np.linalg.eigvalsh(E.table(8).sector_matrix(h1, m).toarray()) for m in range(8)]))
        w2 = np.sort(np.concatenate([np.linalg.eigvalsh(E.table(8).sector_matrix(h2, m).toarray()) for m in range(8)]))
        print('J_to_v ring-8 spectrum diff', np.max(np.abs(w1 - w2)))
    import tk
    for s in make_specs()[:1] + make_specs()[16:17] + make_specs()[28:29]:
        for lam in (0.0, 0.5):
            J, f = path_bond(s, lam)
            p = tk.spiral_parent(J, f)
            print(s['tag'], lam, 'parent Delta', round(p['Delta'], 3), 'theta', round(p['theta'], 3), 'harm', round(p['harm'], 3), 'nsym', len(tk.symmetry_check(J, f)))
    # twisted curve is h_F after the rotation: open 6-site chain spectrum
    J = -R_axis([0, 0, 1], 1.0)
    h = E.bond_from_JB(J, np.zeros(3)); hF = E.bond_from_JB(-np.eye(3), np.zeros(3))
    import itertools
    def openH(h, L):
        H = 0
        for j in range(L - 1):
            H = H + np.kron(np.kron(np.eye(2 ** j), h), np.eye(2 ** (L - j - 2)))
        return H
    print('twisted curve vs h_F, open L=6:', np.max(np.abs(np.linalg.eigvalsh(openH(h, 6)) - np.linalg.eigvalsh(openH(hF, 6)))))
