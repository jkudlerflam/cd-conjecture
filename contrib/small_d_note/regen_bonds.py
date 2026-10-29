"""Regenerate the d = 11 and d = 13 bonds of note.tex and compare them with the repository.

Run from anywhere with
    python3 -B regen_bonds.py
(-B keeps Python from writing __pycache__ into the repository folders).

The script
  (i)  rebuilds both bonds from the move tables written in note.tex (Tables 1 and 2),
       with no repository import,
  (ii) imports codex/work_C84/bond11.py and codex/work_C88/bond13.py read-only,
  (iii) compares the two entry by entry at several generic parameter points and checks
       Hermiticity, conservation of the charge q, the marker projectors and, for d = 13,
       occupied number,
  (iv) prints the off-diagonal moves of the rebuilt d = 11 bond in letter notation.
It writes nothing. Parameters are diagnostic, not the certified schedule.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.dont_write_bytecode = True

MAT = ['L0', 'L1', 'R0', 'R1', 'DA', 'DB', 'XA', 'XB', 'E', 'MA', 'MB']
VAC = ['V0', 'V1']
CHARGE = dict(L0=0, L1=2, R0=0, R1=2, DA=-1, XA=-1, MA=-1, DB=1, XB=1, MB=1, E=0, V0=0, V1=0)
MARKERS = ('E', 'MA', 'MB')


# ---------------------------------------------------------------- note.tex tables
def out_state(a):
    if a in VAC: return int(a[1])
    return 1 if a[0] in 'RDX' else 0

def in_set(a):
    if a in VAC: return (int(a[1]),)
    return {'R0': (0, 1), 'R1': (1,), 'E': (0, 1)}.get(a, (0,))

def allowed13(a, b):
    """Pinless grammar (d = 13): ab legal iff o(a) is in I(b)."""
    return out_state(a) in in_set(b)

def allowed11(a, b):
    """d = 11 grammar: the pinless grammar plus the boundary pin on L1 followed by a marker."""
    return allowed13(a, b) and not (a == 'L1' and b in MARKERS)

def material_table(P, lam):
    """Rows of Table 1 as (output pair, input pair, coefficient); h.c. added later."""
    z = lambda s, a: (-1) ** a if s == 'A' else 1
    phi = dict(A=P['phiA'], B=0.0)
    psi = dict(A=P['psiA'], B=np.pi / 2)
    rows = []
    for a in (0, 1):
        rows.append((('L%d' % a, 'R0'), ('R0', 'R%d' % a), -1.0))
        for s in 'AB':
            rows.append((('L%d' % a, 'D' + s), ('D' + s, 'R%d' % a), -z(s, a)))
            rows.append((('L%d' % a, 'X' + s), ('X' + s, 'R%d' % a), -np.exp(1j * phi[s]) * z(s, a)))
    rows.append((('L0', 'E'), ('R0', 'E'), -1.0))
    for s in 'AB':
        rows.append((('L0', 'M' + s), ('X' + s, 'E'), -np.exp(1j * phi[s])))
        rows.append((('L0', 'M' + s), ('D' + s, 'E'), -P['td'] * np.exp(1j * psi[s])))
        for m in MARKERS:
            rows.append(((m, 'X' + s), (m, 'D' + s), -1.0))
    rows.append((('DB', 'R0'), ('DA', 'R1'), P['g']))
    rows.append((('MA', 'DB'), ('E', 'R0'), lam * np.exp(-1j * np.pi / 4)))
    rows.append((('MB', 'DA'), ('E', 'R0'), lam * np.exp(+1j * np.pi / 4)))
    return rows

def onsite(a, P):
    v = 0.0
    if a in MARKERS: v += P['musel']
    if a[0] in 'DXM': v += P['Delta']
    if a[0] in 'XM': v += P['D']
    if a[0] in 'LR': v += P['field'] * int(a[1])
    return v

def my_b11(P, lam, pinned):
    d = 11; ix = {a: i for i, a in enumerate(MAT)}
    h = np.zeros((d * d, d * d), complex)
    for out, inn, c in material_table(P, lam):
        i = ix[out[0]] * d + ix[out[1]]; j = ix[inn[0]] * d + ix[inn[1]]
        h[i, j] += c; h[j, i] += np.conj(c)
    ok = allowed11 if pinned else allowed13
    for i, a in enumerate(MAT):
        for j, b in enumerate(MAT):
            h[d * i + j, d * i + j] += (onsite(a, P) + onsite(b, P)) / 2 + P['K'] * (not ok(a, b))
    return h

def my_h13(P, lam, J, t, mu):
    let = MAT + VAC; d = 13; ix = {a: i for i, a in enumerate(let)}
    h = np.zeros((d * d, d * d), complex)
    occ = [d * i + j for i in range(11) for j in range(11)]
    h[np.ix_(occ, occ)] = J * my_b11(P, lam, pinned=False)
    for a in MAT:  # -t (|a, V_o(a)><V_q, a| + h.c.), q in I(a)
        for q in in_set(a):
            i = ix[a] * d + ix['V%d' % out_state(a)]; j = ix['V%d' % q] * d + ix[a]
            h[i, j] += -t; h[j, i] += -t
    for i, a in enumerate(let):
        for j, b in enumerate(let):
            if a in VAC or b in VAC:
                h[d * i + j, d * i + j] += J * P['K'] * (not allowed13(a, b))
            h[d * i + j, d * i + j] -= mu * ((a not in VAC) + (b not in VAC)) / 2
    return h


# ---------------------------------------------------------------- repository imports
def import_from(folder, module):
    path = os.path.join(REPO, 'codex', folder)
    sys.path.insert(0, path)
    try:
        for m in ('bond18', 'parent_shift', 'bond11', 'bond13'):
            sys.modules.pop(m, None)
        return __import__(module)
    finally:
        sys.path.remove(path)

def repo_points(p):
    """Diagnostic points; psi_A = p phi_A + pi/2 as in the repository code."""
    rng = np.random.default_rng(7 + p)
    pts = []
    for flux in (np.pi / (2 * p), 0.37):
        pts.append(dict(n=p, D=rng.uniform(3, 9), Delta=rng.uniform(1, 3), g=rng.uniform(.01, .2),
                        td=rng.uniform(1e-3, .1), K=rng.uniform(20, 90), musel=rng.uniform(.1, 2),
                        field=rng.uniform(.01, .3), flux=flux, theta=-np.pi / 4))
    return pts

def to_mine(q):
    return dict(phiA=q['flux'], psiA=q['n'] * q['flux'] + np.pi / 2, td=q['td'], g=q['g'],
                D=q['D'], Delta=q['Delta'], musel=q['musel'], field=q['field'], K=q['K'])


def main():
    b11 = import_from('work_C84', 'bond11')
    b13 = import_from('work_C88', 'bond13')
    q11 = np.diag([CHARGE[a] for a in MAT]).astype(float)
    worst = 0.0
    for p in (2, 3, 12, 14):
        for q in repo_points(p):
            for lam in (1.0, 0.31):
                ref = b11.bond(q, lam=lam); mine = my_b11(to_mine(q), lam, pinned=True)
                worst = max(worst, abs(ref - mine).max())
                Q = np.kron(q11, np.eye(11)) + np.kron(np.eye(11), q11)
                assert abs(mine - mine.conj().T).max() == 0
                assert abs(mine @ Q - Q @ mine).max() < 1e-12
                pm = np.diag([float(a in MARKERS) for a in MAT])  # on-site marker projector
                for Pm in (np.kron(pm, np.eye(11)), np.kron(np.eye(11), pm)):
                    assert abs(mine @ Pm - Pm @ mine).max() < 1e-12
    print('d = 11: max |bond11.bond - note Table 1| =', worst)
    q = repo_points(12)[1]
    print('d = 11: nonzero entries at a generic point =', np.count_nonzero(np.abs(my_b11(to_mine(q), 0.31, True)) > 0))

    let = MAT + VAC
    q13 = np.diag([CHARGE[a] for a in let]).astype(float)
    nocc = np.diag([float(a not in VAC) for a in let])
    worst = 0.0
    for p in (2, 3, 14):
        for q in repo_points(p):
            for (lam, J, t, mu) in ((0.1, 1.0, 1.0, 0.0), (0.37, 0.7, 1.3, 0.21)):
                ref = b13.bond(b13.params(**q), lam=lam, J=J, t=t, mu=mu)
                mine = my_h13(to_mine(q), lam, J, t, mu)
                worst = max(worst, abs(ref - mine).max())
                assert abs(mine - mine.conj().T).max() == 0
                for X in (q13, nocc, np.diag([float(a in MARKERS) for a in let])):
                    XX = np.kron(X, np.eye(13)) + np.kron(np.eye(13), X)
                    assert abs(mine @ XX - XX @ mine).max() < 1e-12
    print('d = 13: max |bond13.bond - note eq. (h13)| =', worst)

    # Pin identity at the bond level: h11(pinned) = h11(pinless) + K Pi_(L1, marker)
    P = to_mine(repo_points(14)[0])
    pin = np.diag([float(a == 'L1' and b in MARKERS) for a in MAT for b in MAT])
    print('pin identity residual =', abs(my_b11(P, .2, True) - my_b11(P, .2, False) - P['K'] * pin).max())

    print('\nRows of Table 1 at phi_A = pi/24, psi_A = pi, t_d = 0.003, g = 0.07, lambda = 1 (each also with h.c.):')
    P = dict(phiA=np.pi / 24, psiA=np.pi, td=0.003, g=0.07, D=0., Delta=0., musel=0., field=0., K=0.)
    h = my_b11(P, 1.0, True); d = 11; ix = {a: i for i, a in enumerate(MAT)}
    offdiag = int(np.count_nonzero(np.abs(h - np.diag(np.diag(h))) > 0))
    for out, inn, c in material_table(P, 1.0):
        i = ix[out[0]] * d + ix[out[1]]; j = ix[inn[0]] * d + ix[inn[1]]
        print('  %-3s %-3s <- %-3s %-3s : %s' % (out + inn + (np.round(h[i, j], 6),)))
    print('  off-diagonal nonzero entries:', offdiag, '(twice the number of rows)')

if __name__ == '__main__':
    main()
