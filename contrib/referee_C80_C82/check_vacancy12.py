"""(3) Vacancy lifts, zero-excitation sector, zero port, written by the referee.

Models (all on rings, vacancy hopping -t on every bond, material moves -J and
penalty J*K on adjacent occupied pairs, as in the referee_C75 F5 test):
  u15 : C75 15-letter, one uncoloured vacancy (letters L, R, M0; bits are spectators)
  c15 : C75 15-letter, two vacancy colours with C82's automaton (d = 17 lift)
  u12 : C80 12-letter, one uncoloured vacancy (letters L0 L1 R0 R1 E0 E1; bits kept)
  c12 : C80 12-letter, two vacancy colours (the d = 14 lift asked about)
Coloured hopping: V_q a -> a V_o(a) for every q in I(a), with h.c.; penalty J*K on
every extended pair (material-material, material-vacancy, vacancy-vacancy)
violating o(left) in I(right). Observable: ground-state weight on words whose
material sequence (vacancies deleted, cyclic) has a bad pair.
Also: (i) the 12- and 15-letter grammars factor through a two-state automaton,
(ii) colouring uniqueness and N_material_bad <= N_bad by exhaustion."""
import json, itertools, os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import my12

out = {}

# ==== automata
def automaton(letters, allowed):
    """Find o: letters->{0,1}, I: letters->subset with allowed(a,b) iff o(a) in I(b)."""
    rows = {}
    for a in letters:
        rows.setdefault(tuple(allowed(a, b) for b in letters), []).append(a)
    if len(rows) > 2:
        return None
    o = {}
    for k, (row, group) in enumerate(sorted(rows.items())):
        for a in group:
            o[a] = k
    I = {b: tuple(sorted(set(o[a] for a in letters if allowed(a, b)))) for b in letters}
    ok = all(allowed(a, b) == (o[a] in I[b]) for a in letters for b in letters)
    return dict(o=o, I=I, factorizes=ok)


L15 = ['L0', 'L1', 'R0', 'R1'] + [f'{k}{s}{a}' for k in 'DX' for s in 'AB' for a in (0, 1)] + ['M0', 'MA', 'MB']
def allowed15(l, r):
    if l[0] in 'LM':
        return True
    return r in ('R0', 'R1', 'M0')

aut12 = automaton(my12.LETTERS, my12.allowed)
aut15 = automaton(L15, allowed15)
out['automaton_12'] = dict(factorizes=aut12['factorizes'], o=aut12['o'], I={k: list(v) for k, v in aut12['I'].items()})
out['automaton_15'] = dict(factorizes=aut15['factorizes'], o=aut15['o'], I={k: list(v) for k, v in aut15['I'].items()})
# C82's table, for comparison (o=1 for R,D,X)
c82_o = {a: int(a[0] in 'RDX') for a in L15}
c82_I = {a: ((0, 1) if (a[0] == 'R' or a == 'M0') else (0,)) for a in L15}
def same_up_to_relabel(aut, o_ref, I_ref):
    for f in (lambda x: x, lambda x: 1 - x):
        if all(f(aut['o'][a]) == o_ref[a] and tuple(sorted(f(q) for q in aut['I'][a])) == I_ref[a] for a in o_ref):
            return True
    return False
out['automaton_15']['matches_C82_table_up_to_relabel'] = same_up_to_relabel(aut15, c82_o, c82_I)
print(out['automaton_12'], out['automaton_15'], flush=True)


def ext_allowed(aut):
    o, I = aut['o'], aut['I']
    def f(a, b):
        oa = int(a[1]) if a in ('V0', 'V1') else o[a]
        Ib = (int(b[1]),) if b in ('V0', 'V1') else I[b]
        return oa in Ib
    return f


# colouring uniqueness and N_mat_bad <= N_bad, exhaustive on short rings
def exhaust(letters, aut, allowed, Lring):
    f = ext_allowed(aut)
    alph = letters + ['V0', 'V1']
    worst = 0; legal_material = 0; uniq_fail = 0
    count = {}
    for w in itertools.product(alph, repeat=Lring):
        mat = [x for x in w if x[0] != 'V']
        if not mat:
            continue
        nb = sum(not f(w[i], w[(i + 1) % Lring]) for i in range(Lring))
        nm = sum(not allowed(mat[i], mat[(i + 1) % len(mat)]) for i in range(len(mat)))
        worst = max(worst, nm - nb)
        if nb == 0:
            key = tuple('V' if x[0] == 'V' else x for x in w)
            count[key] = count.get(key, 0) + 1
            if nm:
                uniq_fail += 1
    multi = sum(1 for v in count.values() if v != 1)
    return dict(ring=Lring, max_Nmat_minus_Next=worst, zero_penalty_coloured_words=sum(count.values()),
                uncoloured_patterns=len(count), patterns_with_nonunique_colouring=multi,
                zero_penalty_but_material_illegal=uniq_fail)

out['exhaust_12'] = [exhaust(my12.LETTERS, aut12, my12.allowed, L) for L in (3, 4, 5)]
out['exhaust_15'] = [exhaust(L15, aut15, allowed15, L) for L in (3, 4)]
print(out['exhaust_12'], out['exhaust_15'], flush=True)


# ==== reduced zero-excitation models
def model(kind):
    if kind in ('u15', 'c15'):
        moves = {}
        for a, b in [(('R', 'R'), ('L', 'R')), (('R', 'M'), ('L', 'M'))]:
            moves.setdefault(a, []).append(b); moves.setdefault(b, []).append(a)
        bad = lambda a, b: a == 'R' and b == 'L'
        o = {'L': 0, 'R': 1, 'M': 0}; I = {'L': (0,), 'R': (0, 1), 'M': (0, 1)}
        def cell(n, bits):
            return ['L'] + ['R'] * (n - 1) + ['M']
    else:
        moves = {}
        for a in (0, 1):
            for x, y in [((f'R0', f'R{a}'), (f'L{a}', 'R0')), ((f'R0', f'E{a}'), (f'L{a}', 'E0'))]:
                moves.setdefault(x, []).append(y); moves.setdefault(y, []).append(x)
        bad = lambda a, b: not my12.allowed(a, b)
        o = {l: aut12['o'][l] for l in ['L0', 'L1', 'R0', 'R1', 'E0', 'E1']}
        I = {l: aut12['I'][l] for l in o}
        def cell(n, bits):  # vacuum at j=1: L_b0 R0 R_b1 ... E_b(n-1)
            return [f'L{bits[0]}', 'R0'] + [f'R{b}' for b in bits[1:n - 1]] + [f'E{bits[n - 1]}']
    return moves, bad, o, I, cell


def run(kind, n, cells, nv, J, K, t=1., bits=None, cap=600000):
    moves, bad, o, I, cell = model(kind)
    coloured = kind[0] == 'c'
    bits = bits or [tuple([1] + [0] * (n - 1))] * cells
    mat = sum((cell(n, bits[c]) for c in range(cells)), [])
    L = len(mat) + nv
    vpos = set(int(round(i * L / nv)) % L for i in range(nv))
    seed = []; it = iter(mat)
    prev = mat[-1]
    for i in range(L):
        if i in vpos:
            seed.append(f'V{o[prev[0] if kind.endswith("15") else prev]}' if coloured else 'V')
        else:
            x = next(it); seed.append(x); prev = x
    key = (lambda l: l[0]) if kind.endswith('15') else (lambda l: l)
    def isV(x):
        return x[0] == 'V'
    def ext_ok(a, b):
        oa = int(a[1]) if isV(a) else o[key(a)]
        Ib = (int(b[1]),) if isV(b) else I[key(b)]
        return oa in Ib
    states = [tuple(seed)]; ix = {states[0]: 0}; r_, c_, v_ = [], [], []; cur = 0
    while cur < len(states):
        w = states[cur]
        def emit(zz, val):
            zz = tuple(zz)
            if zz not in ix:
                ix[zz] = len(states); states.append(zz)
                if len(states) > cap:
                    raise RuntimeError('cap')
            r_.append(ix[zz]); c_.append(cur); v_.append(val)
        diag = 0.
        for i in range(L):
            j = (i + 1) % L; a, b = w[i], w[j]
            if isV(a) != isV(b):
                if not coloured:
                    zz = list(w); zz[i], zz[j] = b, a; emit(zz, -t)
                else:
                    if isV(a):  # V_q b -> b V_o(b) for q in I(b)
                        if int(a[1]) in I[key(b)]:
                            zz = list(w); zz[i], zz[j] = b, f'V{o[key(b)]}'; emit(zz, -t)
                    else:       # a V_o(a) -> V_q a for q in I(a)
                        if int(b[1]) == o[key(a)]:
                            for q in I[key(a)]:
                                zz = list(w); zz[i], zz[j] = f'V{q}', a; emit(zz, -t)
            elif not isV(a):
                for a2, b2 in moves.get((a, b), []):
                    zz = list(w); zz[i], zz[j] = a2, b2; emit(zz, -J)
            if coloured:
                if not ext_ok(a, b):
                    diag += J * K
            elif not isV(a) and not isV(b) and bad(key(a) if kind.endswith('15') else a, key(b) if kind.endswith('15') else b):
                diag += J * K
        if diag:
            r_.append(cur); c_.append(cur); v_.append(diag)
        cur += 1
    N = len(states)
    H = sp.coo_matrix((v_, (r_, c_)), shape=(N, N)).tocsc()
    if N < 2500:
        ev, U = np.linalg.eigh(H.toarray()); E, psi = ev[0], U[:, 0]
    else:
        ev, U = eigsh(H, k=1, which='SA', tol=1e-13); E, psi = ev[0], U[:, 0]
    res = float(np.linalg.norm(H @ psi - E * psi))
    def mbad(w):
        m = [x for x in w if not isV(x)]
        if kind.endswith('15'):
            m = [x[0] for x in m]
        return sum(bad(m[i], m[(i + 1) % len(m)]) for i in range(len(m)))
    nb = np.array([mbad(w) for w in states])
    return dict(kind=kind, n=n, cells=cells, vacancies=nv, J=J, K=K, dim=N, E0=float(E), residual=res,
                illegal_material_weight=float(np.sum(abs(psi[nb > 0]) ** 2)))


rows = []
for kind in ('u15', 'c15', 'u12', 'c12'):
    for J in (1., 1e-2):
        for K in (10., 1e2, 1e3, 1e4):
            r = run(kind, 3, 2, 1, J, K); rows.append(r); print(r, flush=True)
out['leakage_n3_2cells_1vac'] = rows

# bit dependence of the coloured 12-letter ground energy (sector all-0 tail vs all-1 tail)
rows = []
for K in (1e2, 1e3, 1e4):
    e0 = run('c12', 3, 2, 1, 1., K, bits=[(1, 0, 0)] * 2)['E0']
    e1 = run('c12', 3, 2, 1, 1., K, bits=[(1, 1, 1)] * 2)['E0']
    f0 = run('c12', 3, 2, 1, 1., K, t=0., bits=[(1, 0, 0)] * 2)['E0']
    f1 = run('c12', 3, 2, 1, 1., K, t=0., bits=[(1, 1, 1)] * 2)['E0']
    rows.append(dict(K=K, K_times_dE_t1=K * (e0 - e1), K_times_dE_t0=K * (f0 - f1)))
    print(rows[-1], flush=True)
out['c12_bit_dependence'] = rows
json.dump(out, open('check_vacancy12.json', 'w'), indent=1)
