"""(9) Is the C84 boundary pin dynamically needed?
 (a) Frozen chain: for the C84 bond WITHOUT the pin penalty, every literal move and
     port transition on every ring word of length 3, 4, 5 over the 11 letters
     conserves the number of adjacent (L1, marker) pairs. Then Pi_pin commutes with
     H_pinless and H_C84 = H_pinless + K Pi_pin acts identically on Pi_pin = 0.
 (b) Mobile chain: with an uncoloured vacancy V, every material move (on adjacent
     material pairs) and every vacancy hop on rings of length 4, 5 over 12 letters
     conserves the number of squeezed (L1, marker) pairs.
 (c) Two-colour lift of the PINLESS 11-letter grammar (C83's automaton with E_1
     deleted): automaton reproduces the pinless table, unique colouring and
     N_material_bad <= N_bad on rings of 3..5 sites, and K^-2 suppression of
     illegal-material weight for two p = 3 cells with nonzero bits and one vacancy."""
import sys, os, json, itertools
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import my11

M = my11.LETTERS
MK = set(my11.MARKERS)
pm = my11.default_params(3, g=.13, td=.21, phiA=.37, psiA=1.3)
mv = my11.rule_table(my11.move_rules(pm) + my11.port_rules(pm))
out = {}


def pins_ring(w):
    L = len(w)
    return sum(w[i] == 'L1' and w[(i + 1) % L] in MK for i in range(L))


# (a)
fr = {}
for L in (3, 4, 5):
    viol = 0; trans = 0
    for w in itertools.product(M, repeat=L):
        c0 = pins_ring(w)
        for i in range(L):
            j = (i + 1) % L
            for a2, b2, val in mv.get((w[i], w[j]), []):
                z = list(w); z[i] = a2; z[j] = b2
                trans += 1
                if pins_ring(z) != c0:
                    viol += 1
    fr[L] = dict(transitions=trans, pin_count_changes=viol)
    print('frozen', L, fr[L], flush=True)
out['frozen_pin_conservation'] = fr


# (b)
def squeezed_pins(w):
    m = [l for l in w if l != 'V']
    k = len(m)
    return sum(m[i] == 'L1' and m[(i + 1) % k] in MK for i in range(k)) if k else 0


mob = {}
for L in (4, 5):
    viol = 0; trans = 0
    for w in itertools.product(M + ['V'], repeat=L):
        if all(l == 'V' for l in w):
            continue
        c0 = squeezed_pins(w)
        for i in range(L):
            j = (i + 1) % L
            a, b = w[i], w[j]
            outs = []
            if a != 'V' and b != 'V':
                outs = [(a2, b2) for a2, b2, _ in mv.get((a, b), [])]
            elif a == 'V' and b != 'V':
                outs = [(b, 'V')]
            elif a != 'V' and b == 'V':
                outs = [('V', a)]
            for a2, b2 in outs:
                z = list(w); z[i] = a2; z[j] = b2
                trans += 1
                if squeezed_pins(z) != c0:
                    viol += 1
    mob[L] = dict(transitions=trans, squeezed_pin_count_changes=viol)
    print('mobile', L, mob[L], flush=True)
out['mobile_pin_conservation'] = mob

# (c) two-colour lift of the pinless grammar
O2 = {'L0': 0, 'L1': 0, 'E': 0, 'MA': 0, 'MB': 0, 'R0': 1, 'R1': 1, 'DA': 1, 'DB': 1, 'XA': 1, 'XB': 1,
      'V0': 0, 'V1': 1}
I2 = {'L0': (0,), 'L1': (0,), 'DA': (0,), 'DB': (0,), 'XA': (0,), 'XB': (0,), 'MA': (0,), 'MB': (0,),
      'R0': (0, 1), 'R1': (1,), 'E': (0, 1), 'V0': (0,), 'V1': (1,)}
ok2 = lambda a, b: O2[a] in I2[b]
pinless = lambda a, b: my11.allowed(a, b) or (a == 'L1' and b in MK)
out['automaton_reproduces_pinless_table'] = all(ok2(a, b) == pinless(a, b) for a in M for b in M)
LET2 = M + ['V0', 'V1']
isv = lambda l: l[0] == 'V'
ex = {}
for L in (3, 4, 5):
    viol = 0; seen = {}
    for w in itertools.product(LET2, repeat=L):
        occ = [l for l in w if not isv(l)]
        if not occ:
            continue
        nmb = sum(not pinless(occ[i], occ[(i + 1) % len(occ)]) for i in range(len(occ)))
        nb = sum(not ok2(w[i], w[(i + 1) % L]) for i in range(L))
        viol += nmb > nb
        if nb == 0:
            key = tuple('V' if isv(l) else l for l in w)
            seen[key] = seen.get(key, 0) + 1
    legal = 0; bad = 0
    for pat in itertools.product((0, 1), repeat=L):
        k = sum(pat)
        if not k:
            continue
        for mat in itertools.product(M, repeat=k):
            if all(pinless(mat[i], mat[(i + 1) % k]) for i in range(k)):
                legal += 1
                it = iter(mat)
                key = tuple(next(it) if s else 'V' for s in pat)
                bad += seen.get(key, 0) != 1
    ex[L] = dict(legal_squeezed=legal, zero_penalty=sum(seen.values()), N_material_bad_gt_N_bad=int(viol),
                 colouring_not_unique_or_missing=bad)
    print('2col', L, ex[L], flush=True)
out['two_colour_pinless_rings'] = ex


def bond13(pmk, lam, J, t):
    ls = LET2; dd = len(ls); ix = {l: i for i, l in enumerate(ls)}
    pq = dict(pmk); K = pq['K']; pq['K'] = 0.
    hm = my11.bond(pq, lam)
    hm += np.diag([K * (not pinless(a, b)) for a in M for b in M])
    h = np.zeros((dd * dd, dd * dd), complex)
    occ = [dd * i + j for i in range(11) for j in range(11)]
    h[np.ix_(occ, occ)] = J * hm
    for a in M:
        for q in I2[a]:
            src = ix[f'V{q}'] * dd + ix[a]; dst = ix[a] * dd + ix[f'V{O2[a]}']
            h[dst, src] += -t; h[src, dst] += -t
    for a in ls:
        for b in ls:
            if (isv(a) or isv(b)) and not ok2(a, b):
                h[ix[a] * dd + ix[b], ix[a] * dd + ix[b]] += J * K
    return h, ls


def ring_closure(seed, h, ls, cap=400000):
    dd = len(ls); ixx = {l: i for i, l in enumerate(ls)}
    T = {}
    for col in range(dd * dd):
        for row in np.flatnonzero(abs(h[:, col]) > 0):
            if row != col:
                T.setdefault(col, []).append((row, h[row, col]))
    diag = np.real(np.diag(h))
    L = len(seed); states = [tuple(seed)]; idx = {states[0]: 0}; r, c, v = [], [], []; cur = 0
    while cur < len(states):
        w = states[cur]; dg = 0.
        for i in range(L):
            col = ixx[w[i]] * dd + ixx[w[(i + 1) % L]]
            dg += diag[col]
            for row, val in T.get(col, []):
                z = list(w); z[i] = ls[row // dd]; z[(i + 1) % L] = ls[row % dd]; z = tuple(z)
                if z not in idx:
                    idx[z] = len(states); states.append(z)
                    if len(states) > cap:
                        raise RuntimeError('cap')
                r.append(idx[z]); c.append(cur); v.append(val)
        r.append(cur); c.append(cur); v.append(dg); cur += 1
    N = len(states)
    return states, sp.coo_matrix((v, (r, c)), shape=(N, N)).tocsc()


c1 = my11.encode_word(0, 'v', (1, 1)); c2 = my11.encode_word(1, 'v', (0, 1))
seed = list(c1[:2] + ('V',) + c1[2:] + c2)
for i, l in enumerate(seed):
    if l == 'V':
        seed[i] = f'V{O2[seed[i - 1]]}'
seed = tuple(seed)
pmc = my11.default_params(3, D=4., g=.07, td=.003, phiA=np.pi / 6, phiB=0., psiA=np.pi, psiB=np.pi / 2, K=0.)
leak = []
for J, Ks in [(1., (1e2, 1e3, 1e4)), (.1, (1e3, 1e4, 1e5))]:
    for K in Ks:
        pk = dict(pmc); pk['K'] = K
        h, ls = bond13(pk, 0., J, 1.)
        st, H = ring_closure(seed, h, ls)
        sigma = -float(abs(H).sum(axis=1).max()) - 1
        e, U = eigsh(H, k=1, sigma=sigma, which='LM', tol=1e-13)
        v = U[:, 0]
        ill = np.array([any(not pinless(m[i], m[(i + 1) % len(m)]) for i in range(len(m)))
                        for m in ([l for l in w if not isv(l)] for w in st)])
        leak.append(dict(J=J, K=K, JK=J * K, dim=len(st), illegal_material_weight=float(np.sum(abs(v[ill]) ** 2)),
                         residual=float(np.linalg.norm(H @ v - e[0] * v))))
        print(leak[-1], flush=True)
out['two_colour_pinless_leakage'] = dict(seed=seed, runs=leak)
json.dump(out, open(os.path.join(HERE, 'check_pin.json'), 'w'), indent=1)
