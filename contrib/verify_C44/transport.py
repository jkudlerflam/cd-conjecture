"""Direct construction of the marker move and its code overlap, C44 (9).

Two cells (b, b+1) -> (b+1, b). Each factor is built explicitly: the
selector head on its letter configurations, the typed vacancies with the
H_0 -> H_2 conversion, both helix tracks from ED of their parents, and
the biased C39 code. The result is compared with
beta_n = A_b^2 beta_v cos(phi/2)^(2b) and with C44's saved beta values.
"""
import json
import sys
import numpy as np
from core46 import Core, Params
from wrap import Charge, helix_track, plus


def hole_types(config, N):
    occ = set(config)
    first, last = min(config), max(config)
    t = []
    for x in range(N):
        if x in occ:
            t.append('O')
        elif x < first:
            t.append('H0')
        elif x > last:
            t.append('H2')
        else:
            t.append('H1')
    return t


def vacancy_overlap(k, U, wvac, convert=True):
    """<final charge code| move |initial charge code> for the (b,b+1)->(b+1,b) move.

    Initial: left cell N=k+1 (chi_1), right cell N=k+2 (chi_2), right first site H_0.
    Final:   left cell N=k+2 (chi_2), right cell N=k+1 (chi_1).
    The move relabels the transferred site H_0 -> H_2 (if convert) and keeps
    every occupied label. A final configuration is legal only if its typed
    labels match the canonical typing of its cell.
    """
    c1, c2 = Charge(k, 1, U, wvac), Charge(k, 2, U, wvac)
    b = k + 1
    tot = 0.0
    for cL in c1.configs:
        for cR in c2.configs:
            if 0 in cR:
                continue                      # first site of right cell must be a vacancy
            typesL = hole_types(cL, b)
            typesR = hole_types(cR, b + 1)
            assert typesR[0] == 'H0'
            moved = 'H2' if convert else 'H0'
            newL = typesL + [moved]
            cLn = cL                           # occupied positions unchanged in the left cell
            cRn = tuple(x - 1 for x in cR)
            if newL != hole_types(cLn, b + 1) or typesR[1:] != hole_types(cRn, b):
                continue                      # illegal typed word: zero overlap with the code
            amp = c1.chi[c1.index[cL]] * c2.chi[c2.index[cR]]
            amp *= c2.chi[c2.index[cLn]] * c1.chi[c1.index[cRn]]
            tot += amp
    # formula of C44 (9) / C8
    bv = sum(c1.chi[c1.index[tuple(x for x in range(b) if x != j)]] *
             c2.chi[c2.index[tuple(x for x in range(b + 1) if x not in (j, b))]]
             for j in range(b)) ** 2
    return tot, bv, c1, c2


def head_overlap(b):
    """Explicit selector move |R#><#L| on two cells, head ground states."""
    def ground(l):
        H = np.zeros((l - 1, l - 1))
        for i in range(l - 2):
            H[i, i + 1] = H[i + 1, i] = -1.0       # LX <-> XR
        e, v = np.linalg.eigh(H)
        g = v[:, 0] * np.sign(v[:, 0].sum())
        return g                               # index i-1 for head at position i
    gb, gb1 = ground(b), ground(b + 1)
    tot = 0.0
    # initial: left cell length b head at hL in 1..b-1; right cell length b+1 head hR in 1..b
    for hL in range(1, b):
        for hR in range(1, b + 1):
            if hR == 1:
                continue                       # right cell's site 1 must be L for '#L'
            amp_i = gb[hL - 1] * gb1[hR - 1]
            # final: left cell length b+1, head hL (its new last site is R); right cell length b, head hR-1
            amp_f = gb1[hL - 1] * gb[hR - 2]
            tot += amp_i * amp_f
    Ab = sum(np.sqrt(2 / b) * np.sin(np.pi * j / b) * np.sqrt(2 / (b + 1)) * np.sin(np.pi * j / (b + 1))
             for j in range(1, b))
    return tot, Ab


def helix_overlap(b, phi):
    """Both helix tracks over 2b+1 sites, ground states by ED, with the C8 (10) rotation."""
    def gs(N, pin):
        e, v = np.linalg.eigh(helix_track(N, phi, pin))
        assert e[0] < 1e-12 and e[1] > 1e-6
        return v[:, 0]
    Ub = np.diag([1.0, np.exp(1j * b * phi)])
    tot = np.exp(-1j * b * phi)
    for track, pin in (('A', 'first'), ('B', 'last')):
        def pinsite(N):
            return 0 if pin == 'first' else N - 1
        init = np.kron(gs(b, pinsite(b)), gs(b + 1, pinsite(b + 1)))       # sites 0..b-1 | b..2b
        fin = np.kron(gs(b + 1, pinsite(b + 1)), gs(b, pinsite(b)))        # sites 0..b | b+1..2b
        L = 2 * b + 1
        op = np.kron(np.kron(np.eye(2 ** b), Ub), np.eye(2 ** b))          # transferred site b
        ov = np.vdot(fin, op @ init)
        tot *= ov
    # remove the arbitrary global phases of ED ground states: compare magnitudes and
    # recompute the phase with the explicit product states
    prod = np.exp(-1j * b * phi)
    for track in ('A', 'B'):
        if track == 'A':
            ini = [plus(x * phi) for x in range(b)] + [plus(x * phi) for x in range(b + 1)]
            fi = [plus(x * phi) for x in range(b + 1)] + [plus(x * phi) for x in range(b)]
        else:
            ini = [plus((x - b + 1) * phi) for x in range(b)] + [plus((x - b) * phi) for x in range(b + 1)]
            fi = [plus((x - b) * phi) for x in range(b + 1)] + [plus((x - b + 1) * phi) for x in range(b)]
        for x in range(2 * b + 1):
            o = Ub if x == b else np.eye(2)
            prod *= np.vdot(fi[x], o @ ini[x])
    return abs(tot), prod, np.cos(phi / 2) ** (2 * b)


def payload_identity(n, k, wvac, U=1.0):
    """Code-to-code matrix of the move on the full (config x core word) cells."""
    core = Core(Params(n))
    G = core.code()
    tot, bv, c1, c2 = vacancy_overlap(k, U, wvac)
    b = k + 1
    d = 2 ** n
    # initial left: chi_1 x G ; right: chi_2 x G ; move maps (cL, cR) -> (cL, cR-1), words unchanged
    M = np.zeros((d * d, d * d), complex)
    for cL in c1.configs:
        for cR in c2.configs:
            if 0 in cR:
                continue
            cRn = tuple(x - 1 for x in cR)
            if hole_types(cL, b) + ['H2'] != hole_types(cL, b + 1) or hole_types(cR, b + 1)[1:] != hole_types(cRn, b):
                continue
            a = c1.chi[c1.index[cL]] * c2.chi[c2.index[cR]] * c2.chi[c2.index[cL]] * c1.chi[c1.index[cRn]]
            M += a * np.kron(G.conj().T @ G, G.conj().T @ G)
    return M, tot


if __name__ == '__main__':
    out = {}
    ref = json.load(open('../../codex/work_C44/source_results.json'))['cases']
    for n in (1, 2, 3):
        k = n + 2
        b = k + 1
        phi = 1.0 / b
        hd, Ab = head_overlap(b)
        habs, hprod, hform = helix_overlap(b, phi)
        for wvac in (0.16, 0.1):
            vd, vf, *_ = vacancy_overlap(k, 1.0, wvac)
            v0, *_ = vacancy_overlap(k, 1.0, wvac, convert=False)
            M, _ = payload_identity(n, k, wvac)
            beta_direct = hd * vd * habs
            beta_formula = Ab ** 2 * vf * hform
            c44 = [c['beta_length'] for c in ref if c['n'] == n and abs(c['wvac'] - wvac) < 1e-12]
            out['n%d_w%s' % (n, wvac)] = dict(
                head_direct=hd, head_formula=Ab ** 2,
                vacancy_direct=vd, vacancy_formula=vf, vacancy_unconverted=v0,
                helix_direct_abs=habs, helix_product=complex(hprod).real, helix_product_imag=complex(hprod).imag,
                helix_formula=hform, beta_direct=beta_direct, beta_formula=beta_formula,
                beta_C44=c44[0] if c44 else None,
                payload_identity_error=float(np.abs(M - vd * np.eye(M.shape[0])).max()))
            print(n, wvac, json.dumps(out['n%d_w%s' % (n, wvac)]))
    json.dump(out, open('transport.json', 'w'), indent=1)
