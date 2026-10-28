"""Mobile endpoint penalties for a frustration-free core.

Direct local-bond assembly in the invariant ordered-tag sector is compared
with the claimed complete code, ports, gap bound, and cut-shift overlap.
"""
import itertools
import json
import math
import numpy as np
from vacancy_parent import charge_hamiltonian


def core_code(q=3, k=3):
    v = np.zeros((q**k, k+1))
    for bs in itertools.product(range(2), repeat=k):
        weight = sum(bs)
        v[np.ravel_multi_index(bs, (q,)*k), weight] = 1/np.sqrt(math.comb(k, weight))
    return v


def swap_bond(q):
    swap = np.eye(q*q).reshape(q, q, q, q).transpose(1, 0, 2, 3).reshape(q*q, q*q)
    return (np.eye(q*q)-swap)/2


def build_sector(q, k, p):
    n = k+p
    holes = list(itertools.combinations(range(n), p))
    tags = [0]+[1]*(k-2)+[2]  # L, B, R, kept in occupied order.
    states = []
    for hs in holes:
        for payload in itertools.product(range(q), repeat=k):
            seq = iter([1+tag*q+a for tag, a in zip(tags, payload)])
            states.append(tuple(0 if j in hs else next(seq) for j in range(n)))
    return holes, states


def direct_local_h(states, q, U, w, g, a, b):
    index = {s:i for i, s in enumerate(states)}
    h0 = swap_bond(q)
    h = np.zeros((len(states), len(states)), complex)
    for col, s in enumerate(states):
        for j in range(len(s)-1):
            x, y = s[j:j+2]
            if (x == 0) != (y == 0):
                h[col, col] += U
                target = list(s); target[j], target[j+1] = y, x
                h[index[tuple(target)], col] -= w
            if x and y:
                tx, ax = divmod(x-1, q); ty, ay = divmod(y-1, q)
                for cx in range(q):
                    for cy in range(q):
                        coeff = g*h0[cx*q+cy, ax*q+ay]
                        if coeff:
                            target = list(s)
                            target[j], target[j+1] = 1+tx*q+cx, 1+ty*q+cy
                            h[index[tuple(target)], col] += coeff
            # Oriented on-site penalties in each bond, not added separately.
            if x and (x-1)//q == 0:
                ax = (x-1) % q
                for cx in range(q):
                    if a[cx, ax]:
                        target = list(s); target[j] = 1+cx
                        h[index[tuple(target)], col] += a[cx, ax]
            if y and (y-1)//q == 2:
                ay = (y-1) % q
                for cy in range(q):
                    if b[cy, ay]:
                        target = list(s); target[j+1] = 1+2*q+cy
                        h[index[tuple(target)], col] += b[cy, ay]
    return h


def port_matrix(states, q, f, end, tag):
    index = {s:i for i, s in enumerate(states)}
    out = np.zeros((len(states), len(states)), complex)
    for col, s in enumerate(states):
        letter = s[end]
        if not letter or (letter-1)//q != tag:
            continue
        old = (letter-1) % q
        for new in range(q):
            if f[new, old]:
                target = list(s); target[end] = 1+tag*q+new
                out[index[tuple(target)], col] += f[new, old]
    return out


def tensor_dictionary(states1, frame1, states2, frame2):
    out = {}
    for s, row in zip(states1, frame1):
        if not np.any(row):
            continue
        for t, other in zip(states2, frame2):
            if np.any(other):
                out[s+t] = np.kron(row, other)
    return out


def one_case(w):
    q, k, U, g = 3, 3, 1.0, 0.7
    a, b = np.diag([0., 0., .8]), np.diag([0., 0., 1.1])
    h0, v = swap_bond(q), core_code(q, k)
    core = g*(np.kron(h0, np.eye(q))+np.kron(np.eye(q), h0))
    core += np.kron(a, np.eye(q*q))+np.kron(np.eye(q*q), b)
    ec = np.linalg.eigvalsh(core)
    d = v.shape[1]
    gamma = ec[d]
    assert np.linalg.norm(core@v) < 1e-13
    fs = [np.array([[0,1,0],[1,0,0],[0,0,0.]]),
          np.array([[0,-1j,0],[1j,0,0],[0,0,0]]), np.diag([1.,-1.,0.])]
    frames, all_states, chis, charge_states, results = [], [], [], [], []
    for p in (1, 2):
        holes, hc, chi, charge_e = charge_hamiltonian(k, p, U, w)
        check_holes, states = build_sector(q, k, p)
        assert holes == check_holes
        h = direct_local_h(states, q, U, w, g, a, b)
        frame = np.kron(chi[:, None], v)
        ev = np.linalg.eigvalsh(h)
        residual = np.linalg.norm(h@frame-charge_e[0]*frame)
        code_error = np.linalg.norm(frame.T@frame-np.eye(d))
        occupation = 1-sum(chi[i]**2 for i, hs in enumerate(holes) if 0 in hs)
        port_error = 0.
        for f in fs:
            for end, tag, logical in [
                (0, 0, v.T@np.kron(f, np.eye(q*q))@v),
                (-1, 2, v.T@np.kron(np.eye(q*q), f)@v)]:
                phys = port_matrix(states, q, f, end, tag)
                port_error = max(port_error, np.linalg.norm(frame.T@phys@frame-occupation*logical))
        adj = []
        for j in range(k-1):
            prob = 0.
            for weight, hs in zip(chi**2, holes):
                occupied = [l for l in range(k+p) if l not in hs]
                prob += weight*(occupied[j+1] == occupied[j]+1)
            adj.append(prob)
        delta, pmin = charge_e[1]-charge_e[0], min(adj)
        M = g*(k-1)*np.linalg.norm(h0, 2)+np.linalg.norm(a, 2)+np.linalg.norm(b, 2)
        bound = delta*pmin*gamma/((np.sqrt(delta)+np.sqrt(M))**2+pmin*gamma)
        gap = ev[d]-ev[0]
        assert max(residual, code_error, port_error) < 1e-11
        assert abs(ev[0]-charge_e[0]) < 1e-11
        assert ev[d-1]-ev[0] < 1e-11 and gap >= bound-1e-12
        frames.append(frame); all_states.append(states); chis.append(chi); charge_states.append(holes)
        results.append(dict(holes=p, length=k+p, sector_dimension=len(states),
                            code_residual=float(residual), code_error=float(code_error),
                            port_error=float(port_error), ground_degeneracy=d,
                            full_gap=float(gap), lower_gap_bound=float(bound),
                            endpoint_occupation=float(occupation), min_adjacency=float(pmin)))
    before = tensor_dictionary(all_states[0], frames[0], all_states[1], frames[1])
    after = tensor_dictionary(all_states[1], frames[1], all_states[0], frames[0])
    overlap = np.zeros((d*d, d*d))
    for state in before.keys() & after.keys():
        overlap += np.outer(after[state], before[state])
    lookup = {s:i for i, s in enumerate(charge_states[1])}
    amplitude = sum(chis[0][i]*chis[1][lookup[(s[0],k+1)]]
                    for i, s in enumerate(charge_states[0]))
    beta = amplitude**2
    transport_error = np.linalg.norm(overlap-beta*np.eye(d*d))
    assert transport_error < 1e-11
    return dict(q=q, microscopic_dimension=3*q+1, k=k, w=w, U=U, g=g,
                core_gap=float(gamma), cells=results, beta=float(beta),
                transport_error=float(transport_error))


if __name__ == '__main__':
    out = [one_case(w) for w in (.3, .5)]
    with open('tagged_vacancy_checks.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
