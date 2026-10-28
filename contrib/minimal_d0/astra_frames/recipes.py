"""Explicit finite gate lists. Sites 0,...,m-1 are s_1,...,s_m; site m is c.
Gates in each list are chronological, so the last gate is leftmost in W.
"""
import itertools
import json
import time
import numpy as np
from verify import (I, PAULI, NAMES, target_terms, normalizer, inverse_gate,
                    character, eta, circuit_matrix, pauli)

def controlled_frame(m, qlabels):
    """K = controlled-(Q*) on s_1 -> spectators, acting trivially on c."""
    if not any(qlabels):
        return []
    c = m
    gates = [('SWAP', 0, c)]
    for j, a in enumerate(qlabels, 1):
        if a == 1:
            gates.append(('CX', c, j))
        elif a == 2:
            gates.extend([('S', j), ('CX', c, j), ('Sd', j)])
        elif a == 3:
            gates.extend([('H', j), ('CX', c, j), ('H', j)])
    gates.append(('SWAP', 0, c))
    return gates

def make_recipe(m, qlabels, uniform=True):
    qlabels = tuple(qlabels)
    assert len(qlabels) == m-1 and all(a in range(4) for a in qlabels)
    q = 2**m
    base_labels = (0,) * (m-1) if uniform else qlabels
    k = controlled_frame(m, qlabels) if uniform else []
    kinv = [inverse_gate(g) for g in reversed(k)]
    recipe = []
    for r, b, sigma in target_terms(m, base_labels):
        lr, lb = normalizer(r+(0,)), normalizer(b)
        gc = lb + [inverse_gate(g) for g in reversed(lr)]
        for u in itertools.product(range(4), repeat=m):
            if character(r, u) == sigma * eta(r):
                gu = [(NAMES[a], j) for j, a in enumerate(u) if a]
                recipe.append({'weight': 4/q**2, 'gates': kinv+gc+gu+k})
    assert len(recipe) == q*q*(q*q-2)//2
    return recipe

def verify_recipe(m, qlabels, recipe):
    q, d = 2**m, 2*(2**m)**2
    omega = np.eye(q).reshape(q*q) / np.sqrt(q)
    pi0 = np.outer(omega, omega)
    target, local = np.zeros((d,d),complex), np.zeros((d,d),complex)
    for a in (1,2):
        f = pauli((a,)+tuple(qlabels))
        fp = np.kron(f,np.eye(q))
        target += np.kron(fp@pi0+pi0@fp,PAULI[a])
        local += 2/q**2 * np.kron(np.eye(q),np.kron(f.T,PAULI[a]))
    total, cache = np.zeros((d,d),complex), {}
    for entry in recipe:
        W = circuit_matrix(entry['gates'],m+1,cache)
        V = W.conj().T.reshape(2*q,q,2).transpose(1,0,2).reshape(d,2)/np.sqrt(q)
        total += entry['weight'] * V@V.conj().T
    residual = total-(2-4/q**2)*np.eye(d)+local-target
    K = circuit_matrix(controlled_frame(m,qlabels),m+1,cache)
    qs = pauli(qlabels).conj()
    want_k = np.kron(np.diag([1,0]),np.eye(q//2))+np.kron(np.diag([0,1]),qs)
    return {'m':m,'Q':''.join(NAMES[a] for a in qlabels) or 'empty',
            'frames':len(recipe),'max_gates':max(len(x['gates']) for x in recipe),
            'residual_max_abs':float(np.max(np.abs(residual))),
            'residual_fro':float(np.linalg.norm(residual)),
            'relative_fro':float(np.linalg.norm(residual)/np.linalg.norm(target)),
            'controlled_frame_max_abs':float(np.max(np.abs(K-np.kron(want_k,I))))}

if __name__ == '__main__':
    start=time.monotonic()
    rows, export = [], {}
    for m in (1,2,3):
        for labels in itertools.product(range(4),repeat=m-1):
            recipe=make_recipe(m,labels,uniform=True)
            row=verify_recipe(m,labels,recipe)
            rows.append(row)
            print(json.dumps(row),flush=True)
            if m <= 2:
                export[f'm={m},Q={row["Q"]}']=recipe
    with open('explicit_recipes_m1_m2.json','w') as f:
        json.dump(export,f,indent=2)
    with open('uniformity_verification.json','w') as f:
        json.dump({'rows':rows,'elapsed_seconds':time.monotonic()-start},f,indent=2)
    print('Elapsed seconds:',time.monotonic()-start,flush=True)
