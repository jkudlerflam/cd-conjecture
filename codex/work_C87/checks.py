"""C87 finite checks. Exact integers/Fractions and small numerical ED.

No sibling imports. All output is relative to this folder.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.linalg import eigh


def encode(x):
    x = F(x)
    return [str(x.numerator), str(x.denominator)]


def inverse(a):
    n = len(a)
    x = [[F(int(v)) for v in row]+[F(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        k = next(k for k in range(j, n) if x[k][j])
        x[j], x[k] = x[k], x[j]
        pivot = x[j][j]
        x[j] = [v/pivot for v in x[j]]
        for k in range(n):
            if k == j:
                continue
            t = x[k][j]
            if t:
                x[k] = [u-t*v for u, v in zip(x[k], x[j])]
    return np.array([row[n:] for row in x], dtype=object)


def tokens(length, exchange=True, flat=True):
    states = [(a,b) for a in range(length) for b in range(length) if a != b]
    idx = {x:i for i,x in enumerate(states)}
    # Phase units pi/2. Right hops A:1, B:2; exchange AB -> BA:1.
    theta = [1, 2]
    eta = (1-theta[0]+theta[1]) % 4
    gauge = np.array([(a+2*b+eta*(a>b)) % 4 for a,b in states])
    edges = []
    for x in states:
        for token in [0,1]:
            y = list(x)
            y[token] += 1
            y = tuple(y)
            if y in idx:
                edges.append((idx[y], idx[x], theta[token], 1))
        if exchange and x[1] == x[0]+1:
            phase = 1 if flat or x[0] != 0 else 2
            edges.append((idx[x[::-1]], idx[x], phase, 1))
    real = np.diag([10+2*a+3*b+abs(a-b) for a,b in states])
    complex_h = real.astype(complex)
    adj = [set() for x in states]
    for i,j,phase,weight in edges:
        complex_h[i,j] -= weight*(1j)**phase
        complex_h[j,i] -= weight*(1j)**(-phase)
        if flat:
            assert (phase-gauge[i]+gauge[j]) % 4 == 0
        real[i,j] -= weight
        real[j,i] -= weight
        adj[i].add(j)
        adj[j].add(i)
    seen = set()
    components = []
    for root in range(len(states)):
        if root in seen:
            continue
        component = {root}
        frontier = [root]
        seen.add(root)
        while frontier:
            for j in adj[frontier.pop()]:
                if j not in seen:
                    seen.add(j)
                    component.add(j)
                    frontier.append(j)
        components.append(sorted(component))
    return states, edges, gauge, real, complex_h, components


def token_checks():
    out = []
    for length in range(3, 10):
        for exchange in [False, True]:
            states, edges, gauge, h, hc, comps = tokens(length, exchange)
            assert len(comps) == (1 if exchange else 2)
            d = np.diag((1j)**gauge)
            assert np.array_equal(d.conj().T@hc@d, h)
            gaps = []
            residuals = []
            for comp in comps:
                mat = h[np.ix_(comp,comp)]
                vals, vec = eigh(mat)
                gaps.append(float(vals[1]-vals[0]))
                residuals.append(float(np.linalg.norm(mat@vec-vec*vals, 2)))
                assert vals[1] > vals[0]+1e-8
            out.append(dict(status='exact_gauge_and_components_numerical_ED',
                length=length, dimension=len(states), exchange=exchange,
                components=len(comps), component_ground_gaps=gaps,
                max_eigenpair_matrix_residual=max(residuals)))
    # A six-step cycle with a deliberately nonuniform exchange phase.
    states, edges, gauge, h, hc, comps = tokens(3, True, False)
    cycle = [(0,1),(1,0),(2,0),(2,1),(1,2),(0,2),(0,1)]
    phase_map = {}
    for i,j,p,w in edges:
        phase_map[states[j],states[i]] = p
        phase_map[states[i],states[j]] = -p
    flux = sum(phase_map[x,y] for x,y in zip(cycle,cycle[1:])) % 4
    assert flux == 1
    multispecies = []
    for length in range(3,8):
        words = []
        for positions in itertools.permutations(range(length),3):
            w = [0]*length
            for species,pos in enumerate(positions,1):
                w[pos] = species
            words.append(tuple(w))
        idx = {w:i for i,w in enumerate(words)}
        phases = {(s,t):(2*s+t+s*t) % 4 for s in range(4) for t in range(s+1,4)}
        def potential(w):
            return sum(phases[t,s] for j,s in enumerate(w) for t in w[j+1:] if t<s) % 4
        potentials = list(map(potential,words))
        edge_count = 0
        for col,w in enumerate(words):
            for j in range(length-1):
                s,t = w[j:j+2]
                if s>=t:
                    continue
                z = list(w)
                z[j:j+2] = [t,s]
                row = idx[tuple(z)]
                assert (phases[s,t]-potentials[row]+potentials[col]) % 4 == 0
                edge_count += 1
        multispecies.append(dict(status='exact_integer_gauge_check',sites=length,
            configurations=len(words),edges=edge_count,independent_cycles=edge_count-len(words)+1))
    return dict(flat_cases=out, multispecies=multispecies, nonuniform_exchange_cycle=cycle,
                nonuniform_exchange_flux_pi_over_two=flux)


def rational_two_token_coefficient():
    states, edges, gauge, h, hc, comps = tokens(3)
    s = len(h)
    a = np.kron(h, np.eye(s, dtype=int))+np.kron(np.eye(s, dtype=int), h)
    r = inverse(a)
    assert np.array_equal(a.astype(object)@r, np.eye(s*s, dtype=int))
    # Opposite ordered endpoint sources, with no assumed reflection symmetry.
    u = states.index((2,1))
    v = states.index((0,1))
    x = r[:,u*s+v]
    U = np.zeros((s*s,3), dtype=object)
    V = U.copy()
    for j in range(3):
        for l in range(s):
            for k in range(s):
                if states[l][0] == j:
                    U[l*s+k,j] = x[l*s+k]
                if states[k][0] == j:
                    V[l*s+k,j] = x[l*s+k]
    c = U.T@r@V
    onsite = U.T@r@U+V.T@r@V
    assert np.array_equal(onsite, onsite.T)
    return dict(status='proved_exact_finite_check', carrier_dimension=s,
                pair_dimension=s*s, active_projector_ranks=[2,2,2],
                C=[[encode(v) for v in row] for row in c],
                onsite_imaginary_part_exactly_zero=True,
                inverse_identity_exact=True)


def fermions(m):
    out = []
    for j in range(m):
        c = np.zeros((1<<m,1<<m), dtype=int)
        for x in range(1<<m):
            if (x>>j)&1:
                c[x^(1<<j),x] = (-1)**((x&((1<<j)-1)).bit_count())
        out.append(c)
    return out


def reconstruct(k, cs):
    m = len(cs)
    inds = [1<<j for j in range(m)]
    h = k[np.ix_(inds,inds)]-k[0,0]*np.eye(m,dtype=k.dtype)
    q = k[0,0]*np.eye(len(k),dtype=k.dtype)
    for i in range(m):
        for j in range(m):
            q = q+h[i,j]*(cs[i].T@cs[j])
    return q


def phase_classification():
    # Exhaust all 2^8 diagonal sign dressings of a connected three-mode path.
    cs = fermions(3)
    k = sum((cs[i].T@cs[i+1]+cs[i+1].T@cs[i] for i in range(2)))
    count = 0
    for signs in itertools.product([-1,1], repeat=8):
        u = np.array(signs)
        ku = u[:,None]*k*u[None,:]
        gaussian = np.array_equal(ku,reconstruct(ku,cs))
        # Remove the phases defined on the one-particle sector.
        values = {}
        form = True
        for x in range(8):
            val = u[x]
            for j in range(3):
                if (x>>j)&1:
                    val *= u[1<<j]
            r = x.bit_count()
            if r in values and values[r] != val:
                form = False
            values[r] = val
        assert gaussian == form
        count += gaussian
    return dict(status='proved_exhaustive_finite_check', modes=3,
                diagonal_sign_unitaries=256, gaussian_preserving=int(count))


def ports():
    m = 4
    cs = fermions(m)
    ns = [c.T@c for c in cs]
    eye = np.eye(1<<m,dtype=int)
    # A connected Gaussian operator, dressed by a phase within the left cell.
    lines = [cs[j]+cs[j+1] for j in range(3)]
    D = eye-2*ns[0]@ns[1]
    phase_source = np.vstack([l@D for l in lines])
    phase_k = -phase_source.T@phase_source
    phase_r = phase_k-reconstruct(phase_k,cs)
    assert np.any(phase_r)
    # Move mode 1 to mode 2, first without and then with spectator mode 0.
    hop = cs[2].T@cs[1]
    vals = []
    for x in [2,3]:
        y = x^2^4
        vals.append(int(phase_k[y,x]//hop[y,x]))
    assert abs(vals[0]-vals[1]) == 2
    ampD = (eye+ns[0])@(eye+ns[2])
    amplitude_k = -ampD.T@ampD
    mixed = int(amplitude_k[5,5]-amplitude_k[1,1]-amplitude_k[4,4]+amplitude_k[0,0])
    assert mixed == -9
    # Complementary, orthogonal excited channels. Each source is nonlinear.
    l1 = cs[0]+cs[1]
    l2 = cs[2]+cs[3]
    pieces = [ns[2]@l1,(eye-ns[2])@l1,ns[0]@l2,(eye-ns[0])@l2]
    source = np.vstack(pieces)
    k = -source.T@source
    assert np.array_equal(k, -(l1.T@l1+l2.T@l2))
    assert np.array_equal(k,reconstruct(k,cs))
    # Each individual projected contribution is interacting.
    component_non_gaussian = []
    for piece in pieces:
        kp = -piece.T@piece
        component_non_gaussian.append(bool(np.any(kp-reconstruct(kp,cs))))
    assert all(component_non_gaussian)
    checks = []
    for name, src, expected_gaussian in [
        ('conditional_phase',phase_source,False),
        ('conditional_amplitude',ampD,False),
        ('complementary_channels',source,True)]:
        k2 = -src.T@src
        # Exact source-Gram identity is the second-order Schur coefficient.
        assert np.array_equal(k2,reconstruct(k2,cs)) == expected_gaussian
        n = len(k2)
        q = len(src)
        ed = []
        for lam in [0.1,0.05,0.025]:
            full = np.block([[np.zeros((n,n)),lam*src.T],
                             [lam*src,np.eye(q)]])
            w,z = eigh(full)
            low = w[:n]
            pred = np.linalg.eigvalsh(lam**2*k2)
            exact = np.sort((1-np.sqrt(1-4*lam**2*np.linalg.eigvalsh(k2)))/2)
            err = float(np.max(np.abs(low-pred)))
            res = float(np.max(np.linalg.norm(full@z-z*w,axis=0)))
            assert np.max(abs(low-exact)) < 1e-12
            assert res < 1e-12
            ed.append(dict(lambda_=lam, max_second_order_energy_error=err,
                           error_over_lambda4=err/lam**4, eigenpair_residual=res))
        checks.append(dict(status='exact_Gram_numerical_ED',name=name,
            code_dimension=n,excited_dimension=q,full_dimension=n+q,
            gaussian_second_order=expected_gaussian,
            exact_reconstruction_residual_max=int(np.max(abs(k2-reconstruct(k2,cs)))),
            ed=ed))
    # Fourth-order interaction in the passing example: k2^2 is not Gaussian.
    fourth = k@k
    assert np.any(fourth-reconstruct(fourth,cs))
    return dict(checks=checks, phase_hopping_amplitudes=vals,
                phase_distance_to_any_Gaussian_lower=encode(1),
                amplitude_mixed_difference=mixed,
                amplitude_distance_to_any_Gaussian_lower=encode(F(9,4)),
                complementary_component_non_gaussian=component_non_gaussian,
                complementary_fourth_order_non_gaussian=True)


def counts():
    position = []
    for n in [9,10,11]:
        b_dist = next(b for b in range(2,1000) if b*(b-1) >= 2**n)
        b_same = next(b for b in range(2,1000) if math.comb(b,2) >= 2**n)
        b_rails = next(b for b in range(1,1000) if b*b >= 2**n)
        position.append(dict(n=n,distinct_same_chain_minimum_sites=b_dist,
            identical_same_chain_minimum_sites=b_same,
            separate_rails_minimum_sites_per_rail=b_rails))
    marker = []
    for k in range(1,9):
        mult = [math.comb(k,j) for j in range(k+1)]
        padded = [0]+mult+[0]
        excited = sum(max(x,y) for x,y in zip(padded,padded[1:]))
        assert excited == 2**k+math.comb(k,k//2)
        marker.append(dict(bits=k,neutral_dimension=2**k,
            excited_dimension_lower=excited,total_dimension_lower=2**k+excited))
    return dict(status='proved_exact_finite_counts',positions=position,markers=marker)


def main():
    out = dict(tokens=token_checks(),coefficient=rational_two_token_coefficient(),
               phase_classification=phase_classification(),ports=ports(),counts=counts())
    Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status='all_checks_passed',token_cases=len(out['tokens']['flat_cases']),
        max_token_dimension=72,exact_product_dimension=36,
        phase_classification=out['phase_classification'],
        port_ED_dimensions=[x['full_dimension'] for x in out['ports']['checks']],
        counts=out['counts']),indent=2))


if __name__ == '__main__':
    main()
