"""Small literal-ring ED and exact projection checks, all output local.

Run with python3 -B checks.py. No source program outside this folder runs
as a script. Read-only import of C75 is used for the independent bond audit.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import sys
sys.dont_write_bytecode = True
import itertools
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import bond16 as b

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'work_C75'))
import bond15 as old
import bond18 as parent


def maximum(x):
    return float(np.max(np.abs(x))) if np.size(x) else 0.


def ring_matrix(states, h):
    """Literal physical bond sum, including the oriented closing bond."""
    ix = {w: i for i, w in enumerate(states)}
    channels = {i: [(int(j), h[j, i]) for j in np.flatnonzero(h[:, i])]
                for i in range(256)}
    rows, cols, vals = [], [], []
    for col, w in enumerate(states):
        for x in range(len(w)):
            y = (x+1) % len(w)
            for row, val in channels[w[x]*16+w[y]]:
                z = list(w)
                z[x], z[y] = divmod(row, 16)
                z = tuple(z)
                if z not in ix:
                    raise AssertionError(('noninvariant basis', w, z))
                rows.append(ix[z]); cols.append(col); vals.append(val)
    return sp.coo_matrix((vals, (rows, cols)), shape=(len(states),)*2).tocsr()


def physical_states(L, words):
    pos = list(itertools.combinations(range(L), len(words[0])))
    states = []
    for x in pos:
        for w in words:
            z = [0]*L
            for i, a in zip(x, w):
                z[i] = a+1
            states.append(tuple(z))
    return states, pos


def squeezed_parts(L, words, positions, h):
    """Independent ordered-coordinate hopping and contact implementation."""
    N, nw = len(words[0]), len(words)
    wi = {w: i for i, w in enumerate(words)}
    pi = {x: i for i, x in enumerate(positions)}
    rr, cc, vv = [], [], []
    vr, vc, va = [], [], []
    channels = {i: [(int(j), h[j, i]) for j in np.flatnonzero(h[:, i])]
                for i in range(225)}
    for ip, x in enumerate(positions):
        occupied = set(x)
        for iw, w in enumerate(words):
            col = ip*nw+iw
            for k, site in enumerate(x):
                for direction in (-1, 1):
                    dest = (site+direction) % L
                    if dest in occupied:
                        continue
                    y = tuple(sorted((occupied-{site}) | {dest}))
                    word = w
                    if site == L-1 and dest == 0:
                        word = (w[-1],)+w[:-1]
                    if site == 0 and dest == L-1:
                        word = w[1:]+(w[0],)
                    rr.append(pi[y]*nw+wi[word]); cc.append(col); vv.append(-1.)
            for j in range(N):
                k = (j+1) % N
                if (x[j]+1) % L != x[k]:
                    continue
                for row, val in channels[w[j]*15+w[k]]:
                    z = list(w); z[j], z[k] = divmod(row, 15)
                    vr.append(ip*nw+wi[tuple(z)]); vc.append(col); va.append(val)
    shape = (len(positions)*nw,)*2
    return (sp.coo_matrix((vv, (rr, cc)), shape=shape).tocsr(),
            sp.coo_matrix((va, (vr, vc)), shape=shape).tocsr())


def full_two_particle_check(h):
    L, N = 4, 2
    words = list(itertools.product(range(15), repeat=N))
    states, positions = physical_states(L, words)
    T, V = squeezed_parts(L, words, positions, h)
    direct_t = ring_matrix(states, b.lift(h, J=0))
    direct_v = ring_matrix(states, b.lift(h, t=0, J=1))
    basis_errors = [maximum((T-direct_t).data), maximum((V-direct_v).data)]
    # All cyclic eigenspaces, both pi-flux sea branches included.
    swap = np.zeros((225, 225))
    for a, c in words:
        swap[c*15+a, a*15+c] = 1
    ev, internal = la.eigh(swap)
    Hs = h+swap@h@swap
    columns, pred_blocks = [], []
    alpha = N/L-np.sin(np.pi*N/L)**2/(N*L*np.sin(np.pi/L)**2)
    sea_energies = []
    for sign in (-1, 1):
        u = internal[:, np.abs(ev-sign)<1e-8]
        # Compress T in one representative internal vector.
        embed = np.kron(np.eye(len(positions)), u[:, :1])
        tc = embed.conj().T@(T@embed)
        e, z = la.eigh(tc)
        seas = z[:, np.abs(e-e[0])<1e-10]
        for k in range(seas.shape[1]):
            columns.append(np.kron(seas[:, k:k+1], u))
            pred_blocks.append(alpha*u.conj().T@Hs@u)
            sea_energies += [float(e[0])]*u.shape[1]
    P = np.column_stack(columns)
    expected = la.block_diag(*pred_blocks)
    projection_error = maximum(P.conj().T@(V@P)-expected)
    assert max(basis_errors) < 1e-13 and projection_error < 1e-11
    records = []
    # The lowest 120 levels all descend from the zero-twist sea.
    for J in (.0008, .0004, .0002):
        H = (T+J*V).toarray()
        e, u = la.eigh(H, subset_by_index=(0, 119), driver='evr')
        effective = np.diag(sea_energies)+J*expected
        pred = la.eigvalsh(effective)[:120]
        err = maximum(e-pred)
        residual = maximum(np.linalg.norm(H@u-u*e, axis=0))
        assert residual < 1e-10
        records.append(dict(J=J, retained=120, max_error=err,
                            error_over_J2=err/J**2, residual=residual))
    return dict(L=L, N=N, dimension=len(states), sea_dimension=P.shape[1],
                basis_errors=basis_errors, alpha=alpha,
                projection_error=projection_error, ed=records)


def alternating_words():
    # Two n=1 cells, complete local states, total Q=2. No grammar truncation.
    data = [i for i, x in enumerate(b.LET) if x[0] != 'M']
    markers = [i for i, x in enumerate(b.LET) if x[0] == 'M']
    out = []
    for phase in (0, 1):
        pools = [markers if j % 2 == phase else data for j in range(4)]
        for w in itertools.product(*pools):
            if sum(b.charge(b.LET[a]) for a in w) == 2:
                out.append(w)
    return out


def cell_ring_checks(h):
    words = alternating_words()
    records = []
    for L in (5, 6):
        states, pos = physical_states(L, words)
        T, V = squeezed_parts(L, words, pos, h)
        literal = ring_matrix(states, b.lift(h, J=.01))
        identity_error = maximum((literal-T-.01*V).data)
        assert identity_error < 1e-12
        for J in (.04, .02, .01):
            H = T+J*V
            e, u = eigsh(H, k=6, which='SA', tol=2e-11, ncv=60,
                         v0=np.cos(np.arange(H.shape[0])+.123), maxiter=10000)
            order = np.argsort(e); e, u = e[order], u[:, order]
            resid = maximum(np.linalg.norm(H@u-u*e, axis=0))
            assert resid < 1e-8
            records.append(dict(L=L, N=4, Q=2, markers=2,
                                material_words=len(words), dimension=len(states),
                                J=J, lowest_six=e.tolist(), residual=resid,
                                literal_squeezed_error=identity_error))
    return records


def main():
    p = b.params(1)
    lam = .13
    h = b.occupied_bond(p, lam)
    def oldlabel(x):
        k,s,a=x
        return ('e','',None) if (k,s)==('M','0') else ({'D':'d','X':'x','M':'m'}.get(k,k),s,a)
    permutation = [old.IX[oldlabel(x)] for x in b.LET]
    pairperm = [15*a+c for a in permutation for c in permutation]
    pp = parent.params(1, D=p['D'], Delta=p['Delta'], g=p['g'], td=p['td'],
                       K=p['K'], musel=p['musel'], field=p['field'])
    inherited = old.bond(pp, lam)[np.ix_(pairperm,pairperm)]
    H = b.lift(h, t=1.2, J=.013, mu=.037)
    charges = dict(N=np.r_[0,np.ones(15)],
                   Q=np.r_[0,[b.charge(x) for x in b.LET]],
                   M=np.r_[0,[x[0]=='M' for x in b.LET]])
    comm = {}
    for name, q in charges.items():
        total = np.add.outer(q,q).ravel()
        comm[name] = maximum(H*(total[:,None]-total[None,:]))
    out = dict(status='numerical_resolved_fixtures_not_the_n15_hierarchy',
               parameters=p, lambda_=lam, dimension=16, bond_shape=list(H.shape),
               nonzeros=int(np.count_nonzero(H)),
               independent_C75_error=maximum(h-inherited),
               hermiticity=maximum(H-H.conj().T), charge_commutators=comm,
               port_norm=float(la.norm(b.parts(p)[3],2)))
    assert out['independent_C75_error']<1e-13
    assert max(comm.values())==0 and out['hermiticity']==0
    out['two_particle'] = full_two_particle_check(h)
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)
    out['two_cell_rings'] = cell_ring_checks(h)
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out['two_cell_rings'],indent=2),flush=True)


if __name__ == '__main__':
    main()
