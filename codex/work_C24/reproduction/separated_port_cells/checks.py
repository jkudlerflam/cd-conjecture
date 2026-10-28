"""Independently realize left/right spectral weights with separated ports.

NumPy only. Imports the graded weight prescription, but builds a different
physical cell: two local boundary registers sharing one logical register.
No microscopic spatial Hamiltonian or CFT is claimed.
"""
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location(
    "graded", Path(__file__).parents[1] / "graded_route_a" / "checks.py")
graded = importlib.util.module_from_spec(spec)
spec.loader.exec_module(graded)


def complement(a):
    _, s, vh = np.linalg.svd(a.conj().T, full_matrices=True)
    return vh[np.count_nonzero(s > 1e-10):].conj().T


def one_port(weights, energies, parity):
    n = len(parity)
    k = len(energies)
    anc = (3*k + 4)//4
    boundary = 4*anc
    bell = np.zeros(boundary, complex)
    bell[0] = bell[3*anc] = 1/np.sqrt(2)
    g = np.kron(bell[:, None], np.eye(n))
    fs = [np.kron(f, np.eye(2*anc*n)) for f in (graded.Z, graded.X, graded.Y)]
    s = np.hstack([f@g for f in fs])
    bp = np.kron([1., -1., -1., 1.], np.ones(anc))
    pp = np.kron(bp, np.diag(parity).real)
    sp = np.kron([1., -1., -1.], np.diag(parity).real)
    dp = np.tile(sp, k)
    roots = []
    for w in weights:
        ev, u = np.linalg.eigh(w)
        roots.append((u*np.sqrt(ev))@u.conj().T)
    w = np.vstack(roots)
    j = np.zeros((boundary*n, 3*k*n), complex)
    for sign in (1., -1.):
        ix = np.flatnonzero(pp == sign)
        iy = np.flatnonzero(dp == sign)
        iz = np.flatnonzero(sp == sign)
        ig = np.flatnonzero(np.diag(parity).real == sign)
        ws = w[np.ix_(iy, iz)]
        ss = s[np.ix_(ix, iz)]
        gs = g[np.ix_(ix, ig)]
        from_basis = np.hstack([ws, complement(ws)])
        to_basis = np.hstack([ss, complement(np.hstack([gs, ss]))])
        j[np.ix_(ix, iy)] = to_basis[:, :len(iy)]@from_basis.conj().T
    h = (j*np.repeat(energies, 3*n))@j.conj().T
    h += energies[-1]*(np.eye(boundary*n)-g@g.conj().T-j@j.conj().T)
    h = (h+h.conj().T)/2
    realization = 0.
    for e in range(k):
        je = j[:, e*3*n:(e+1)*3*n]
        realization = max(realization, np.linalg.norm(s.conj().T@je@je.conj().T@s-weights[e]))
    return h, bell, bp, float(realization)


def run_case(name, ops_a, ops_b, grades, coupling):
    parity = graded.Z
    weights, energies, _, expected, weights_check = graded.spectral_weights(
        ops_a, ops_b, grades, coupling, parity)
    n = len(parity)
    hl, bell, bp, rl = one_port([w[:3*n, :3*n] for w in weights], energies, parity)
    hr, _, _, rr = one_port([w[3*n:, 3*n:] for w in weights], energies, parity)
    b = len(bell)
    d = b*b*n
    # Order is (left boundary register, right boundary register, logical).
    full_l = np.einsum('aibj,cd->acibdj', hl.reshape(b,n,b,n), np.eye(b)).reshape(d,d)
    full_r = np.einsum('cidj,ab->acibdj', hr.reshape(b,n,b,n), np.eye(b)).reshape(d,d)
    h = full_l+full_r
    g = np.kron(np.kron(bell[:,None],bell[:,None]),np.eye(n))
    pl = np.kron(np.kron(np.outer(bell,bell.conj()),np.eye(b)),np.eye(n))
    pr = np.kron(np.kron(np.eye(b),np.outer(bell,bell.conj())),np.eye(n))
    q = np.eye(d)-g@g.conj().T
    pp = np.kron(np.kron(bp,bp),np.diag(parity).real)
    left = [np.kron(f,np.eye(d//2)) for f in (graded.Z,graded.X,graded.Y)]
    right = [np.kron(np.eye(b),np.kron(f,np.eye(b*n//2)))
             for f in (graded.Z,graded.X,graded.Y)]
    vals,u = np.linalg.eigh(h)
    den = vals[:,None]+vals[None,:]
    inv = np.zeros_like(den)
    inv[den>.5] = 1/den[den>.5]
    amplitudes = sum(np.kron(u.conj().T@fr@g,u.conj().T@fl@g)
                     for fr,fl in zip(right,left))
    second = -amplitudes.conj().T@(inv.ravel()[:,None]*amplitudes)
    error = float(np.linalg.norm(second-expected))
    result = dict(
        name=name, logical_dimension=n, boundary_dimension=b,
        cell_dimension=d, two_cell_dimension=d*d,
        single_port_realization_error=max(rl,rr),
        code_error=float(np.linalg.norm(h@g)),
        code_degeneracy=int(sum(abs(vals)<1e-8)),
        gap=float(vals[n]),
        lower_bound_min_eigenvalue=float(np.linalg.eigvalsh(h-q).min()),
        parity_commutator=float(np.linalg.norm(pp[:,None]*h-h*pp[None,:])),
        other_code_projector_commutator=max(float(np.linalg.norm(full_l@pr-pr@full_l)),
                                          float(np.linalg.norm(full_r@pl-pl@full_r))),
        gadget_commutator=float(np.linalg.norm(full_l@full_r-full_r@full_l)),
        first_order_port_error=max(float(np.linalg.norm(g.conj().T@f@g)) for f in left+right),
        second_order_error=error,
        second_order_error_relative_to_exchange=error/weights_check['rho_squared'])
    assert result['code_degeneracy']==n
    assert result['gap']>.999999
    assert result['lower_bound_min_eigenvalue']>-1e-8
    assert error<1e-9
    return result


def main():
    start=time.time()
    cases=[run_case('odd_rank_one',[graded.X],[graded.Y],[1],np.ones((1,1))),
           run_case('mixed_parity',[graded.X,graded.Z],[graded.Y,graded.Z],[1,0],np.diag([1.,.3]))]
    out=dict(scope='Abstract cells only; no spatial compiler or central-charge estimate.',
             cases=cases,runtime_seconds=time.time()-start)
    Path(__file__).with_name('checks_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
