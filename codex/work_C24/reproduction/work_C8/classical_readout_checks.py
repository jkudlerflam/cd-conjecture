"""C8 nonnegative joint ports and clean-readout finite histories.

Run with OPENBLAS_NUM_THREADS=1 python3 classical_readout_checks.py.
Only finite-dimensional algebra is tested, not a microscopic CFT.
"""
import itertools
import json
import numpy as np
from scipy import sparse as sp
from scipy.linalg import eigvalsh


def opnorm(a):
    a = a.toarray() if sp.issparse(a) else a
    return float(np.linalg.norm(a, 2))


def pairs_and_ops(q, n, graded=False):
    pairs = list(itertools.combinations(range(q), 2))
    if graded:
        pairs = [(a, b) for a, b in pairs if a < q // 2 <= b]
    assert len(pairs) >= 2 * n
    pairs = np.array(pairs[:2*n]).reshape(n, 2, 2)
    left, right = [], []
    for f in range(n):
        x, y = np.zeros((q, q), complex), np.zeros((q, q), complex)
        for neg, (a, b) in enumerate(pairs[f]):
            s = 1 - 2*neg
            x[a, b] = x[b, a] = s
            y[b, a], y[a, b] = 1j*s, -1j*s
        left.extend([sp.csr_matrix(x), sp.csr_matrix(y)])
        right.extend([sp.csr_matrix(y), sp.csr_matrix(-x)])
    return pairs, left, right


def bits(x, n):
    return [(x >> (n-1-j)) & 1 for j in range(n)]


def erased_bit(x, n, f):
    bs = bits(x, n)
    z = 0
    for j, b in enumerate(bs):
        if j != f:
            z = 2*z + b
    return z


def port_value(x, n, f, side, pairs):
    bs = bits(x, n)
    neg = sum(bs[:f] if side == 0 else bs[f+1:]) % 2
    return pairs[f, neg, bs[f]]


def cliffords(n):
    d = 2**n
    rows, cols, xvals, yvals = [], [], [], []
    out = []
    for f in range(n):
        rows, cols, xvals, yvals = [], [], [], []
        for x in range(d):
            bs = bits(x, n)
            y = x ^ (1 << (n-1-f))
            s = (-1)**sum(bs[:f])
            rows.append(y); cols.append(x)
            xvals.append(s); yvals.append(1j*s*(-1)**bs[f])
        out.extend([sp.csr_matrix((vals, (rows, cols)), shape=(d, d))
                    for vals in (xvals, yvals)])
    parity = sp.diags([(-1)**sum(bits(x, n)) for x in range(d)], dtype=float)
    return out, [1j*a@parity for a in out], parity


def static_code(q, n, pairs):
    d, e = 2**n, n*2**n
    rows, cols = [], []
    for x in range(d):
        for side in range(2):
            for f in range(n):
                env = (side*n+f)*(d//2)+erased_bit(x, n, f)
                v = port_value(x, n, f, side, pairs)
                l, r = (v, 0) if side == 0 else (0, v)
                rows.append((l*q+r)*e+env); cols.append(x)
    v = sp.csr_matrix((np.full(len(rows), 1/np.sqrt(2*n)), (rows, cols)),
                      shape=(q*q*e, d))
    return v, e


def static_check(q, n, graded=False):
    pairs, left, right = pairs_and_ops(q, n, graded)
    v, e = static_code(q, n, pairs)
    gs, rs, parity = cliffords(n)
    d = 2**n
    errs = [opnorm(v.T@v-sp.eye(d))]
    for ops, targets, side in [(left, gs, 0), (right, rs, 1)]:
        for f, target in zip(ops, targets):
            phys = sp.kron(sp.kron(f, sp.eye(q)), sp.eye(e)) if side == 0 else \
                   sp.kron(sp.kron(sp.eye(q), f), sp.eye(e))
            errs.append(opnorm(v.T@phys@v-target/(2*n)))
    # The flat parent is stoquastic because its off-diagonal entries are -VV^T.
    p = v@v.T
    off = p - sp.diags(p.diagonal())
    min_positive_projector_entry = float(off.data.min())
    # Validate the entire physical bond, using compressed factors.
    compressed_bond = sum(sp.kron((v.T@sp.kron(sp.kron(sp.eye(q), b), sp.eye(e))@v),
                                 (v.T@sp.kron(sp.kron(a, sp.eye(q)), sp.eye(e))@v))
                          for a, b in zip(left, right))
    target_bond = sum(sp.kron(r, g) for r, g in zip(rs, gs))/(2*n)**2
    bond_err = float(sp.linalg.norm(compressed_bond-target_bond))
    assert max(errs) < 2e-13 and bond_err < 2e-12
    assert min_positive_projector_entry > 0
    result = dict(q=q, n=n, graded=graded, logical_dimension=d,
                  physical_dimension=v.shape[0], alpha=1/(2*n),
                  max_operator_error=max(errs), bond_frobenius_error=bond_err,
                  projector_min_offdiagonal=min_positive_projector_entry)
    if graded:
        omega = np.array([1 if a < q//2 else -1 for a in range(q)])
        env_parity = np.zeros(e)
        for side in range(2):
            for f in range(n):
                for x in range(d):
                    env = (side*n+f)*(d//2)+erased_bit(x, n, f)
                    val = (-1)**sum(bits(x, n))*omega[port_value(x, n, f, side, pairs)]*omega[0]
                    if env_parity[env] != 0:
                        assert env_parity[env] == val
                    env_parity[env] = val
        phys_parity = sp.diags(np.kron(np.kron(omega, omega), env_parity))
        result['grading_error'] = opnorm(phys_parity@v-v@parity)
        assert result['grading_error'] < 1e-13
    return result


def history_check(q, n):
    """Explicit full permutations, global clock, and initialization penalty.

    These are deliberately NOT claimed to be a nearest-neighbor realization.
    """
    pairs, left, right = pairs_and_ops(q, n)
    d, qp = 2**n, q+1
    k = qp*qp*d
    base = np.array([(q*qp+q)*d+x for x in range(d)])
    perms = []
    history_rows = [base]
    for side in range(2):
        for f in range(n):
            read = []
            for x in range(d):
                val = port_value(x, n, f, side, pairs)
                l, r = (val, q) if side == 0 else (q, val)
                read.append((l*qp+r)*d+erased_bit(x, n, f))
            read = np.array(read)
            per = np.arange(k)
            per[base] = read; per[read] = base
            u = sp.csr_matrix((np.ones(k), (per, np.arange(k))), shape=(k, k))
            perms.extend([u, u])
            history_rows.extend([read, base])
    c = len(history_rows)
    rows = np.concatenate([t*k+r for t, r in enumerate(history_rows)])
    cols = np.tile(np.arange(d), c)
    v = sp.csr_matrix((np.full(c*d, 1/np.sqrt(c)), (rows, cols)), shape=(c*k, d))
    blocks = [[None]*c for _ in range(c)]
    for t in range(c):
        diag = 0.5 if t in (0, c-1) else 1
        blocks[t][t] = diag*sp.eye(k)
    init = np.ones(k); init[base] = 0
    blocks[0][0] += sp.diags(init)
    for t, u in enumerate(perms):
        blocks[t+1][t] = -u/2
        blocks[t][t+1] = -u.T/2
    h = sp.bmat(blocks, format='csr')
    gs, rs, _ = cliffords(n)
    errs = [opnorm(v.T@v-sp.eye(d)), float(sp.linalg.norm(h@v))]
    for operators, targets, side in [(left, gs, 0), (right, rs, 1)]:
        for f, target in zip(operators, targets):
            f = sp.block_diag((f, sp.csr_matrix((1, 1))))
            phys = sp.kron(sp.kron(f, sp.eye(qp)), sp.eye(d)) if side == 0 else \
                   sp.kron(sp.kron(sp.eye(qp), f), sp.eye(d))
            errs.append(opnorm(v.T@sp.kron(sp.eye(c), phys)@v-target/c))
    spectrum = eigvalsh(h.toarray())
    theoretical = np.concatenate([
        np.repeat(1-np.cos(np.pi*np.arange(c)/c), d),
        np.repeat(1-np.cos(np.pi*(np.arange(c)+0.5)/c), k-d)])
    spectrum_error = float(np.max(np.abs(spectrum-np.sort(theoretical))))
    off = h-sp.diags(h.diagonal())
    assert off.data.max() <= 0
    assert max(errs) < 2e-13 and spectrum_error < 5e-13
    return dict(q=q, n=n, clock_length=c, total_dimension=c*k,
                alpha=1/c, max_code_error=max(errs), spectrum_error=spectrum_error,
                measured_ground_degeneracy=int(np.sum(abs(spectrum)<1e-10)),
                gap=float(spectrum[d]), predicted_gap=float(1-np.cos(np.pi/(2*c))))


def colored_readout_check(q, n):
    assert q % 2 == 0
    s, d = q//2, 2**n
    pairs, left, right = pairs_and_ops(q, n, True)
    e = n*d
    dims = (q, s+1, e, s+1, q)
    rows, cols, coordinates = [], [], []
    for x in range(d):
        for side in range(2):
            for f in range(n):
                bs = bits(x,n)
                neg = sum(bs[:f] if side == 0 else bs[f+1:]) % 2
                a,b = pairs[f,neg]
                color = (a+b-s) % s
                val = port_value(x,n,f,side,pairs)
                env = (side*n+f)*(d//2)+erased_bit(x,n,f)
                coord = (val,color,env,s,0) if side == 0 else (0,s,env,color,val)
                rows.append(np.ravel_multi_index(coord,dims)); cols.append(x)
                coordinates.append(coord)
    v = sp.csr_matrix((np.full(len(rows),1/np.sqrt(2*n)),(rows,cols)),
                      shape=(int(np.prod(dims)),d))
    gs,rs,_ = cliffords(n)
    port_error = 0.
    rest = v.shape[0]//q
    for a,b,g,r in zip(left,right,gs,rs):
        port_error = max(port_error,
            opnorm(v.T@sp.kron(a,sp.eye(rest))@v-g/(2*n)),
            opnorm(v.T@sp.kron(sp.eye(rest),b)@v-r/(2*n)))
    masks = []
    density_error = 0.
    for f in range(n):
        mask = np.zeros((q,s+1))
        for a,b in pairs[f]:
            color = (a+b-s) % s
            mask[a,color],mask[b,color] = 1.,-1.
        masks.append(mask)
        values = np.array([mask[c[0],c[1]] for c in coordinates])
        acted = sp.csr_matrix((values/np.sqrt(2*n),(rows,cols)),shape=v.shape)
        z = sp.diags([(-1)**bits(x,n)[f] for x in range(d)],dtype=float)
        density_error=max(density_error,opnorm(v.T@acted-z/(2*n)))
    overlap = max(float(np.max(abs(a*b))) for i,a in enumerate(masks)
                  for j,b in enumerate(masks) if i!=j)
    mu = np.linspace(-.7,.9,n)
    combined=sum(t*a for t,a in zip(mu,masks))
    assert port_error<1e-13 and density_error<1e-13 and overlap==0
    assert abs(np.max(abs(combined))-np.max(abs(mu)))<1e-13
    return dict(q=q,n=n,physical_dimension=v.shape[0],neighbor_colors=s,
                max_port_error=port_error,max_density_error=density_error,
                distinct_density_support_overlap=overlap,
                diagonal_field_norm=float(np.max(abs(combined))))


def main():
    result = {'status': 'finite code checks, not CFT numerics',
              'static': [static_check(q, n, g) for q, n, g in
                         [(3, 1, False), (4, 3, False), (5, 5, False),
                          (6, 6, False), (4, 2, True), (6, 4, True)]],
              'history': [history_check(3, 1), history_check(4, 2)],
              'density_readout': [colored_readout_check(4,2),colored_readout_check(6,4)]}
    with open('classical_readout_checks.json', 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
