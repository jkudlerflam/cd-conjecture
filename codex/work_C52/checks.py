"""C52 finite algebra controls. No microscopic IR or C42-margin computation."""
import json
from pathlib import Path
import numpy as np


def fermions(n):
    size = 1 << n
    cs = []
    for j in range(n):
        c = np.zeros((size, size), complex)
        for b in range(size):
            if (b >> j) & 1:
                c[b ^ (1 << j), b] = (-1) ** ((b & ((1 << j)-1)).bit_count())
        cs.append(c)
    return cs


def schmidt_rank(h, dl, dr):
    reshuffled = h.reshape(dl, dr, dl, dr).transpose(0, 2, 1, 3).reshape(dl*dl, dr*dr)
    return int(np.linalg.matrix_rank(reshuffled, tol=1e-10))


def transfer_matrix(h, n):
    d = 1 << n
    pairs = [(1 << j, 0) for j in range(n)] + [(0, 1 << j) for j in range(n)]
    return np.array([[h[a*d+c, b*d+e] for c, e in pairs] for a, b in pairs])


def hopping_controls():
    rng = np.random.default_rng(52)
    results = []
    for n in (1, 2, 3):
        cs = fermions(2*n)
        d = 1 << n
        raw = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
        t = raw + 3*np.eye(n)
        h = np.zeros((d*d, d*d), complex)
        for a in range(n):
            for b in range(n):
                term = t[a,b] * cs[a].conj().T @ cs[n+b]
                h += term + term.conj().T
        # Site-local terms and normal-ordered density quartics do not affect transfer.
        h_extra = h.copy()
        ns = [c.conj().T @ c for c in cs]
        for a in range(n):
            h_extra += .17*(a+1)*ns[a] + .23*ns[a]@ns[n+a]
        block = transfer_matrix(h, n)
        extra_error = np.linalg.norm(block-transfer_matrix(h_extra, n))
        assert extra_error < 1e-12
        assert np.linalg.matrix_rank(block, tol=1e-10) == 2*n
        assert schmidt_rank(h, d, d) == 2*n
        results.append(dict(n=n, hopping_rank=int(np.linalg.matrix_rank(t)),
                            transfer_rank=2*n, hopping_operator_schmidt_rank=schmidt_rank(h,d,d),
                            quartic_and_onsite_extraction_error=float(extra_error)))
    return results


def boundary_control():
    # Two logical fermions occupy each two-qubit cell. Rotate the one-particle
    # code sector, but preserve charge. One physical cut still carries rank one.
    angle = .43
    g = np.eye(4)
    g[np.ix_([1,2],[1,2])] = [[np.cos(angle), -np.sin(angle)],
                              [np.sin(angle), np.cos(angle)]]
    cs = fermions(4)
    h = cs[1].conj().T @ cs[2]
    h = h + h.conj().T
    gg = np.kron(g,g)
    hc = gg.conj().T @ h @ gg
    # Ordering of little-endian fermion basis and tensor cells exchanges labels,
    # which preserves the rank under test.
    t = hc[np.ix_([4,8],[1,2])]
    assert np.linalg.matrix_rank(t, tol=1e-10) == 1
    assert schmidt_rank(hc,4,4) == 2
    # A rank-one 3-band hopping has determinant Fourier degree <= 1.
    rng = np.random.default_rng(521)
    a = rng.normal(size=(3,3)); a = a+a.T+np.diag([-.8,.2,1.3])
    u = rng.normal(size=3)+1j*rng.normal(size=3)
    v = rng.normal(size=3)+1j*rng.normal(size=3)
    tr = np.outer(u,v.conj())
    kk = 2*np.pi*np.arange(64)/64
    vals = np.array([np.linalg.det(a+np.exp(1j*k)*tr+np.exp(-1j*k)*tr.conj().T) for k in kk])
    coeff = np.fft.fft(vals)/len(vals)
    tail = max(abs(coeff[2:-1]))
    assert tail < 1e-10
    return dict(encoded_hopping_rank=1, encoded_boundary_schmidt_rank=2,
                rank_one_determinant_forbidden_fourier_error=float(tail))


def adjacent_commutator_controls():
    rng = np.random.default_rng(522)
    rows = []
    for n in (1,2):
        d = 1 << n
        t = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))+2*np.eye(n)
        h = np.zeros((d*d,d*d),complex)
        for a in range(n):
            for b in range(n):
                h[(1<<a)*d,1<<b] = t[a,b]
                h[1<<b,(1<<a)*d] = np.conj(t[a,b])
        # Arbitrary number-conserving changes outside the vacuum/one-particle
        # sectors test that the argument does not require Gaussian higher sectors.
        for number in range(2,2*n+1):
            ids = [i for i in range(d*d) if i.bit_count()==number]
            raw = rng.normal(size=(len(ids),len(ids)))+1j*rng.normal(size=(len(ids),len(ids)))
            h[np.ix_(ids,ids)] = raw+raw.conj().T
        tensor = h.reshape(d,d,d,d)
        left = np.einsum('abcb->ac',tensor)/d
        right = np.einsum('abad->bd',tensor)/d
        connected = h-np.kron(left,np.eye(d))-np.kron(np.eye(d),right)+np.trace(h)/(d*d)*np.eye(d*d)
        k12 = np.kron(connected,np.eye(d)); k23 = np.kron(np.eye(d),connected)
        comm = k12@k23-k23@k12
        outer = comm[np.ix_([(1<<i)*d*d for i in range(n)], [1<<i for i in range(n)])]
        error = np.linalg.norm(outer-t@t)
        assert error < 1e-12
        rows.append(dict(n=n, outer_commutator_identity_error=float(error),
                         outer_commutator_norm=float(np.linalg.norm(outer,2))))
    return rows


def head_controls():
    rows = []
    for b in (4,5,8,16):
        short = np.sqrt(2/b)*np.sin(np.pi*np.arange(1,b)/b)
        long = np.sqrt(2/(b+1))*np.sin(np.pi*np.arange(1,b+1)/(b+1))
        grow = np.r_[short,0.]
        shrink = long[1:]
        out = np.kron(grow,shrink)
        target = np.kron(long,short)
        overlap = float(target @ out)
        a = float(short @ long[:-1])
        s = 1-float(long[0]**2)
        leakage = float(np.linalg.norm(out-overlap*target)**2)
        formula = s-a**4
        lower = s*(1-s)
        assert abs(overlap-a*a) < 2e-14
        assert abs(leakage-formula) < 2e-14
        assert leakage >= lower-2e-14 and lower > 0
        rows.append(dict(b=b, projected_amplitude=overlap, output_norm_squared=float(out@out),
                         head_leakage_squared=leakage, lower_bound=lower,
                         formula_error=abs(leakage-formula)))
    return rows


def resolvent_control():
    # A rank-one source can produce arbitrarily high Schmidt rank after a
    # non-product resolvent. A Cauchy matrix gives a finite check of this fact.
    rows = []
    for n in (2,3,4):
        aa = np.arange(1,n+1,dtype=float)
        bb = np.arange(1,n+1,dtype=float)+.3
        coefficients = 1/(aa[:,None]+bb[None,:])
        h = np.diag(coefficients.ravel())
        rank = schmidt_rank(h,n,n)
        assert rank == n
        rows.append(dict(cell_code_dimension=n, effective_operator_schmidt_rank=rank))
    return rows


if __name__ == '__main__':
    result = dict(scope='Finite algebra controls only; no central charge or microscopic matching claim.',
                  hopping=hopping_controls(), boundary=boundary_control(),
                  marker_heads=head_controls(), adjacent_commutators=adjacent_commutator_controls(),
                  virtual_resolvent=resolvent_control())
    Path(__file__).with_name('checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
