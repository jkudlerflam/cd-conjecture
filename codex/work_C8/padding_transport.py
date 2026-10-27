"""Two-length payload compatibility, without a local-parent claim."""
import json
import numpy as np
from scipy.sparse import csr_matrix, kron
from c8_resume_checks import joint_code


def embed_and_pad(V, q, k, p):
    """One vacancy symbol, and p vacancies coherently on either side."""
    r = q+1
    D = V.shape[1]
    base_ids = []
    for i in range(q**k):
        digits = np.unravel_index(i, (q,)*k)
        base_ids.append(np.ravel_multi_index(tuple(x+1 for x in digits), (r,)*k))
    base_ids = np.asarray(base_ids)
    row = np.concatenate([np.repeat(base_ids*r**p, D), np.repeat(base_ids, D)])
    col = np.tile(np.arange(D), 2*len(base_ids))
    data = np.tile(V.reshape(-1), 2)/np.sqrt(2)
    return csr_matrix((data, (row, col)), shape=(r**(k+p), D))


def endpoint_compression(W, F, r, length, side):
    v = W.toarray().reshape((r,)*length+(W.shape[1],))
    j = 0 if side == 'left' else length-1
    fv = np.tensordot(F, v, axes=(1, j))
    fv = np.moveaxis(fv, 0, j)
    return W.conj().T @ fv.reshape(W.shape)


def run():
    q = 2; m = 2; k = 5
    fs, gs, rs, parity, alpha, V, _, _ = joint_code(q, m)
    D = len(parity)
    # L,R,E1,E2,E3 -> L,E1,E2,E3,R, so the ports are the endpoints.
    V = V.reshape((2,)*5+(D,)).transpose(0, 2, 3, 4, 1, 5).reshape(2**5, D)
    W = [embed_and_pad(V, q, k, p) for p in [1, 2]]
    err = 0.
    iso = 0.
    for p, w in zip([1, 2], W):
        iso = max(iso, np.linalg.norm((w.conj().T@w).toarray()-np.eye(D)))
        for f, g, rr in zip(fs, gs, rs):
            F = np.zeros((q+1, q+1), complex); F[1:, 1:] = f
            err = max(err, np.linalg.norm(endpoint_compression(w, F, q+1, k+p, 'left')-alpha*g/2))
            err = max(err, np.linalg.norm(endpoint_compression(w, F, q+1, k+p, 'right')-alpha*rr/2))
    before = kron(W[0], W[1], format='csr')
    after = kron(W[1], W[0], format='csr')
    projected = (after.conj().T@before).toarray()
    transport_err = np.linalg.norm(projected-np.eye(D**2)/4)
    # Endpoint bond is independent of the two choices of cell length.
    bond_err = 0.
    target = alpha**2/4*sum(np.kron(a, b) for a, b in zip(rs, gs))
    for p, v in zip([1, 2], W):
        for pp, w in zip([1, 2], W):
            bond = np.zeros((D**2, D**2), complex)
            for f in fs:
                F = np.zeros((q+1, q+1), complex); F[1:, 1:] = f
                a = endpoint_compression(v, F, q+1, k+p, 'right')
                b = endpoint_compression(w, F, q+1, k+pp, 'left')
                bond += np.kron(a, b)
            bond_err = max(bond_err, np.linalg.norm(bond-target))
    assert max(iso, err, transport_err, bond_err) < 1e-10
    out = dict(core_length=k, port_dimension=q, padded_site_dimension=q+1,
               logical_dimension=D, lengths=[k+1, k+2],
               endpoint_scale=alpha/2, isometry_error=float(iso),
               endpoint_compression_error=float(err),
               transport_amplitude=.25, transport_error=float(transport_err),
               length_independent_bond_error=float(bond_err),
               local_parent_constructed=False)
    with open('padding_transport_checks.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    run()
