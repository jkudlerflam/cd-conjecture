"""U(1) diagnostic: |<b_i b_j>| (two annihilations, b = |0><1|) and <b+_i b_j> on the window of a checkpoint."""
import sys, json, numpy as np
sys.path.insert(0, '.')
import idmrg196b as B
cfg = json.load(open(sys.argv[1])); ck = sys.argv[2]; chi = int(sys.argv[3])
mpo, W, e_off = B.build_model(cfg)
st = B.CellIDMRG(mpo, W, cfg['q'], e_off=e_off)
st.load(ck)
st.iterate(chi)
M, n = st.M, st.n
b = np.array([[0, 1], [0, 0]], dtype=complex)   # |0><1|
bd = b.T.copy(); I2 = np.eye(2)
Rs = [None] * (n + 1)
Rs[n] = np.eye(M[-1].shape[2], dtype=complex)
for j in range(n - 1, -1, -1):
    Rs[j] = sum(M[j][:, s, :] @ Rs[j + 1] @ M[j][:, s, :].conj().T for s in (0, 1))
def tr(Lx, A, O):
    return sum(O[sp, s] * A[:, sp, :].conj().T @ Lx @ A[:, s, :] for sp in (0, 1) for s in (0, 1) if O[sp, s] != 0)
nrm = np.trace(Rs[0] @ np.eye(M[0].shape[0])).real
def corr(i, Oi, Oj):
    Lx = np.eye(M[0].shape[0], dtype=complex)
    for k in range(i):
        Lx = tr(Lx, M[k], I2)
    Lx = tr(Lx, M[i], Oi)
    out = []
    for j in range(i + 1, n):
        out.append(np.sum(tr(Lx, M[j], Oj) * Rs[j + 1].T) / nrm)
        Lx = tr(Lx, M[j], I2)
    return np.array(out)
i0 = 4
bb = corr(i0, b, b); bdb = corr(i0, bd, b)
for r in (1, 2, 3, 5, 8, 12, 16, 24, 32, 40, 48, 56):
    if r - 1 < len(bb):
        print('r=%2d |<b b>|=%.3e  |<b+ b>|=%.3e' % (r, abs(bb[r - 1]), abs(bdb[r - 1])))
