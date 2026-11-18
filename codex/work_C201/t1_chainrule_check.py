"""Check of the decoupling estimate (14) in the refereed form of Theorem A, on ring ground states, with
single-site windows (t=1, q=2) and the octahedral rank-one POVM {P_k/3}, dual X_k = 3 P_k - 1 (||X||_1 = 2q-1 = 3).
For each j:  || rho_{U W_j} - sum_v p_v w_U^v (x) w_{W_j}^v ||_1  <=  (2q-1) sqrt(2 I(U:Y_j|Y_<j)),
and sum_j I(U:Y_j|Y_<j) <= S(U).   v = outcomes on W_1..W_{j-1}.
Claude Fable 5.2 checking Astra gpt-6.2-astra, work_C201."""
import numpy as np, itertools
from t1_lemma_checks import ring_gs, rdm, tn

def vn(r):
    w = np.linalg.eigvalsh(r); w = w[w > 1e-14]; return float(-(w*np.log(w)).sum())
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1., -1]).astype(complex)
projs = [(np.eye(2)+s*P)/2 for P in (X, Y, Z) for s in (1, -1)]
# reconstruction identity
g = np.array([[0.7, 0.2-0.1j], [0.2+0.1j, 0.3]])
assert np.allclose(sum(np.trace(g@(P/3))*(3*P-np.eye(2)) for P in projs), g)
for kind in ('heis', 'xx', 'tfim', 'ghz'):
    L = 12; psi = ring_gs(L, kind); K = 4
    SU = vn(rdm(psi, L, [0])); tot = 0
    for j in range(1, K+1):
        pre = list(range(1, j)); Wj = j
        rho = rdm(psi, L, [0, Wj]+pre).reshape([2]*(2*(j+1)))    # U, W_j, prefix
        full = rdm(psi, L, [0, Wj]+pre)
        dimp = 2**(j-1)
        eta = np.zeros((4, 4), complex); Icond = 0.0
        for outs in itertools.product(range(6), repeat=j-1):
            F = np.array([[1.0]])
            for o in outs: F = np.kron(F, projs[o]/3)
            blk = np.einsum('aibj,ji->ab', full.reshape(4, dimp, 4, dimp), F) if j > 1 else full
            p = np.trace(blk).real
            if p < 1e-14: continue
            w = blk/p
            wU = np.einsum('abcb->ac', w.reshape(2, 2, 2, 2)); wW = np.einsum('abac->bc', w.reshape(2, 2, 2, 2))
            eta += p*np.kron(wU, wW)
            # I(U:Y_j) in the conditional state: S(U) - sum_y p_y S(U|y)
            Iy = vn(wU)
            for P in projs:
                c = np.einsum('abcd,db->ac', w.reshape(2, 2, 2, 2), P/3); py = np.trace(c).real
                if py > 1e-14: Iy -= py*vn(c/py)
            Icond += p*Iy
        lhs = tn(rdm(psi, L, [0, Wj])-eta); rhs = 3*np.sqrt(2*max(Icond, 0)); tot += Icond
        print(f"  {kind:5s} j={j}: ||rho_UW - eta||_1 = {lhs:.5f} <= 3 sqrt(2 I) = {rhs:.5f}   I(U:Y_j|Y_<j) = {Icond:.5f}   {'ok' if lhs <= rhs+1e-9 else 'FAIL'}")
    print(f"  {kind:5s} sum_j I = {tot:.5f} <= S(U) = {SU:.5f}: {tot <= SU+1e-9}")
