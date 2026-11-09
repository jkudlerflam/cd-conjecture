"""C176 validation, part 3: pair correlations against ED (c frame, fermion operators with explicit signs)."""
import numpy as np
from fibdmrg import ground
from validate import ed
from ring_ed import popcount
def apply_cc(b, psi, i, j, dag):
    """(c_i c_j) psi if not dag, (c_i c_j)^dag psi = c_j^+ c_i^+ psi if dag; JW from site 0; returns dict state->amp."""
    out = {}
    for s, a in zip(b, psi):
        s = int(s); amp = a
        ops = [('d', i), ('d', j)] if dag else [('a', j), ('a', i)]   # rightmost acts first
        # dag: c_j^+ c_i^+ : first c_i^+, then c_j^+ ; not dag: c_i c_j : first c_j, then c_i
        seq = [('c', i), ('c', j)] if dag else [('a', j), ('a', i)]
        ok = True
        for (t, site) in seq:
            occ = (s >> site) & 1
            sign = (-1) ** popcount(np.array([s & ((1 << site) - 1)]))[0]
            if t == 'a':
                if not occ: ok = False; break
                s ^= (1 << site); amp *= sign
            else:
                if occ: ok = False; break
                s ^= (1 << site); amp *= sign
        if ok:
            out[s] = out.get(s, 0) + amp
    return out
for (L, N, phi, g) in ((20, 4, np.pi/2, 1), (18, 4, 1.45, 0)):
    w = np.exp(1j*phi)
    b, e0, psi, e1 = ed(L, N, w)
    dm, hist = ground(L, N, w, chis=[16, 32, 64, 128], ch_kw=dict(gauge=g), tol=1e-12)
    i0 = L//2 - 3
    m = dm.measure(gwin=(0, L), refs=[i0])
    for dd in (2, 3):
        row = np.array(m['Dpair'][f"{dd}_{i0}"]); Dm = row[:, 0] + 1j*row[:, 1]
        ref = []
        for jj in range(i0 + dd + 1, L - dd):
            A = apply_cc(b, psi, jj, jj + dd, False)          # c_j c_{j+d} |psi>
            B = apply_cc(b, psi, i0, i0 + dd, False)          # c_i0 c_{i0+d} |psi>
            ref.append(sum(np.conj(B[k]) * A[k] for k in A if k in B))
        ref = np.array(ref)
        print(f"L={L} N={N} phi={phi:.3f} d={dd}: max|D_dmrg - D_ed| = {np.max(np.abs(Dm - ref)):.1e}, max|D| = {np.max(np.abs(ref)):.2e}")
