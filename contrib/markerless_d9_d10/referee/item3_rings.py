"""Item 3, faithfulness of the padding: genuinely periodic rings built literally.

ring3: three cells on a ring.  The one-particle block of the lam^2 g^2 logical
Hamiltonian must be block-circulant with diagonal block mu^L + mu^R and
nearest-cell block T taken from the padded two-cell computation; in the
N-particle sector the seam hop carries the Jordan-Wigner sign (-1)^(N-1).
ring2: two cells on a ring (shown only to exhibit the wrap-around artefact).
"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import second_order as so

res = []
def rec(d):
    res.append(d); print(json.dumps(d), flush=True)

def gauss_twisted(M2, h_open, wrap_blocks, nm, n):
    """residual of M2 against const + dGamma(h) with seam sign (-1)^(N-1) per particle-number sector"""
    c = so.fermion_ops(nm)
    dim = 2 ** nm
    const = M2[0, 0]
    Nop = sum(cc.T @ cc for cc in c)
    par = np.diag(np.exp(1j * np.pi * (np.diag(Nop) - 1))).real      # (-1)^(N-1), evaluated after the hop (N conserved)
    G = const * np.eye(dim, dtype=complex)
    for m in range(nm):
        for mp_ in range(nm):
            if h_open[m, mp_] != 0:
                G += h_open[m, mp_] * (c[m].T @ c[mp_])
            if wrap_blocks[m, mp_] != 0:
                G += wrap_blocks[m, mp_] * (c[m].T @ c[mp_]) @ par
    return float(abs(M2 - G).max())

for (n, a, D, td) in [(1, 2.6, 3.0, 0.07), (1, 2.2, 0.4, 0.3), (2, 2.6, 3.0, 0.07), (2, 2.2, 0.4, 0.3)]:
    t0 = time.time()
    o2, h2, e2, s2 = so.analyse(n, a, D, td)                       # padded two cells
    T, muL, muR = h2[:n, n:], h2[:n, :n], h2[n:, n:]
    sc = abs(T).max()
    o3, h3, e3, s3 = so.analyse(n, a, D, td, geometry="ring3", ncell=3)
    pred = np.zeros((3 * n, 3 * n), complex)
    for x in range(3):
        y = (x + 1) % 3
        pred[x*n:(x+1)*n, x*n:(x+1)*n] += muL + muR
        pred[x*n:(x+1)*n, y*n:(y+1)*n] += T
        pred[y*n:(y+1)*n, x*n:(x+1)*n] += T.conj().T
    # split the prediction into the part with ordinary JW strings and the seam part
    seam = np.zeros_like(pred)
    seam[2*n:3*n, 0:n] = T; seam[0:n, 2*n:3*n] = T.conj().T
    openpart = pred - seam
    c3 = e3["M2"][0, 0]
    rec(dict(test="ring of three cells vs padded two cells", n=n, a=a, D=D, td=td, dim_ring=o3["dim"], F_max=o3["F_max"],
             first_order=o3["first_order"], code_resid=o3["code_resid"],
             M0_ring_minus_3xM0_pad=float(abs(o3["M0_const"] - 3 * o2["M0_const"])), M0_bitdep_over_M0=o3["M0_bit_dependence"] / abs(o3["M0_const"]),
             M1_max=o3["M1_max"],
             const_ring_minus_3xconst_pad_over_T=float(abs(c3 - 3 * e2["M2"][0, 0]) / sc),
             one_particle_block_error_over_T=float(abs(h3 - pred).max() / sc),
             gauss_resid_over_T_with_seam_sign=gauss_twisted(e3["M2"], openpart, seam, 3 * n, n) / sc,
             gauss_resid_over_T_naive_no_seam_sign=o3["gauss_resid"] / sc,
             T_scale=float(sc), seconds=round(time.time() - t0, 1)))

# two-cell ring: wrap-around artefact
for (n, a, D, td) in [(2, 2.6, 3.0, 0.07), (2, 2.2, 0.4, 0.3)]:
    o2, h2, e2, s2 = so.analyse(n, a, D, td)
    T, muL, muR = h2[:n, n:], h2[:n, :n], h2[n:, n:]
    sc = abs(T).max()
    oR, hR, eR, sR = so.analyse(n, a, D, td, geometry="ring2", ncell=2)
    naive = np.block([[muL + muR, T + T.conj().T], [T + T.conj().T, muL + muR]])
    rec(dict(test="ring of two cells (artefact expected)", n=n, a=a, D=D, td=td, dim=oR["dim"],
             first_order=oR["first_order"], M0_bitdep_over_M0=oR["M0_bit_dependence"] / abs(oR["M0_const"]),
             M0_ring=str(oR["M0_const"]), two_x_M0_pad=str(2 * o2["M0_const"]),
             one_particle_block_deviation_from_naive_over_T=float(abs(hR - naive).max() / sc),
             gauss_resid_over_T=oR["gauss_resid"] / sc))
json.dump(res, open(os.path.join(HERE, "item3_rings.json"), "w"), indent=1)
