"""Item 3 extension: literal padded two-cell model at n = 4 (the author's brute-force check stops at n = 3)."""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.dirname(HERE))
import mpmath as mp
import second_order as so, readme_scalar as rs, scalar9
mp.mp.dps = 40
res = []
for (n, a, D, td) in [(4, 2.2, 0.4, 0.3), (4, 2.6, 3.0, 0.07)]:
    t0 = time.time()
    out, h, e, s = so.analyse(n, a, D, td)
    al4 = (2 / (n + 2)) ** 2 * np.sin(np.pi / (n + 2)) ** 4
    c = lambda M: np.array(M.tolist(), dtype=complex) * al4
    Ta, mLa, mRa = [c(x) for x in scalar9.blocks(n, a=mp.mpf(a), D=mp.mpf(D), td=mp.mpf(td))]
    Tr, mLr, mRr = [x * al4 for x in rs.blocks_np(n, a, D, td)]
    T, mL, mR = h[:n, n:], h[:n, :n], h[n:, n:]; sc = abs(T).max()
    d = dict(n=n, a=a, D=D, td=td, dim=out["dim"], first_order=out["first_order"], M1_max=out["M1_max"],
             M0_bitdep_over_M0=out["M0_bit_dependence"] / abs(out["M0_const"]), gauss_resid_over_T=out["gauss_resid"] / sc,
             dT_author=float(abs(T - Ta).max() / sc), dmuL_author=float(abs(mL - mLa).max() / sc), dmuR_author=float(abs(mR - mRa).max() / sc),
             dT_readme=float(abs(T - Tr).max() / sc), dmuL_readme=float(abs(mL - mLr).max() / sc), dmuR_readme=float(abs(mR - mRr).max() / sc),
             seconds=round(time.time() - t0, 1))
    res.append(d); print(json.dumps(d), flush=True)
    json.dump(res, open(os.path.join(HERE, "item3_n4.json"), "w"), indent=1)
