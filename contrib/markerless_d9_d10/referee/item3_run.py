"""Item 3 driver: literal second-order logical Hamiltonian (second_order.py) against
 (i) the author's scalar9.blocks, (ii) the referee's reading of the README formula.
Also t_d = 0, alternative paddings, and periodic rings."""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(1, os.path.dirname(HERE))
import mpmath as mp
import ref_bonds as rb
import second_order as so
import readme_scalar as rs
import scalar9

mp.mp.dps = 40
def c(M): return np.array(M.tolist(), dtype=complex)
def al4(n): return (2 / (n + 2)) ** 2 * np.sin(np.pi / (n + 2)) ** 4
res = []
def rec(d):
    res.append(d); print(json.dumps(d), flush=True)

pts = [(2, 2.6, 3.0, 0.07), (2, 2.2, 0.4, 0.3), (3, 2.6, 3.0, 0.07), (3, 2.2, 0.4, 0.3)]
if len(sys.argv) > 1 and sys.argv[1] == "n2":
    pts = pts[:2]
store = {}
for (n, a, D, td) in pts:
    for tdv in (td, 0.0):
        t0 = time.time()
        out, h, e, s = so.analyse(n, a, D, tdv)
        T_bf, muL_bf, muR_bf = h[:n, n:], h[:n, :n], h[n:, n:]
        store[(n, a, D, tdv)] = (T_bf, muL_bf, muR_bf)
        Ta, mLa, mRa = [c(x) * al4(n) for x in scalar9.blocks(n, a=mp.mpf(a), D=mp.mpf(D), td=mp.mpf(tdv))]
        Tr, mLr, mRr = [x * al4(n) for x in rs.blocks_np(n, a, D, tdv)]
        sc = abs(T_bf).max()
        d = dict(test="padded two cells [R0|cell|cell|L0]", **{k: (str(v) if isinstance(v, complex) else v) for k, v in out.items()})
        d.update(T_scale=float(sc), muL_over_T=float(abs(muL_bf).max() / sc), muR_over_T=float(abs(muR_bf).max() / sc),
                 mu_sum_over_T=float(abs(muL_bf + muR_bf).max() / sc),
                 gauss_resid_over_T=out["gauss_resid"] / sc, M0_bitdep_over_M0=out["M0_bit_dependence"] / abs(out["M0_const"]),
                 M1_over_T=out["M1_max"] / sc,
                 dT_author=float(abs(T_bf - Ta).max() / sc), dmuL_author=float(abs(muL_bf - mLa).max() / sc),
                 dmuR_author=float(abs(muR_bf - mRa).max() / sc),
                 dT_readme=float(abs(T_bf - Tr).max() / sc), dmuL_readme=float(abs(muL_bf - mLr).max() / sc),
                 dmuR_readme=float(abs(muR_bf - mRr).max() / sc),
                 seconds=round(time.time() - t0, 1))
        rec(d)

# alternative paddings at n = 2 (the neighbour letter of a real ring can be R_1, or a head)
n, a, D, td = 2, 2.6, 3.0, 0.07
ref = store[(n, a, D, td)]
# (a head is not a frozen letter, so head paddings are not a valid test; the ring test covers that case)
for (pl, pr, name) in [(rb.R1, rb.L1, "R1|..|L1"), (rb.R1, rb.L0, "R1|..|L0"), (rb.R0, rb.L1, "R0|..|L1")]:
    out, h, e, s = so.analyse(n, a, D, td, padL=pl, padR=pr)
    sc = abs(ref[0]).max()
    rec(dict(test="padding " + name, dim=out["dim"], F_max=out["F_max"],
             dT=float(abs(h[:n, n:] - ref[0]).max() / sc), dmuL=float(abs(h[:n, :n] - ref[1]).max() / sc),
             dmuR=float(abs(h[n:, n:] - ref[2]).max() / sc), gauss_resid_over_T=out["gauss_resid"] / sc,
             code_resid=out["code_resid"]))

json.dump(res, open(os.path.join(HERE, "item3_padded.json"), "w"), indent=1)
