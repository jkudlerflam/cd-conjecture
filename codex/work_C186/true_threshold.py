"""C186: true (off-grid) dressed valley minima of the mixed model at the ring-tuned (phi, mu), from the infinite-line
few-body code (mixed_fewbody.py, N_max = 5, box R1 = 26; the dressing cloud is compact), so that the ring pair level
can be compared with the infinite-volume threshold rather than the grid one."""
import numpy as np, json, sys
import mixed_fewbody as mf
for f in sys.argv[1:]:
    rows = [json.loads(l) for l in open(f)]
    d = rows[-1]
    m = mf.Model(d['u'], d['phi'], d['v'])
    m.mu = d['mu']
    ka0, kb0 = mf.bare_valleys(m)[:2]
    out = dict(v=d['v'], Nmax_pair=d['Nmax'], Epair=d['Epair'])
    for name, k0, kgrid in (("a", ka0, None), ("b", kb0, None)):
        k, e = mf.dressed_min(m, 5, 26, k0, 0.15)
        out["k" + name] = float(k); out["e" + name] = float(e)
    # grid values with the same few-body code, for consistency with the ring E1
    L = d['L']
    out["ea_grid_fb"] = float(m.lowest(2*np.pi*round(out["ka"] % (2*np.pi) * L/(2*np.pi))/L, 5, 26)[0])
    out["T_true"] = out["ea"] + out["eb"]
    out["B_true"] = out["T_true"] - d['Epair']
    out["E1a_ring"] = d['E1a']; out["E1b_ring"] = d['E1b']
    print(json.dumps(out), flush=True)
