"""C164: print the one-particle band of a bond at several ring sizes (lattice values, interpolated valleys,
second level, two-particle floor check)."""
import json, sys, numpy as np, qband, est, qopt
def show(v, Ls=(12, 14, 16)):
    v = np.asarray(v, float); v /= np.linalg.norm(v)
    J, f = est.J_f(v); h = est.bond_from_JB(J, f)
    for L in Ls:
        b = qband.band(h, L, nev=3)
        ks, Ek, iA, iB, bar = qopt.valleys(b['E1'], 1440)
        iso = qopt.iso_discrete(b['E1'])
        s = f"L={L} m0={b['m0']} A=({ks[iA]:.3f},{Ek[iA]:.4f})"
        if iB is not None:
            s += f" B=({ks[iB]:.3f},{Ek[iB]:.4f}) barrier={bar:.4f}"
        print(s, flush=True)
        print('   k/(2pi/L) E1:', ' '.join(f'{x:.3f}' for x in b['E1']))
        print('   2nd level   :', ' '.join(f'{x:.3f}' for x in b['lev'][:, 1] - b['E0']))
        print('   isolated    :', ''.join('1' if x else '0' for x in iso))
if __name__ == '__main__':
    v = json.loads(sys.argv[1]); Ls = [int(x) for x in sys.argv[2:]] or (12, 14, 16)
    show(v, Ls)
