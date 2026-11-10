import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import idmrg196 as I, vumps196 as Vm
cfg = eval(sys.argv[1]); Ds = [int(x) for x in sys.argv[2].split(',')]; nit = int(sys.argv[3])
mpo, W, eo = I.build_model(cfg); ex = I.exact_energy(cfg)
v = Vm.VUMPS(W, mpo, D=Ds[0], seed=0)
t0 = time.time(); pts = []
for D in Ds:
    v.expand(D)
    for k in range(nit):
        r = v.step()
        if r['eps'] < 1e-9 and k > 10: break
    x = v.transfer_spectrum(k=6)
    de = (r['e'] + eo - ex['e']) if ex else float('nan')
    pts.append((D, r['S'], x['xi_even'], x['xi_odd']))
    so = x['spec'][1]['abs']; se = x['spec'][0]['abs']
    print(D, 'it', k, 'e %.12f' % (r['e'] + eo), 'de %.1e' % de, 'eps %.1e' % r['eps'], 'n %.5f' % r['n'], 'S %.5f' % r['S'],
          'xi_e %.2f xi_o %.3g' % (x['xi_even'], x['xi_odd']), 'even', np.round(se[:4], 8).tolist(), 'odd', np.round(so[:4], 8).tolist(), '%.0fs' % (time.time() - t0), flush=True)
pts = np.array(pts)
for col, name in ((2, 'xi_even'),):
    print('local c from', name, np.round(6 * np.diff(pts[:, 1]) / np.diff(np.log(pts[:, col])), 3).tolist())
