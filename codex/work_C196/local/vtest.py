import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import idmrg196 as I, vumps196 as Vm, validate196 as V
cfg = eval(sys.argv[1]); D = int(sys.argv[2]); nit = int(sys.argv[3])
mpo, W, eo = I.build_model(cfg); ex = I.exact_energy(cfg)
v = Vm.VUMPS(W, mpo, D=D, seed=0, cmode=(sys.argv[4] if len(sys.argv) > 4 else 'eig'))
t0 = time.time()
for k in range(nit):
    r = v.step()
    if k % max(1, nit // 25) == 0 or k == nit - 1:
        de = (r['e'] + eo - ex['e']) if ex else float('nan')
        print(k, '%.12f' % (r['e'] + eo), '%.2e' % de, 'eR-e %.1e' % (r['eR'] - r['e']), 'eps %.1e' % r['eps'], 'nA', r['nmvA'], 'nC', r['nmvC'], 'gm', r['gmvL'], r['gmvR'], '%.0e' % r['gmL'], 'n %.5f' % r['n'], 'S %.5f' % r['S'], '%.0fs' % (time.time() - t0), flush=True)
x = v.transfer_spectrum(); print('lam0', x['lam0'], 'xi_e', x['xi_even'], 'xi_o', x['xi_odd'], ex)
