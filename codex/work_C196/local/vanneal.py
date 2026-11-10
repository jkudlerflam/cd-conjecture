import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import idmrg196 as I, vumps196 as Vm, models as Mo
b = Mo.bond_family('p153')
mus = [float(x) for x in sys.argv[1].split(',')]; D = int(sys.argv[2]); nit = int(sys.argv[3])
v = None; t0 = time.time()
for mu in mus:
    hz = Mo.hz_for_mu_eff(b, mu)
    mpo, W, eo = I.build_model(dict(model='p153', hz=hz))
    if v is None:
        v = Vm.VUMPS(W, mpo, D=D, seed=0)
    else:
        v.Wh = np.ascontiguousarray(W); v.mpo = mpo; v.Le = v.Re = None
    for k in range(nit):
        r = v.step()
        if r['eps'] < 1e-8 and k > 10: break
    print('mu_eff %.4f hz %.8f it %d e %.10f eps %.1e (last5 max %s) n %.4f S %.4f nA %d nC %d %.0fs' % (mu, hz, k, r['e'] + eo, r['eps'], '', r['n'], r['S'], r['nmvA'], r['nmvC'], time.time() - t0), flush=True)
x = v.transfer_spectrum(); print('xi_e', x['xi_even'], 'xi_o', x['xi_odd'], np.round(x['spec'][1]['abs'][:3], 6).tolist())
