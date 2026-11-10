import sys, json, time, numpy as np
sys.path.insert(0, '.')
import idmrg196b as B, vumps196b as VB
cfg = json.load(open(sys.argv[1])); ck = sys.argv[2]; nst = int(sys.argv[3])
mpo, W, e_off = B.build_model(cfg)
st = B.CellIDMRG(mpo, W, cfg['q'], e_off=e_off)
st.load(ck)
chi = int(sys.argv[4]) if len(sys.argv) > 4 else 32
r = st.iterate(chi)
print('idmrg e', r['e_site'], 'S', r['Sbar'], flush=True)
Ms, n0 = VB.cell_from_idmrg(st)
vm = VB.CellVUMPS(mpo, W, e_off=e_off)
eps = vm.from_cell(Ms, n0)
print('start eps', eps, 'D', vm.D[:4], flush=True)
for k in range(nst):
    rec = vm.step()
    print(k, 'e=%.11f eR=%.11f eps=%.2e res=%.1e mv=%.0f gm=%.0e/%.0e S0=%.5f N=%.4f t=%.0fs' % (rec['e_site'], rec['eR'], rec['eps'], rec['res_max'], rec['mv_mean'], rec['gmL'], rec['gmR'], rec['S'][0], sum(vm.density()), rec['secs']), flush=True)
