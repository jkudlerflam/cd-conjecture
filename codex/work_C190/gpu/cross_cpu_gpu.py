# C5 cross-check between the CPU runs (z2dmrg.py) and the GPU runs (z2dmrg_xp.py) of the same configuration.
# usage: python3 -E -s cross_cpu_gpu.py <cpu folder> <gpu folder>. Only runs with status done in both are compared.
# Registered rule (RESULT.md C5 and the GPU amendment): where the energies agree within 1e-7, Sbar6 must agree within 0.01.
import json,sys,os,glob
C,G=sys.argv[1],sys.argv[2]
def fin(f):
    d=json.load(open(f))
    if d.get('status')!='done': return None
    h=d['hist']; a=h[-1]; chi=a['chi']; n=len([x for x in h if x['chi']==chi])
    return dict(chi=chi,n=n,tot=len(h),E=a.get('E_phys',a['E']),S1=a['Sbar1'],S6=a['Sbar6'],S12=a['Sbar12'],nm=a['n_mid'],res=a['res_max'],tr=a['trunc'])
for f in sorted(glob.glob(os.path.join(C,'*_L*_P*.json'))):
    t=os.path.basename(f)[:-5]; g=os.path.join(G,t+'.json')
    if not os.path.exists(g): continue
    a,b=fin(f),fin(g)
    if a is None or b is None: continue
    dE=a['E']-b['E']; dS=[a[k]-b[k] for k in ('S1','S6','S12')]
    v='ok' if abs(dS[1])<=0.01 else 'FAIL'
    if abs(dE)>1e-7: v='energies differ by more than 1e-7, lower one is used ('+('CPU' if dE<0 else 'GPU')+')'
    print(f"{t}: chi {a['chi']}/{b['chi']} sweeps {a['tot']}/{b['tot']} E_cpu={a['E']:.10f} E_gpu={b['E']:.10f} dE={dE:+.2e} "
          f"Sbar(1,6,12) cpu=({a['S1']:.4f},{a['S6']:.4f},{a['S12']:.4f}) gpu=({b['S1']:.4f},{b['S6']:.4f},{b['S12']:.4f}) dS6={dS[1]:+.1e} "
          f"n_mid {a['nm']:.4f}/{b['nm']:.4f} res {a['res']:.1e}/{b['res']:.1e} trunc {a['tr']:.1e}/{b['tr']:.1e} C5={v}")
