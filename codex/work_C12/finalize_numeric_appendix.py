import json,numpy as np,platform,scipy
from qutrit_observables import cfit
lines=['## Final numerical appendix','']
base=json.load(open('ed_L24_results.json'))
extra=json.load(open('ed_L24_remaining_sectors.json'))
levels=base['levels']+extra
checked=sorted(set(r['k'] for r in levels));missing=sorted(set(range(13))-set(checked))
best=min(levels,key=lambda r:r['energies'][0])
if not missing:
 lines.append('At L=24, all momentum sectors k=0,...,12 were checked. Reality of h gives the remaining sectors by complex conjugation. The lowest state is in k=%d, with E=%.11f. The largest reported eigenvector residual in this check is %.2g.'%(best['k'],best['energies'][0],max(r['residual'] for r in levels)))
else:
 lines.append('At L=24, momentum indices %s were checked; %s remain unchecked. The lowest state found is in k=%d. The quoted L=24 entropies therefore refer to that state, without a complete global-sector guarantee.'%(checked,missing,best['k']))
lines+=['','The completed periodic DMRG checkpoints are listed below. These are raw finite-bond-dimension fits, not claims for the thermodynamic central charge.','',
'| L | Largest completed chi | Energy | Half-chain S | Half-chain S2 | Shape fit, vN | Shape fit, Renyi2 | Density range |',
'|---:|---:|---:|---:|---:|---:|---:|---|']
for L in [24,36,48]:
 r=json.load(open('periodic_L'+str(L)+'_highchi.json'))[-1]
 p=[x for x in r['entropy'] if x['cut']<=L//2];fit=cfit(p,L);ns=np.array(r['density'])[:,0]
 lines.append('| %d | %d | %.11f | %.6f | %.6f | %.4f | %.4f | %.4f to %.4f |'%(L,r['chi'],r['E'],p[-1]['S'],p[-1]['S2'],fit['c_vn'],fit['c_renyi2'],ns.min(),ns.max()))
 if L==24:
  error=r['E']-base['levels'][0]['energies'][0]
  Serror=base['entropy'][-1]['S']-p[-1]['S']
  comparison=(error,Serror)
lines+=['','At L=24, the last checkpoint has energy error %.3g and half-chain entropy deficit %.3g relative to ED. Density modulation and continuing bond-dimension dependence prevent a reliable large-L extrapolation. The interrupted runs and wall-clock stopping condition are recorded in PROGRESS.md.'%comparison]
try:
 rr=json.load(open('ed_delta075_results.json'))
 lines+=['','The nearby detuning delta=0.75 is also affected by momentum and size effects.']
 for r in rr:
  lines.append('At L=%d, the lowest tested momentum is k=%d, the shape estimates are %.4f and %.4f, and half-chain S=%.6f.'%(r['L'],r['ground_k'],r['fit']['c_vn'],r['fit']['c_renyi2'],r['half']['S']))
 lines.append('All sectors were checked at L=21; only k=0,1,4,8,12 were checked at L=24 in this auxiliary detuning scan. Its molecular minima have already disappeared, so this does not rescue the original three-fluid argument.')
except FileNotFoundError:pass
lines+=['','Numerical environment: Python %s, NumPy %s, SciPy %s. BLAS was limited to one thread per Python process. All random scans use fixed seeds in their scripts.'%(platform.python_version(),np.__version__,scipy.__version__),'']
with open('NUMERICAL_APPENDIX.md','w') as f:f.write('\n'.join(lines))
print('\n'.join(lines))
