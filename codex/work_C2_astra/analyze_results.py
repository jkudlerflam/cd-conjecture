import json, numpy as np
r=json.load(open('results_extended.json'))
labels=['TFIM','Heisenberg','OF0.2','OF0.35','J1J2']
def fit(xs,ys,powers):
 X=np.array([[n**(-p) for p in powers] for n in xs]); b=np.linalg.lstsq(X,ys,rcond=None)[0]
 err=np.sqrt(np.mean((X@b-ys)**2));return b,err
for label in labels:
 sub=[x for x in r if x['label']==label]
 print('\n'+label)
 print('| N | c2 | m0 N² | m1 N³ | Var(J)/N | cMV |')
 print('|---:|---:|---:|---:|---:|---:|')
 for x in sub:
  a=x['centered'];print(f"|{x['N']}|{a['c2']:.6f}|{a['m0N2']:.6f}|{a['m1N3']:.6f}|{a['JvarN']:.8f}|{a['cMV']:.6f}|")
 for lo in [12,14,16]:
  subfit=[x for x in sub if x['N']>=lo]; ns=[x['N'] for x in subfit]
  if len(ns)<3:continue
  for key in ['c2','cMV']:
   y=[x['centered'][key] for x in subfit]
   print('fit',lo,ns[-1],key,'c+a/N²+b/N⁴',fit(ns,y,[0,2,4]))
  if 'Aminus2' in subfit[0]['centered']:
   co={}
   for key in ['Aminus2','Aminus1','JvarN']:
    co[key]=fit(ns,[x['centered'][key] for x in subfit],[0,2,4])[0][0]
   co['a']=16*np.pi**2*co['Aminus2'];co['d']=16*np.pi**2*co['Aminus1']
   co['K']=128*np.pi**2*co['Aminus2']**3/co['Aminus1']**2
   print('UV extrapolation',co)
print('max ground residual',max(x['ground_residual'] for x in r))
print('max translation residual',max(x['translation_residual'] for x in r))
print('min uniqueness gap',min(x['gapcheck'] for x in r))
