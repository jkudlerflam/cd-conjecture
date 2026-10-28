import json
from momentum_ed import OrbitHamiltonian
results=[]
for V12,V22 in [(20,-2),(20,0),(20,2),(20,5),(20,20),(20,100),(5,0),(5,5),(5,20),(0,0),(0,20)]:
 row=dict(V12=V12,V22=V22,levels=[])
 for L in [18,24,30]:
  ob=OrbitHamiltonian(L,4,V12=V12,V22=V22)
  energies={}
  for k in [0,L//3]:
   es,vs,r=ob.solve(k,nev=1,tol=1e-11);energies[str(k)]=float(es[0])
  row['levels'].append(dict(L=L,opposite=energies['0'],same=energies[str(L//3)],difference=energies['0']-energies[str(L//3)]))
 print(json.dumps(row),flush=True);results.append(row)
 with open('scattering_screen_results.json','w') as f:json.dump(results,f,indent=2)
