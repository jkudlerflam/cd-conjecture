from ks_test import *
rows=[]
for label,model,c in [('TFIM','ising',0),('Heisenberg','spin',0)]:
 for N in [8,10,12,14,16,18,20,22,24]:
  sec='even' if model=='ising' else 'sz0';k=0 if model=='ising' or N%4==0 else N//2
  o=Orbits(N,sec);d=density(N,model,c); ht=total(d,N)
  H=o.matrix(ht,k);e,g=lowest(H,1);Hp=o.matrix(ht,(k+1)%N)
  hv=o.matrix(total(d,N,1),k,(k+1)%N)@g[:,0]
  m0=np.vdot(hv,hv).real;m1=np.vdot(hv,Hp@hv-e[0]*hv).real
  row=dict(label=label,N=N,m0N3=float(m0*N**3),m1N3=float(m1*N**3));rows.append(row)
  print(row,flush=True)
  json.dump(rows,open('mode_one_results.json','w'),indent=2)
