from periodic_dmrg import PeriodicDMRG
from simple_dmrg import qutrit_model
from qutrit_ed import solve
L,N=8,3;h,o=qutrit_model(.9,mu=0)
e,v,s,r=solve(L,N,k=1,U=-14/9+.9)
obj=PeriodicDMRG(L,N,h,o);res=obj.run([8,16]+[27]*10,output='periodic_validation',tol=1e-11)
print('ED / periodic DMRG',e[0],res[-1]['E'],'error',res[-1]['E']-e[0])
assert abs(res[-1]['E']-e[0])<1e-8
