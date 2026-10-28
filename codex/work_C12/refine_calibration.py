import numpy as np
from simple_dmrg import DMRG
swap=np.zeros((9,9))
for a in range(3):
 for b in range(3):swap[3*a+b,3*b+a]=1
for L in [24,48]:
 obj=DMRG(L,swap,np.zeros((3,3)))
 z=np.load('calibration_su3_L'+str(L)+'_mps.npz')
 obj.mps=[z[str(i)] for i in range(L)];obj.right_canon();obj.rebuild_right()
 obj.run([128,128,192,192],output='calibration_su3_refined_L'+str(L),tol=1e-10)
