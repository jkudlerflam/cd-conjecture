import numpy as np,argparse,signal,os,json
from periodic_dmrg import PeriodicDMRG
from simple_dmrg import qutrit_model

def load_state(obj,prefix,N):
 z=np.load(prefix+'_mps.npz');obj.mps=[z[str(i)] for i in range(obj.L)]
 obj.q=[np.array([0],int)]
 for M in obj.mps:
  l,d,r=M.shape
  ix=np.argmax(abs(M.reshape(l*d,r)),axis=0)
  qr=obj.q[-1][ix//d]+ix%d
  invalid=(obj.q[-1][:,None,None]+np.arange(d)[None,:,None])!=qr[None,None,:]
  assert np.linalg.norm(M[invalid])<1e-8
  obj.q.append(qr)
 assert obj.q[-1][0]==N
 obj.rebuild_right()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--L',type=int);p.add_argument('--N',type=int);p.add_argument('--input');p.add_argument('--output');p.add_argument('--schedule',default='128,192,192,256,256');p.add_argument('--seconds',type=int,default=1500);args=p.parse_args()
 signal.alarm(args.seconds)
 with open(args.output+'_pid.txt','w') as f:f.write(str(os.getpid()))
 h,o=qutrit_model(.9,mu=0);obj=PeriodicDMRG(args.L,args.N,h,o)
 load_state(obj,args.input,args.N)
 obj.run([int(x) for x in args.schedule.split(',')],output=args.output,tol=1e-8)
