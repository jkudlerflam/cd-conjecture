"""Full-marginal (k=L) reference SDPs on tiny rings, with exact certificates."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from bootstrap import bonds,ring,certify


def run():
    folder=Path('ring_certificates'); folder.mkdir(exist_ok=True)
    results=[]
    for name,d,h in bonds():
        for L in ([4,6] if d==2 else [3,4]):
            H=ring(h,d,L); D=len(H)
            energies,U=eigh(H)
            lower,cert=certify(H,np.zeros((0,D,D)),np.zeros((0,D,D)),np.zeros((0,0,D,D)),
                               np.zeros(0),np.zeros(0),np.zeros((0,0)))
            v=np.array([int(a) for a in np.rint(U[:,0]*2**24)],dtype=object)
            H16=np.array([[int(a) for a in row] for row in H*16],dtype=object)
            un=int(v@H16@v); ud=16*int(v@v)
            cert.update({'model':name,'d':d,'L':L,'bond':h.tolist(),'trial_vector_integer':v.tolist(),
                         'upper_numerator':un,'upper_denominator':ud})
            (folder/(name+'_L'+str(L)+'.json')).write_text(json.dumps(cert,indent=2)+'\n')
            results.append({'model':name,'d':d,'L':L,'lower':lower,'upper':un/ud,'interval_width':un/ud-lower})
    Path('ring_intervals.json').write_text(json.dumps(results,indent=2)+'\n')


def verify():
    # Verify the operator and residual with Python integers.
    for path in sorted(Path('ring_certificates').glob('*.json')):
        c=json.loads(path.read_text()); d=c['d']; L=c['L']; D=d**L
        h16=[[int(round(a*16)) for a in row] for row in c['bond']]
        assert all(a*16==b for row,ir in zip(c['bond'],h16) for a,b in zip(row,ir))
        H=np.zeros((D,D),dtype=object)
        for col in range(D):
            word=[(col//(d**(L-1-i)))%d for i in range(L)]
            for i in range(L):
                j=(i+1)%L
                for a in range(d):
                    for b in range(d):
                        row=col+(a-word[i])*d**(L-1-i)+(b-word[j])*d**(L-1-j)
                        H[row,col]+=h16[a*d+b][word[i]*d+word[j]]
        K=np.array(c['K_integer'],dtype=object); Q=np.array(c['Q_integer'],dtype=object)
        assert np.array_equal(K,H*2**28)
        R=K*2**16-c['e_trial_integer']*2**16*np.eye(D,dtype=object)-Q@Q.T
        delta=max(sum(abs(a) for a in row) for row in R)
        assert delta==c['residual_row_bound_integer']
        assert c['bound_numerator']==c['e_trial_integer']*2**16-delta
        v=np.array(c['trial_vector_integer'],dtype=object)
        assert c['upper_numerator']==v@H@v
        assert c['upper_denominator']==16*(v@v)
        assert c['bound_numerator']*c['upper_denominator']<=c['upper_numerator']*c['bound_denominator']
        print('EXACT RING INTERVAL PASS',path)

if __name__=='__main__':
    run(); verify()
