"""Number conserving Gaussian LTI reduction, polynomial size; see PROOF.md."""
from core import *
def solve_cov(n,mu=.5):
    c=cp.Variable(n); C=cp.bmat([[c[abs(i-j)] for j in range(n)] for i in range(n)])
    p=cp.Problem(cp.Minimize(c[1]+mu*(c[0]-.5)),[C>>0,np.eye(n)-C>>0])
    p.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10)
    return p.value,c.value
if __name__=='__main__':
    from itertools import permutations,combinations
    def det(a):
        n=len(a);s=0
        for p in permutations(range(n)):
            sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n));v=sign
            for i in range(n):v*=int(a[i,p[i]])
            s+=v
        return s
    out={}
    for n in [3,4]:
        e,c=solve_cov(n); print(n,e,c,flush=True)
        out[str(n)]=dict(e=e,c=c.tolist())
    c=np.array(out['4']['c']);c=.998*c;c[0]+=.001
    ints=np.round(c*10000).astype(int);C=np.array([[ints[abs(i-j)] for j in range(4)] for i in range(4)])
    minors=[]
    for A in [C,10000*np.eye(4,dtype=int)-C]:
        minors.append([det(A[np.ix_(ii,ii)]) for k in range(1,5) for ii in combinations(range(4),k)])
    out['rational']=dict(denominator=10000,coefficients=ints.tolist(),principal_minors=minors,energy=float((ints[1]+.5*ints[0])/10000-.25))
    assert all(x>0 for row in minors for x in row)
    print(out['rational']);write_json('free_covariance.json',out)
