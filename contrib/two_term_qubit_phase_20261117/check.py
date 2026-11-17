"""Independent matrix and exact integer Pauli-algebra checks."""
import itertools,json,time
from pathlib import Path
import numpy as np
from model import Block,physical_full,pauli_full,compass_phases


def multiply(a,b):
    phase=1;out=[]
    for x,y in zip(a,b):
        if x=='I':out.append(y)
        elif y=='I':out.append(x)
        elif x==y:out.append('I')
        else:
            phase*=1j if (x,y) in [('X','Y'),('Y','Z'),('Z','X')] else -1j
            out.append(({'X','Y','Z'}-{x,y}).pop())
    return ''.join(out),phase


def comm(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            z,p=multiply(x,y);_,q=multiply(y,x)
            out[z]=out.get(z,0)+c*d*(p-q)
    return {x:c for x,c in out.items() if c}


def add(a,b):
    out=a.copy()
    for x,c in b.items():out[x]=out.get(x,0)+c
    return {x:c for x,c in out.items() if c}


def run():
    start=time.perf_counter();errors=[]
    # Four exact integer samples certify this polynomial identity of degree <=3.
    for g in (0,1,2,3):
        a={'XYI':2,'ZZI':g};b={'IXY':2,'IZZ':g}
        actual=comm(add(a,b),comm(a,b))
        want={'IXY':8*g*g+32,'XYI':-8*g*g-32,'IZZ':16*g,'ZZI':-16*g,'YIY':32*g,'XIX':-32*g}
        assert actual=={x:c for x,c in want.items() if c}
        # Arbitrary range <=2 correction cannot produce a separated XIX/YIY.
        for word in itertools.product('IXYZ',repeat=2):
            c=''.join(word)
            ans=comm(add(a,b),{c+'I':1,'I'+c:1})
            assert ans.get('XIX',0)==ans.get('YIY',0)==0
    center='IXYI';leaves=['XYII','IIXY','ZZII']
    assert all(comm({center:1},{x:1}) for x in leaves)
    assert all(not comm({x:1},{y:1}) for x,y in itertools.combinations(leaves,2))
    for L in (4,6,8):
        H=physical_full(L,.5).toarray();u=compass_phases(L);R=u[:,None]*H*u.conj()[None,:];vals=[]
        blocks={}
        for p,q in ((1,1),(1,-1),(-1,1),(-1,-1)):
            b=Block(L,p,q);A,D=b.operators();K=A+.5*D;blocks[p,q]=K
            U=np.array([b.expand(v) for v in np.eye(len(b.states))]).T
            errors.append(float(np.max(np.abs(U.conj().T@R@U-K.toarray()))))
            if L<=6:
                for j in range(L):
                    for letters,op in [('xy',b.A[j]),('yx',b.B[j]),('zz',b.D[j])]:
                        P=pauli_full(L,{j:letters[0],(j+1)%L:letters[1]}).toarray()
                        P=u[:,None]*P*u.conj()[None,:]
                        errors.append(float(np.max(np.abs(U.conj().T@P@U-op.toarray()))))
                other=Block(L,p,-q)
                V=np.array([other.expand(v) for v in np.eye(len(other.states))]).T
                Z=pauli_full(L,{0:'z'}).toarray()
                errors.append(float(np.max(np.abs(V.conj().T@Z@U-np.diag(b.z[:,0])))))
            T=b.translation(1 if p==1 else 2)
            errors.append(float(np.max(np.abs((T@K-K@T).toarray()))))
            vals.extend(np.linalg.eigvalsh(K.toarray()))
        errors.append(float(np.max(np.abs(np.sort(vals)-np.linalg.eigvalsh(H)))))
        b=Block(L,-1,1);T=b.translation()
        errors.append(float(np.max(np.abs((T@blocks[-1,1]-blocks[-1,-1]@T).toarray()))))
        # Physical four-site sign duality: I,Y,Z,X, repeated.
        W=pauli_full(L,{j:['i','y','z','x'][j%4] for j in range(L) if j%4})
        for periodic in ((False,True) if L%4==0 else (False,)):
            e=W@physical_full(L,.5,periodic)@W-physical_full(L,-.5,periodic)
            errors.append(float(np.max(np.abs(e.toarray()))))
        # The full periodic boost-charge commutator has the advertised obstruction.
        h=[2*pauli_full(L,{j:'x',(j+1)%L:'y'})+pauli_full(L,{j:'z',(j+1)%L:'z'}) for j in range(L)]
        H1=sum(h);Q=sum(h[j]@h[(j+1)%L]-h[(j+1)%L]@h[j] for j in range(L))
        target=32*sum(pauli_full(L,{j:'y',(j+2)%L:'y'})-pauli_full(L,{j:'x',(j+2)%L:'x'}) for j in range(L))
        errors.append(float(np.max(np.abs((H1@Q-Q@H1-target).toarray()))))
    assert max(errors)<1e-10
    return dict(status='passed',max_matrix_error=max(errors),lengths=[4,6,8],
        exact_polynomial_samples=[0,1,2,3],seconds=time.perf_counter()-start,
        checks=['literal physical spectrum equals all four real blocks','individual physical XY/YX/ZZ and parity-changing density observables','physical translation and paired odd-parity blocks','four-site g-sign duality','claw obstruction','local double commutator','no range-two repair','full periodic boost-charge obstruction'])


if __name__=='__main__':
    result=run();(Path(__file__).parent/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
