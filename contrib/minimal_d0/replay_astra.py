"""Independent replay of Astra's explicit positive-weight recipes (astra_frames/explicit_recipes_m1_m2.json).
Own gate matrices, own Bell state and own target; only the gate lists are taken from Astra's export."""
import json, numpy as np
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)
H=np.array([[1,1],[1,-1]])/np.sqrt(2); S=np.diag([1,1j]); Sd=S.conj().T
one={'H':H,'S':S,'Sd':Sd,'X':X,'Y':Y,'Z':Z}
P={0:I2,1:X,2:Y,3:Z}
def op1(U,j,n):
    o=np.eye(1)
    for k in range(n): o=np.kron(o,U if k==j else I2)
    return o
def cx(c,t,n):
    P0=np.diag([1,0]);P1=np.diag([0,1])
    return op1(P0,c,n)+op1(P1,c,n)@op1(X,t,n)
def swap(a,b,n): return cx(a,b,n)@cx(b,a,n)@cx(a,b,n)
def gate(g,n):
    if g[0]=='CX': return cx(g[1],g[2],n)
    if g[0]=='SWAP': return swap(g[1],g[2],n)
    return op1(one[g[0]],g[1],n)
def kron(*a):
    o=np.eye(1)
    for x in a:o=np.kron(o,x)
    return o
data=json.load(open('astra_frames/explicit_recipes_m1_m2.json'))
worst=0
for key,rec in data.items():
    m=int(key.split(',')[0][2:]); q=2**m; n=m+1
    lab=key.split('Q=')[1]; Ql=[] if lab=='empty' else ['IXYZ'.index(ch) for ch in lab]
    Om=np.eye(q).reshape(-1)/np.sqrt(q); Pi=np.kron(np.outer(Om,Om),I2)
    Dm=2*q*q; tot=np.zeros((Dm,Dm),complex)
    entries=rec if isinstance(rec,list) else rec.get('recipe',rec.get('frames'))
    for e in entries:
        W=np.eye(2**n,dtype=complex)
        for g in e['gates']: W=gate(g,n)@W        # chronological: later gates on the left
        Wf=np.kron(np.eye(q),W)
        tot+=e['weight']*Wf.conj().T@Pi@Wf
    T=np.zeros((Dm,Dm),complex); Lq=np.zeros((Dm,Dm),complex)
    Om_=np.outer(Om,Om)
    for a in (1,2):
        F=kron(P[a],*[P[x] for x in Ql]); Fp=np.kron(F,np.eye(q))
        T+=np.kron(Fp@Om_+Om_@Fp,P[a])
        Lq+=2/q**2*kron(np.eye(q),F.T,P[a])
    res=tot-(2-4/q**2)*np.eye(Dm)-(T-Lq)
    worst=max(worst,np.abs(res).max())
    print(key,len(entries),'frames, max residual',f'{np.abs(res).max():.2e}')
print('worst',worst)
