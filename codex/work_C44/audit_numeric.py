"""C44 wrapped C39 source audit. NumPy/SciPy only, bounded small instances.
Uses attributed local snapshots, core_literal modified to expose each bond.
All charge and reachable helix levels enter the source contraction.
"""
import json,time,itertools
from pathlib import Path
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import core_literal,core_model
import wrapper_tools as wt

def norm(A): return float(la.norm(A,2))
def tl(A): return A-np.trace(A)/len(A)*np.eye(len(A))
def core(n):
 c,_,_,_=core_literal.make(n,a=4.,g=.4)
 c['H']=c['H'].toarray()
 return c

def buscore(c,bus):
 d=dict(c); d['ex']=c['typ']==bus
 d['F']={side:{f:c['F'][bus,side,i].toarray() for i,f in enumerate(('0','+','-'))} for side in ('L','R')}
 return d

def contract(me1,me2,h1,h2,phase):
 e,aa,_=me1;f,bb,_=me2;A=aa['R'];B=bb['L']
 C=np.zeros((len(e),len(f)))
 for h,a in zip(*h1):
  for hh,b in zip(*h2): C+=a*b/(e[:,None]+f[None,:]+h+hh)
 terms=[('0','0',1.),('+','-',np.exp(1j*phase)),('-','+',np.exp(-1j*phase))]
 N=A['0'].shape[1];H=np.zeros((N*N,N*N),complex)
 for ar,al,v in terms:
  for br,bl,w in terms:
   X=np.einsum('ei,ej->eij',A[ar].conj(),A[br]).reshape(len(e),-1)
   Y=np.einsum('ei,ej->eij',B[al].conj(),B[bl]).reshape(len(f),-1)
   H-=v.conjugate()*w*(X.T@C@Y).reshape(N,N,N,N).transpose(0,2,1,3).reshape(N*N,N*N)
 return H

def source_scan():
 start=time.monotonic();out=[];saved={}
 hzero=(np.array([0.]),np.array([1.]))
 for n in (1,2):
  c=core(n);N=2**n;eta=.6
  K=core_model.exact_bond(c,eta);saved[f'bare_n{n}']=K
  first=max(norm(c['G'].T@F@c['G']) for F in c['F'].values())
  for wvac in (.16,.10):
   cs={p:wt.charge(c['k'] if 'k' in c else n+2,p,wvac) for p in (1,2)}
   b=n+3;phi=1/b;comp=wt.compensator(cs[1]['r'],cs[2]['r'],b,phi)
   beta=wt.transport(cs[1],cs[2]);ab=sum(np.sqrt(2/b)*np.sin(np.pi*j/b)*np.sqrt(2/(b+1))*np.sin(np.pi*j/(b+1)) for j in range(1,b))
   flat={p:dict(c,H=c['ends']+np.einsum('j,jab->ab',cs[p]['avgD'],c['parts'])) for p in cs}
   ref={(p,q):np.zeros_like(K) for p in cs for q in cs}
   for bus,phase,weight in [(1,np.pi/2,1),(2,np.pi,eta)]:
    mm={p:wt.measure(buscore(flat[p],bus)) for p in cs}
    for p,q in ref:ref[p,q]+=weight*contract(mm[p],mm[q],hzero,hzero,phase)
   for ratio in (.03,.003):
    s=ratio*min(ch['gap'] for ch in cs.values());ks={pq:np.zeros_like(K) for pq in ref};res=0.
    hm={p:wt.helix(ch['N'],comp,phi,s) for p,ch in cs.items()}
    for bus,phase,weight in [(1,np.pi/2,1),(2,np.pi,eta)]:
     mm={p:wt.measure(buscore(c,bus),ch,s) for p,ch in cs.items()}
     res=max(res,*[mm[p][2] for p in cs])
     for p,q in ks:ks[p,q]+=weight*contract(mm[p],mm[q],hm[p],hm[q],phase)
    # All K values are s*K2, not lambda^2*K2.
    key=f'n{n}_w{wvac}_t{ratio}'
    for (p,q),v in ks.items():saved[f'{key}_{p}{q}']=v
    Dint=ks[2,2]-ks[2,1]-ks[1,2]+ks[1,1]
    cc=float(np.trace(Dint).real/len(K))
    rec=dict(n=n,wvac=wvac,ratio=ratio,s=s,core_dimension=c['d'],
       charge_dimensions=[ch['d'] for ch in cs.values()],charge_gaps=[ch['gap'] for ch in cs.values()],
       source_dimensions_per_bus=[ch['d']*3*(2*n+2)*N*(ch['N']+1) for ch in cs.values()],
       first_order_port_norm=first,code_residual=norm(c['H']@c['G']),
       adjacency=[ch['avgD'].tolist() for ch in cs.values()],endpoint_weights=[ch['r'] for ch in cs.values()],
       beta_length=float(beta*ab**2*np.cos(phi/2)**(2*b)),compensator_norm=norm(comp),
       graph_residual=res,bare_traceless_norm=norm(tl(K)),
       bare_error=max(norm(v-K) for v in ks.values()),
       length_error=max(norm(v-ks[1,1]) for v in ks.values()),
       fast_error=max(norm(ks[pq]-ref[pq]) for pq in ks),
       limiting_bare_error=max(norm(v-K) for v in ref.values()),
       limiting_length_error=max(norm(v-ref[1,1]) for v in ref.values()),
       scalar_density_interaction=cc,mixed_density_payload_norm=norm(tl(Dint)),
       max_hermiticity=max(norm(v-v.conj().T) for v in ks.values()))
    assert rec['max_hermiticity']<1e-12 and res<1e-8
    out.append(rec);print(key,rec['length_error'],rec['fast_error'],cc,flush=True)
    Path('source_results.json').write_text(json.dumps(dict(status='numerical',cases=out,seconds=time.monotonic()-start),indent=2)+'\n')
    np.savez('source_bonds.npz',**saved)


def filter_parts(k):
 # Untwisted and rotated filter. Both boundary pins are |1><1|.
 sing=np.array([0,1,-1,0])/np.sqrt(2);h=np.outer(sing,sing);pin=np.diag([0,1])
 parts=[sp.kron(sp.eye(2**j),sp.kron(sp.csr_matrix(h),sp.eye(2**(k-j-2))),format='csr') for j in range(k-1)]
 ends=sp.kron(pin,sp.eye(2**(k-1)))+sp.kron(sp.eye(2**(k-1)),pin)
 return parts,ends.tocsr()

def cell_ed():
 start=time.monotonic();out=[]
 n=1;k=n+2;c=core(n);fp,fe=filter_parts(k);df=2**k;dc=c['d'];sf=.2
 kap=1/((k+2)*(2*k+3));bet=np.sin(np.pi/(k+2))**2
 eps=kap*bet/((1+np.sqrt(kap))**2+bet);jf=2/eps
 for p in (1,2):
  ch=wt.charge(k,p,.20);dch=ch['d'];raw_e=float(ch['chi']@ch['H']@ch['chi'])
  H=sp.kron(sp.csr_matrix(ch['H']-raw_e*np.eye(dch)),sp.eye(dc*df),format='csr')
  ends=sf*sp.kron(c['ends'],sp.eye(df))+jf*sp.kron(sp.eye(dc),fe)
  H+=sp.kron(sp.eye(dch),ends)
  for j in range(k-1):
   B=sf*sp.kron(c['parts'][j],sp.eye(df))+jf*sp.kron(sp.eye(dc),fp[j])
   H+=sp.kron(sp.diags(ch['D'][:,j]),B)
  f0=np.zeros((df,1));f0[0]=1
  G=np.kron(ch['chi'][:,None],np.kron(c['G'],f0))
  # Shift inverse resolves the small charge gap without slow Krylov convergence.
  ev,V=sla.eigsh(H,k=4,sigma=-1e-7,which='LM',tol=1e-10)
  order=np.argsort(ev);ev=ev[order];V=V[:,order]
  residual=max(la.norm(H@V[:,j]-ev[j]*V[:,j]) for j in range(4))
  Vcode=la.qr(V[:,:2],mode='economic')[0]
  out.append(dict(n=n,p=p,dimension=H.shape[0],filter_scale=jf,core_scale=sf,
   eigenvalues=ev.tolist(),eigenpair_residual=residual,code_residual=norm(H@G),
   code_subspace_error=norm(G-Vcode@(Vcode.conj().T@G)),charge_gap=ch['gap'],
   scope='All legal initialized core states, every vacancy configuration, and all occupation-filter states. Endpoint helix tracks factor exactly and are checked separately.'))
  print('cell',out[-1],flush=True)
  assert out[-1]['code_residual']<1e-10 and out[-1]['code_subspace_error']<1e-7
 Path('cell_ed_results.json').write_text(json.dumps(dict(cases=out,seconds=time.monotonic()-start),indent=2)+'\n')

if __name__=='__main__':
 import sys
 if len(sys.argv)>1 and sys.argv[1]=='cell':cell_ed()
 else:source_scan()
