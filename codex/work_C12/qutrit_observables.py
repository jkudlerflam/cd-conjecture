import numpy as np
from scipy.linalg import svdvals

def occupations(states,L):
    return np.stack([((states>>np.uint64(2*j))&np.uint64(3)).astype(np.int8) for j in range(L)],axis=1)

def entropy_profile(states,v,L,cuts=None):
    occ=occupations(states,L);cs=np.cumsum(occ,axis=1)
    if cuts is None:cuts=range(1,L//2+1)
    out=[]
    for l in cuts:
        lc=states&np.uint64((1<<(2*l))-1);rc=states>>np.uint64(2*l)
        probs=[]
        for q in np.unique(cs[:,l-1]):
            ix=np.flatnonzero(cs[:,l-1]==q)
            aa,ai=np.unique(lc[ix],return_inverse=True);bb,bi=np.unique(rc[ix],return_inverse=True)
            mat=np.zeros((len(aa),len(bb)),dtype=v.dtype);mat[ai,bi]=v[ix]
            probs.extend(svdvals(mat,check_finite=False)**2)
        p=np.array(probs);p=p[p>1e-14]
        out.append(dict(cut=int(l),S=float(-sum(p*np.log(p))),S2=float(-np.log(sum(p*p))),norm=float(sum(p))))
    return out

def cfit(profile,L,exclude=2):
    dat=[r for r in profile if r['cut']>=exclude]
    x=np.log(L/np.pi*np.sin(np.pi*np.array([r['cut'] for r in dat])/L))
    c1=3*np.polyfit(x,[r['S'] for r in dat],1)[0]
    c2=4*np.polyfit(x,[r['S2'] for r in dat],1)[0]
    return dict(c_vn=float(c1),c_renyi2=float(c2),exclude=exclude)

def describe(states,es,vs,L):
    shifted=((states<<np.uint64(2))&np.uint64((1<<(2*L))-1))|(states>>np.uint64(2*(L-1)))
    dest=np.searchsorted(states,shifted)
    Tsub=vs[dest,:].T.conj()@vs
    occ=occupations(states,L)
    out=[]
    first=0
    while first<len(es):
        last=first+1
        while last<len(es) and abs(es[last]-es[first])<1e-7:last+=1
        ts=Tsub[first:last,first:last]
        if last-first>1:
            te,tv=np.linalg.eig(ts.astype(complex))
            mom=np.angle(te).tolist();mod=np.abs(te).tolist()
        else:
            # A real standing wave can have |<T>|<1 if its partner was not returned.
            te=ts[0,0]
            mom=[float(np.arccos(np.clip(np.real(te),-1,1)))];mod=[float(abs(te))]
        doublons=[float(np.sum(abs(vs[:,a])**2*np.sum(occ==2,axis=1))) for a in range(first,last)]
        out.append(dict(E=float(es[first]),multiplicity=last-first,momenta=mom,translation_moduli=mod,doublons=doublons))
        first=last
    return out

def momentum_eigenstate(states,v,L):
    shifted=((states<<np.uint64(2))&np.uint64((1<<(2*L))-1))|(states>>np.uint64(2*(L-1)))
    dest=np.searchsorted(states,shifted)
    tv=np.empty_like(v);tv[dest]=v
    cosine=float(np.real(np.vdot(v,tv)))
    if abs(cosine)>1-1e-8:
        return v,dict(cosine=cosine,momentum=float(np.arccos(np.clip(cosine,-1,1))),translation_residual=float(np.linalg.norm(tv-np.sign(cosine)*v)))
    sine=np.sqrt(1-cosine*cosine)
    w=(tv-cosine*v)/sine
    vc=(v+1j*w)/np.sqrt(2)
    vc/=np.linalg.norm(vc)
    tvc=np.empty_like(vc);tvc[dest]=vc
    phase=np.vdot(vc,tvc)
    return vc,dict(cosine=cosine,momentum=float(np.angle(phase)),translation_residual=float(np.linalg.norm(tvc-phase*vc)))
