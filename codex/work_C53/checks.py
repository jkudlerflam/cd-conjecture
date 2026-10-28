"""Small exact diagonalization for C53. No infrared inference from finite sizes."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy.integrate import quad

I = np.eye(2)
Z = np.diag([1., -1.])
SP = np.array([[0., 1.], [0., 0.]])
SM = SP.T

def kron(*xs):
    r = np.ones((1, 1))
    for x in xs: r = np.kron(r, x)
    return r

def at(o, j, n):
    return kron(*[o if k == j else I for k in range(n)])

def perm(o, positions):
    n = len(positions)
    inv = list(np.argsort(positions))
    return o.reshape([2]*(2*n)).transpose(inv+[n+i for i in inv]).reshape(2**n, 2**n)

def ingredients():
    e = np.eye(4)
    g, a, u, d = (e[0]+e[3])/np.sqrt(2), (e[0]-e[3])/np.sqrt(2), e[1], e[2]
    hp = kron(np.eye(4)-np.outer(g,g), I)
    r = .6*(kron(np.outer(a,u),SP)-kron(np.outer(a,d),SM))
    hp += r+r.T
    Ghalf = kron(g[:,None], I)
    A = kron(hp,np.eye(8))+kron(np.eye(8),hp)
    V = 2*(kron(at(SP,0,3),at(SM,0,3))+kron(at(SM,0,3),at(SP,0,3)))-kron(at(Z,0,3),at(Z,0,3))
    G = kron(Ghalf,Ghalf)
    cell = kron(hp,np.eye(4))+perm(kron(hp,np.eye(4)),[4,3,2,0,1])
    Gcell = np.column_stack([np.kron(np.kron(g,I[:,c]),g) for c in (0,1)])
    return A,V,G,cell,Gcell

def main():
    A,V,G,cell,Gcell = ingredients()
    alpha = F(9,16); t = F(-9,32); u = F(243,11648); a0 = F(-20505,11648)
    hop = kron(SP,SM)+kron(SM,SP)
    zz = kron(Z,Z)
    zsum = kron(Z,I)+kron(I,Z)
    ev,U = eigh(A)
    keep = ev > 1e-10
    R = (U[:,keep]/ev[keep])@U[:,keep].T
    h2 = -G.T@V@R@V@G
    target2 = float(a0)*np.eye(4)+float(t)*hop+float(u)*zz
    err = np.max(np.abs(h2-target2))
    assert err < 1e-12
    rows=[]
    for coupling in (F(0),F(1,50),F(-1,50)):
        q = alpha*coupling/4-u
        assert 4*(u+q)/alpha == coupling
        field = F(1,3)
        W = float(q)*kron(at(Z,2,3),at(Z,2,3))-float(alpha*field/2)*(at(Z,2,6)+at(Z,5,6))
        target = target2 + float(q)*zz-float(alpha*field/2)*zsum
        assert np.max(np.abs(h2+G.T@W@G-target)) < 1e-12
        residuals=[]
        for lam in (.02,.01,.005,.0025):
            H = A+lam*V+lam**2*W
            w,vec=eigh(H,subset_by_index=[0,3])
            overlap=G.T@vec
            left,sv,right=np.linalg.svd(overlap)
            rot=left@right
            Heff=rot@np.diag(w)@rot.T
            rem=Heff-lam**2*target
            residuals.append({'lambda':lam,'error':float(np.linalg.norm(rem,2)),
                              'error_over_lambda3':float(np.linalg.norm(rem,2)/lam**3),
                              'min_overlap':float(sv.min())})
        assert all(residuals[i+1]['error'] < residuals[i]['error']/7 for i in range(3))
        rows.append({'g':str(coupling),'diagonal_counterterm':str(q),'residuals':residuals})

    # Complete 32^2 physical cell space, not just the 64-dimensional excursion.
    lam=.1; scale=float(alpha)*lam**2; coupling=.02; field=1/3
    A2=kron(cell,np.eye(32))+kron(np.eye(32),cell)
    cplus=at(SP,2,5); cz=at(Z,2,5)
    logical=-.5*hop+coupling/4*zz-field/2*zsum
    Lfull=-.5*(kron(cplus,cplus.T)+kron(cplus.T,cplus))+coupling/4*kron(cz,cz)-field/2*(kron(cz,np.eye(32))+kron(np.eye(32),cz))
    H=A2+scale*Lfull
    GG=kron(Gcell,Gcell)
    intertwiner=float(np.max(np.abs(H@GG-GG@(scale*logical))))
    ew=eigh(H,eigvals_only=True)
    lw=eigh(scale*logical,eigvals_only=True)
    specerr=float(np.max(np.abs(ew[:4]-lw)))
    threshold=float(ew[4]-lw[0])
    assert intertwiner<1e-12 and specerr<1e-12 and threshold>=.4-1e-12
    thermal=[]
    for beta in (5.,10.,20.):
        # Common shift avoids underflow, cancels from the ratio.
        ratio=float(np.exp(-beta*(ew-ew[0])).sum()/np.exp(-beta*(lw-ew[0])).sum())
        bound=(1+15*np.exp(-.4*beta))**2
        assert 1-1e-12 <= ratio <= bound+1e-12
        thermal.append({'beta':beta,'Z_full_over_Z_logical':ratio,'bound':float(bound)})

    # Fixed-number logical spin ED. Physical periodic spin boundary, with the
    # parity-dependent fermionic twist included in the free comparison.
    sectors=[]
    for M in (8,10,12):
        for coupling in (0.,.02):
            energies=[]
            free_errors=[]
            for N in range(M+1):
                states=[sum(1<<j for j in occ) for occ in combinations(range(M),N)]
                index={s:j for j,s in enumerate(states)}
                rr=[]; cc=[]; vv=[]
                for i,s in enumerate(states):
                    diag=-field*(N-M/2)
                    for j in range(M):
                        k=(j+1)%M
                        nj=(s>>j)&1; nk=(s>>k)&1
                        diag+=coupling*(nj-.5)*(nk-.5)
                        if nj!=nk:
                            rr.append(i); cc.append(index[s^(1<<j)^(1<<k)]); vv.append(-.5)
                    rr.append(i); cc.append(i); vv.append(diag)
                mat=csr_matrix((vv,(rr,cc)),shape=(len(states),len(states)))
                if len(states)<=2: E=float(eigh(mat.toarray(),eigvals_only=True)[0])
                else:
                    vals,vec=eigsh(mat,k=1,which='SA',tol=1e-11,v0=np.random.default_rng(53+M+N).normal(size=len(states)))
                    E=float(vals[0]); assert np.linalg.norm(mat@vec[:,0]-E*vec[:,0])<1e-9
                energies.append(E)
                if coupling==0:
                    shift=0. if N%2 else .5
                    ks=2*np.pi*(np.arange(M)+shift)/M
                    Efree=float(np.sort(-np.cos(ks)-field)[:N].sum()+field*M/2)
                    free_errors.append(abs(E-Efree))
            n=int(np.argmin(energies))
            sectors.append({'M':M,'g':coupling,'ground_N':n,'density':n/M,
                            'ground_energy':energies[n],
                            'max_free_spectrum_error':max(free_errors) if free_errors else None})
            if free_errors: assert max(free_errors)<1e-10
            assert n>M/2

    # Free thermal coefficient, dimensionless logical energy units.
    kf=np.arccos(-field); vf=np.sin(kf)
    e0=(-np.sin(kf)-field*kf)/np.pi+field/2
    thermal_free=[]
    for T in (.04,.02,.01,.005):
        f=field/2-quad(lambda k:T*np.logaddexp(0,(np.cos(k)+field)/T),0,np.pi,
                      epsabs=1e-12,points=[kf])[0]/np.pi
        ctherm=-6*vf*(f-e0)/(np.pi*T*T)
        thermal_free.append({'T':T,'thermal_c_estimator':float(ctherm)})
    assert abs(thermal_free[-1]['thermal_c_estimator']-1)<1e-4
    out={'second_order_max_error':float(err),'tuned_chain':rows,
         'exact_redesign':{'dimension':1024,'lambda':lam,'scale':scale,'g':.02,
            'intertwiner_max_error':intertwiner,'low_spectrum_max_error':specerr,
            'first_noncode_energy_above_ground':threshold,'thermal_comparison':thermal},
         'logical_ED':sectors,'free_density':float(kf/np.pi),'free_velocity':float(vf),
         'free_thermal_coefficient':thermal_free,
         'limitations':'Finite ED checks algebra and scale only. g=.02 is not certified inside the existential RG radius.'}
    Path('checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
