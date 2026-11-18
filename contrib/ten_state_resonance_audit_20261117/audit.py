"""Reproduce bounded literal-bond, vertex, root and one-loop checks."""
import itertools,json,time,os
from pathlib import Path
import numpy as np
from mechanism import scalar,pair_coefficients,normal_order,projected_controls,fermions,resonant_seed,fermi_modes,vertex,bloch
from literal import parts,cell,literal_coefficients,CHARGE,MARKER


def rg_record(T,W,center,radius):
    seed=resonant_seed(T,center,radius);m=fermi_modes(T,seed['M'],seed['cos_roots'])
    v1=m[0]['v'];v2=m[2]['v'];assert min(v1,v2)>0
    f=lambda a,b,c,d:vertex(W['W'],W['pairs'],m,a,b,c,d)
    g1=f(0,3,0,3).real;g2=f(1,2,1,2).real;gx=f(0,1,0,1).real
    gs=f(0,2,0,2).real;J=f(0,3,1,2)
    D=g1/(2*np.pi*v1)+g2/(2*np.pi*v2)-2*gx/(np.pi*(v1+v2))
    B=1/(2*np.pi**2*v1*v2)+2/(np.pi**2*(v1+v2)**2)
    Y=abs(J)*np.sqrt(B)
    pairs=list(itertools.combinations_with_replacement(range(4),2));collisions=[];separations=[]
    for a,b in itertools.combinations(pairs,2):
        delta=abs(np.exp(1j*sum(m[i]['k'] for i in a))-np.exp(1j*sum(m[i]['k'] for i in b)))
        if delta<1e-10:collisions.append([a,b])
        else:separations.append(delta)
    assert collisions==[[(0,3),(1,2)]]
    assert max(x['residual'] for x in m)<1e-12
    return dict(center=center,radius=radius,seed=seed,modes=m,forward_g1=g1,forward_g2=g2,
        forward_gx=gx,same_chiral_g=gs,pair_transfer=J,D_coefficient=D,Y_coefficient=Y,
        stable_margin_coefficient=D-Y,KT_invariant_coefficient=D*D-Y*Y,
        pair_sum_collisions=collisions,other_pair_sum_chord_margin=min(separations))


def encode(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x))


def main():
    assert os.environ.get('OPENBLAS_NUM_THREADS')=='1'
    started=time.perf_counter();out=dict(physical_d=9,workers=1,checks=[],interpretation='Finite frozen-sector coefficient audit and perturbative one-loop evidence; no all-n or ground-sector CFT certificate.')
    q=(CHARGE[:,None]+CHARGE[None,:]).ravel();ml=np.repeat(MARKER,9);mr=np.tile(MARKER,9)
    for name,H in parts().items():
        errors=[np.max(abs(H-H.T)),np.max(abs(H*(q[:,None]-q[None,:]))),np.max(abs(H*(ml[:,None]-ml[None,:]))),np.max(abs(H*(mr[:,None]-mr[None,:])))]
        assert max(errors)==0
        out['checks'].append(dict(part=name,hermiticity_charge_marker_error=max(errors),literal_nonzero_entries=int(np.count_nonzero(H))))
    out['flavours']=[]
    for n in (2,3):
        start=time.perf_counter();a=pair_coefficients(n);b,dim=literal_coefficients(n);c=cell(n);G=c['G'];f=fermions(n)
        left,swap=projected_controls(n);errors=[]
        for name,h in [('number',np.eye(n)),('left',np.diag(left)),('swap',swap)]:
            target=sum(h[i,j]*(f[i].T@f[j]).toarray() for i in range(n) for j in range(n))
            errors.append(float(np.max(abs(G.T@c[name]@G-target))))
        errors.append(float(np.max(abs(c['walk']@G-c['Ev']*G))))
        q2=normal_order(a[2],2*n);q4=normal_order(a[4],2*n);H,R,T,A=scalar(n)
        errors.extend([float(np.max(abs(a[k]-b[k]))) for k in a])
        errors.extend([q2['gaussian_residual'],q4['higher_residual'],np.max(abs(q2['h'][:n,n:]-T)),np.max(abs(q2['h'][:n,:n]+q2['h'][n:,n:]))])
        assert max(errors)<1e-12
        ks=np.linspace(0,np.pi,1001);ev=np.array([np.linalg.eigvalsh(bloch(T,np.zeros_like(T),k)) for k in ks])
        out['flavours'].append(dict(n=n,physical_sites=n+2,cell_code_dimension=1<<n,**dim,
            max_matching_error=max(errors),seconds=time.perf_counter()-start,T=T,
            bare_band_grid_ranges=np.c_[ev.min(0),ev.max(0)],quartic_W=q4['W'],quartic_pairs=q4['pairs']))
    T=scalar(2)[2];W=normal_order(pair_coefficients(2)[4],4)
    out['primary']=rg_record(T,W,1/5,2/3)
    out['neighbourhood']=[rg_record(T,W,c,r) for c in (3/20,1/5,1/4) for r in (3/5,2/3,7/10)]
    assert min(r['stable_margin_coefficient'] for r in out['neighbourhood'])>0
    out['seconds']=time.perf_counter()-started
    (Path(__file__).parent/'audit.json').write_text(json.dumps(out,indent=2,default=encode)+'\n')
    p=out['primary'];print(json.dumps(dict(seconds=out['seconds'],literal_checks='passed',primary={k:p[k] for k in ('D_coefficient','Y_coefficient','stable_margin_coefficient','KT_invariant_coefficient','other_pair_sum_chord_margin')},minimum_sampled_margin=min(r['stable_margin_coefficient'] for r in out['neighbourhood'])),indent=2,default=encode))


if __name__=='__main__':main()
