#!/usr/bin/env python3
"""Summarize saved quantum runs and compare a fixed-walker variational state."""
from quantum_selector import *

def scalar_energy(pattern,holes,mu,J):
    N=sum(pattern);gaps,ix,Hg,Ec,alpha=charge(N,holes)
    source=np.zeros(N);start=0
    for length in pattern:
        v=np.sqrt(2/(length+1))*np.sin(math.pi*np.arange(1,length+1)/(length+1))
        for r in range(length-1):source[start+r]-=2*v[r]*v[r+1]
        marker=start+length-1
        source[(marker-1)%N]+=mu/2;source[marker]+=mu/2
        start+=length
    potential=(np.asarray(gaps)==0)@source
    H=Hg+J*sp.diags(potential)
    e,vec=eigsh(H,k=1,which='SA',v0=np.ones(H.shape[0]),tol=1e-12)
    return float(e[0]),float(np.linalg.norm(H@vec[:,0]-e[0]*vec[:,0])),source.tolist()

def main():
    allruns=[json.loads((HERE/f'{s}_results.json').read_text()) for s in ('small','large','extended')]
    records=[r for x in allruns for r in x['results']];comparisons=[r for x in allruns for r in x['comparisons']]
    wrong=[r for r in comparisons if r['wrong_cells']]
    scalar=[]
    for r in records:
        if (r['selected_b'],r['N'],tuple(r['pattern'])) not in [(2,8,(2,2,2,2)),(3,6,(3,3)),(3,9,(3,3,3)),(3,9,(2,3,4)),(3,9,(9,))]:continue
        for x in r['sweep']:
            if x['J'] not in (.1,.01,.002):continue
            es,res,source=scalar_energy(tuple(r['pattern']),r['L']-r['N'],r['mu_selector'],x['J'])
            scalar.append(dict(b=r['selected_b'],N=r['N'],pattern=r['pattern'],J=x['J'],quantum_energy=x['energy'],fixed_walker_energy=es,quantum_minus_fixed= x['energy']-es,quantum_minus_fixed_over_J2=(x['energy']-es)/x['J']**2,scalar_residual=res,scalar_sources=source))
    summary=dict(sectors=len(records),quantum_diagonalizations=sum(len(r['sweep']) for r in records),max_dimension=max(r['dimension'] for r in records),max_residual=max(x['residual'] for r in records for x in r['sweep']),negative_sector_excess_count=sum(r['exact_excess']<-1e-10 for r in wrong),minimum_wrong_sector_gap=min(wrong,key=lambda r:r['exact_excess']),most_negative_remainder_over_J2=min(wrong,key=lambda r:r['remainder_over_J2']),largest_remainder_over_J2=max(wrong,key=lambda r:r['remainder_over_J2']),total_seconds=sum(x['seconds'] for x in allruns),fixed_walker_comparison=scalar)
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
