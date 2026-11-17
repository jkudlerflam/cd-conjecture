#!/usr/bin/env python3
"""Full constrained-space identities, kernels, Schur jets and local witnesses."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits
import leading as Ld
from higher_jets import series,ROOT,A

HERE=Path(__file__).resolve().parent


def err(M):return float(np.max(abs(M.data),initial=0)) if sparse.issparse(M) else float(np.max(abs(M),initial=0))


def main():
    start=time.time();graphs=[];jets=[];witnesses=[];kernels=[]
    with threadpool_limits(limits=1):
        for q,lengths in ((1,(12,16,20,24)),(3,(12,15,18,21))):
            m=Ld.parameters(q)['m']
            for L in lengths:
                states,H=Ld.matrix(L,q);lap,V,diff=Ld.laplacian_and_penalty(L,q,states,H)
                assert err(diff)<1e-10
                off=H-sparse.diags(H.diagonal());ncomp,labels=connected_components(off,directed=False)
                words=[Ld.residue_word(int(s),L,m) for s in states]
                coo=off.tocoo();bad=sum(words[i]!=words[j] for i,j in zip(coo.row,coo.col));assert bad==0
                zero_components=[k for k in range(ncomp) if max(V[labels==k],default=0)<1e-10]
                assert len(zero_components)==L+1
                assert all(len(set(words[i]))<=1 for i in range(len(states)) if labels[i] in zero_components)
                row=dict(q=q,L=L,dimension=len(states),graph_components=ncomp,
                    laplacian_identity_error=err(diff),residue_word_violations=bad,
                    zero_components=len(zero_components),expected_kernel=L+1)
                graphs.append(row);print('graph',row,flush=True)
                if L<=16:
                    for alpha in (0.,1.,1j):
                        _,Ha=Ld.matrix(L,q,states,alpha=alpha)
                        ev=eigvalsh(Ha.toarray());null=int(np.count_nonzero(abs(ev)<1e-8))
                        assert min(ev)>-1e-8 and null==(L+1 if alpha==0 else m+1)
                        kernels.append(dict(q=q,L=L,alpha=[complex(alpha).real,complex(alpha).imag],
                            minimum=float(ev[0]),kernel=null,first_positive=float(ev[null])))
        for q in (1,3):
            for L in (10,12,13):
                for alpha,beta in ((0,0),(1,0),(1j,0),(.3+.4j,.7-.2j)):
                    states,K,_=series(L,q,alpha,beta,order=6)
                    _,H=Ld.matrix(L,q,states,alpha=alpha,shift=False)
                    e=err(K[2]-H);assert e<1e-8 and err(K[1])<1e-8
                    row=dict(q=q,L=L,alpha=[complex(alpha).real,complex(alpha).imag],
                        beta=[complex(beta).real,complex(beta).imag],dimension=len(states),
                        fifth_norm=err(K[1]),sixth_full_matrix_error=e)
                    jets.append(row);print('jet',row,flush=True)
            for L in (12,13,14,15):
                states,K,_=series(L,q,order=7);pos={int(s):i for i,s in enumerate(states)}
                right=complex(K[3][pos[1|(1<<5)],pos[1]])
                left=complex(K[3][pos[1|(1<<(L-5))],pos[1]])
                # Projection onto all pure-residue uniform fixed-N states.
                if L%Ld.parameters(q)['m']==0:
                    cols=[]
                    for n in range(L//Ld.parameters(q)['m']+1):
                        for c in range(Ld.parameters(q)['m'] if n else 1):
                            col=np.array([float(int(s).bit_count()==n and all(j%Ld.parameters(q)['m']==c for j in Ld.bits(int(s),L))) for s in states])
                            if np.linalg.norm(col):cols.append(col/np.linalg.norm(col))
                    U=np.array(cols).T;proj=err(U.conj().T@K[3]@U)
                else:proj=None
                row=dict(q=q,L=L,right_create_at_5=[right.real,right.imag],left_create_at_5=[left.real,left.imag],
                         seventh_hermiticity_error=err(K[3]-K[3].conj().T),zero_manifold_projection=proj)
                witnesses.append(row);print('seventh',row,flush=True)
                assert abs(right+left)<1e-8 and abs(right)>1e-5
    out=dict(graph_checks=graphs,alpha_kernel_checks=kernels,sixth_jet_checks=jets,seventh_witnesses=witnesses,
             runtime_seconds=time.time()-start,source_commit='b29602b9',
             input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [ROOT/'contrib/qubit_cancellation_classification_20261116/algebra.py',ROOT/'contrib/qubit_cancellation_classification_20261116/DILUTE_SIXTH.md']})
    (HERE/'checks.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
