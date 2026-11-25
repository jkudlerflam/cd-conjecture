"""Small dense checks only; no many-body optimization required."""
import argparse
import tempfile
import time
import numpy as np
import h5py
from common import source_hashes, write_json
import ims
import ofd
import meas_quartet as old_measure
import measure_matrix
from tenpy.algorithms.exact_diag import ExactDiag
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from tenpy.tools import hdf5_io


def dense_term(ops, L, site_ops):
    out = np.ones((1, 1), complex)
    for j in range(L):
        local = np.eye(2, dtype=complex)
        for name, k in ops:
            if j == k:
                local = local @ site_ops[name]
        out = np.kron(out, local)
    return out


def checks():
    t0 = time.time()
    errors = {}
    L = 6
    P = ims.params({'model': 'target', 'd': '0.1875', 'h': '1.783'})
    model = ofd.fmodel(P, L)
    site = model.lat.mps_sites()[0]
    O = {n: site.get_op(n).to_ndarray() for n in ('Sigmax','Sigmay','Sigmaz','Sp','Sm')}
    H = np.zeros((2**L, 2**L), complex)
    for j in range(L):
        H -= 1.783*dense_term([('Sigmaz',j)], L, O)
    for j in range(L-1):
        for c,a,b in ((19/32,'Sigmax','Sigmay'),(13/32,'Sigmay','Sigmax'),(1,'Sigmaz','Sigmaz')):
            H += c*dense_term([(a,j),(b,j+1)], L, O)
    U = np.ones(1, complex)
    for j in range(L):
        U = np.kron(U, np.exp(-.5j*j*np.pi/2*np.diag(O['Sigmaz'])))
    Hrot = U.conj()[:,None]*H*U[None,:]
    mpo = ExactDiag(model); mpo.build_full_H_from_mpo()
    Hm = mpo.full_H.split_legs().to_ndarray().reshape(2**L,2**L)
    errors['literal_rotated_mpo_max_abs'] = float(np.max(abs(Hm-Hrot)))
    orig = ofd.fmodel({**P,'frame':'orig'}, L)
    ed = ExactDiag(orig); ed.build_full_H_from_mpo()
    Ho = ed.full_H.split_legs().to_ndarray().reshape(H.shape)
    errors['literal_original_mpo_max_abs'] = float(np.max(abs(Ho-H)))
    errors['hermiticity_max_abs'] = float(np.max(abs(H-H.conj().T)))
    errors['field_diagonal_max_abs'] = float(np.max(abs(np.diag(H)-np.diag(Hm))))
    energies = {}
    zprod = np.ones(1)
    for _ in range(L): zprod = np.kron(zprod, np.diag(O['Sigmaz']).real)
    for parity in (0,1):
        ii = np.where(zprod == (-1)**parity)[0]
        energies[str(parity)] = float(np.linalg.eigvalsh(H[np.ix_(ii,ii)])[0])

    # Complex random state in a fixed parity sector exercises all quadratures.
    L = 12
    M = ofd.fmodel(P,L)
    sites = M.lat.mps_sites()
    rng = np.random.default_rng(20261125)
    v = rng.normal(size=2**L) + 1j*rng.normal(size=2**L)
    parity = np.array([i.bit_count()%2 for i in range(2**L)])
    v[parity == 1] = 0; v /= np.linalg.norm(v)
    arr = npc.Array.from_ndarray(v.reshape([2]*L), [s.leg for s in sites],
                                 labels=[f'p{i}' for i in range(L)])
    psi = MPS.from_full(sites,arr,bc='finite',normalize=True)
    def ev(ops):
        w=v.reshape([2]*L)
        for name,j in reversed(ops):
            w=np.moveaxis(np.tensordot(O[name],w,axes=([1],[j])),0,j)
        return np.vdot(v,w.reshape(-1))
    def connected(left,right,i,j):
        a=[(n,i+k) for n,k in left]; b=[(n,j+k) for n,k in right]
        return ev(a+b)-ev(a)*ev(b)
    # Shortened shapes for dense matrix validation, retaining all monomials.
    original_ops = measure_matrix.operators
    test_ops = original_ops()
    for name, offsets in (('Q25',(0,1,2,3)),('Q36',(0,1,3,4))):
        test_ops[name+'m']=[('Sm',i) for i in offsets]
        test_ops[name+'p']=[('Sp',i) for i in offsets]
    measure_matrix.operators=lambda:test_ops
    m=measure_matrix.measure(psi,[0,1],rmin=5,rmax=7,progress=False)
    err=0.; count=0
    for i,row in m['rows'].items():
        C=np.array(row['C']); C=C[...,0]+1j*C[...,1]
        for a,left in enumerate(m['names']):
            for b,right in enumerate(m['names']):
                for t,r in enumerate(row['r']):
                    err=max(err,abs(C[t,a,b]-connected(test_ops[left],test_ops[right],int(i),int(i)+r)))
                    count+=1
    measure_matrix.operators=original_ops
    errors['extended_matrix_dense_max_abs']=float(err)
    q=old_measure.measure(psi,L,[0,1],[(2,1),(1,3)])
    err=0.; nold=0
    for kind,items in (('Q',[(f'{d}_{s}',old_measure.qterm, (d,s)) for d,s in [(2,1),(1,3)]]),
                       ('P',[(str(d),old_measure.pterm,(d,)) for d in (1,2)])):
        for key,maker,args in items:
            for i in (0,1):
                if str(i) not in q[kind][key]: continue
                row=q[kind][key][str(i)]
                for lab,a,b in [('md','Sm','Sp'),('dm','Sp','Sm'),('mm','Sm','Sm'),('dd','Sp','Sp'),('mz','Sm',None)]:
                    left=maker(a,*args)
                    right=maker(b,*args) if b else [('Sigmaz',0)]
                    for t,y in enumerate(row[lab]):
                        err=max(err,abs(complex(*y)-connected(left,right,i,row['j0']+t))); nold+=1
    errors['baseline_correlators_dense_max_abs']=float(err)
    with tempfile.TemporaryDirectory() as d:
        path=d+'/roundtrip.h5'
        with h5py.File(path,'w') as f: hdf5_io.save_to_hdf5(f,{'psi':psi})
        with h5py.File(path,'r') as f: loaded=hdf5_io.load_from_hdf5(f)['psi']
        errors['hdf5_dense_vector_max_abs']=float(np.max(abs(loaded.get_theta(0,L).to_ndarray().reshape(-1)-v)))
        errors['hdf5_energy_abs']=float(abs(M.H_MPO.expectation_value(psi)-M.H_MPO.expectation_value(loaded)))
    assert max(errors.values()) < 2e-11, errors
    return {'passed':True,'errors':errors,'L_hamiltonian':6,'L_dense_vector':L,
            'dense_matrix_entries':count,'baseline_entries':nold,'sector_ground_energies_L6':energies,
            'seed':20261125,'seconds':time.time()-t0,'source_sha256':source_hashes(),
            'matrix_note':'Shortened quartet offsets in the L12 random-state test; production-offset check recorded separately.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',required=True); a=p.parse_args()
    result=checks(); write_json(a.output,result); print(result,flush=True)
