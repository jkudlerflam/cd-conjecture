"""C116 full-trace checks. Small mode is local; larger ED requires SLURM.
Every spectrum includes all on-site charge sectors. No sparse low-level truncation.
"""
from pathlib import Path
import argparse
import json
import os
import time
from collections import defaultdict
import numpy as np
from scipy.linalg import eigvalsh
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parent

def logz_xx(length, beta):
    # Jordan-Wigner boundary condition: even occupancy antiperiodic, odd periodic.
    answers = []
    for shift, parity in [(0.5, 0), (0.0, 1)]:
        eps = np.cos(2*np.pi*(np.arange(length)+shift)/length)
        even, odd = 0.0, -np.inf
        for energy in eps:
            even, odd = (np.logaddexp(even, odd-beta*energy),
                         np.logaddexp(odd, even-beta*energy))
        answers.append((even, odd)[parity])
    return float(np.logaddexp(*answers))

def spectrum(model, length):
    d = 3 if model == 'SU3' else 2
    sectors = defaultdict(list)
    powers = d**np.arange(length, dtype=np.int64)
    for state in range(d**length):
        digits = (state // powers) % d
        counts = tuple(int(np.count_nonzero(digits == a)) for a in range(d-1))
        sectors[counts].append(state)
    energies = []
    maxdim = 0
    trace = 0.0
    trace2 = 0.0
    for basis in sectors.values():
        maxdim = max(maxdim, len(basis))
        index = {state: j for j, state in enumerate(basis)}
        h = np.zeros((len(basis), len(basis)))
        for col, state in enumerate(basis):
            digits = (state // powers) % d
            for i in range(length):
                j = (i+1) % length
                a, b = int(digits[i]), int(digits[j])
                if model == 'SU3':
                    swapped = state + (b-a)*int(powers[i]) + (a-b)*int(powers[j])
                    h[index[swapped], col] += 1.0
                else:
                    if model == 'Heisenberg':
                        h[col, col] += (0.25 if a == b else -0.25)
                    if a != b:
                        swapped = state + (b-a)*int(powers[i]) + (a-b)*int(powers[j])
                        h[index[swapped], col] += 0.5
        assert np.max(np.abs(h-h.T)) < 1e-13
        ev = eigvalsh(h, overwrite_a=False, check_finite=False)
        trace += float(np.trace(h))
        trace2 += float(np.sum(h*h))
        energies.extend(ev.tolist())
    energies = np.sort(energies)
    assert len(energies) == d**length
    assert abs(np.sum(energies)-trace) < 1e-8*max(1,abs(trace))
    assert abs(np.dot(energies, energies)-trace2) < 1e-8*max(1,trace2)
    return energies, maxdim

def lz(energies, beta):
    return float(logsumexp(-beta*energies))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cluster', action='store_true')
    args = parser.parse_args()
    if args.cluster and not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Larger exact diagonalizations require a SLURM allocation.')
    start = time.time()
    lengths = {'XX':[4,6,8], 'Heisenberg':[4,6,8], 'SU3':[3,6]}
    if args.cluster:
        lengths = {'XX':[4,6,8,10,12], 'Heisenberg':[4,6,8,10,12], 'SU3':[3,6,9]}
    constants = {'XX':(2,0.5,-1/np.pi,1.0),
                 'Heisenberg':(2,0.5,0.25-np.log(2),np.pi/2),
                 'SU3':(3,1.0,1-np.log(3)-np.pi/(3*np.sqrt(3)),4*np.pi/3)}
    out = {'mode':'cluster' if args.cluster else 'small', 'models':{}, 'xx_analytic':[], 'massive_spectator':[]}
    xx_error = 0.0
    for model, sizes in lengths.items():
        d, m, e_inf, target = constants[model]
        data = {'d':d,'M':m,'e_inf':e_inf,'A_target':target,'rings':[], 'ratios':[]}
        spectra = {}
        for length in sizes:
            ev, blockdim = spectrum(model,length)
            spectra[length] = ev
            data['rings'].append({'L':length,'dimension':len(ev),'max_block':blockdim,
                                  'E0':float(ev[0]), 'A_energy':float(6*length*(e_inf*length-ev[0])/np.pi)})
            for tau in [0.25,0.5,1.0,2.0,4.0,8.0]:
                beta = tau*length
                if model == 'XX':
                    xx_error = max(xx_error,abs(lz(ev,beta)-logz_xx(length,beta)))
        for length in sizes:
            if 2*length not in spectra:
                continue
            ev, ev2 = spectra[length], spectra[2*length]
            for tau in [0.25,0.5,1.0,2.0,4.0,8.0]:
                beta = tau*length
                r = 2*lz(ev,beta)-lz(ev2,beta)
                q = 4*lz(ev,beta)-lz(ev2,2*beta)
                assert abs(r) <= 4*beta*m+1e-9
                data['ratios'].append({'L':length,'tau':tau,'R':r,'Q':q,
                    'A_R':4*r/(np.pi*tau),'A_Q':2*q/(np.pi*tau),
                    'gluing_bound':4*beta*m,'rank_bound_log_d2':2*np.log(d)})
        out['models'][model] = data
    assert xx_error < 1e-10
    out['xx_ED_vs_full_trace_max_error'] = xx_error
    for length in [4,8,16,32,64,128,256]:
        for tau in [0.25,0.5,1.0,2.0,4.0,8.0]:
            beta=tau*length
            r=2*logz_xx(length,beta)-logz_xx(2*length,beta)
            q=4*logz_xx(length,beta)-logz_xx(2*length,2*beta)
            out['xx_analytic'].append({'L':length,'tau':tau,'R':r,'Q':q,
                                      'A_R':4*r/(np.pi*tau),'A_Q':2*q/(np.pi*tau)})
    # Massive ferromagnetic classical Ising spectator, h=-J Z tensor Z.
    # log Z = L log(2 cosh(beta J)) + log(1+tanh(beta J)**L).
    # The extensive part cancels exactly from R, leaving the expression below.
    for jcoupling in [0.5,0.05,0.005]:
        for length in [4,16,64,256,1024]:
            tau=4.0
            t=np.tanh(tau*length*jcoupling)
            r=2*np.log1p(t**length)-np.log1p(t**(2*length))
            out['massive_spectator'].append({'J':jcoupling,'L':length,'tau':tau,
                                            'R_massive':float(r),'bound_log_2':float(np.log(2))})
            assert -1e-12 <= r <= np.log(2)+1e-12
    # A two-site gate is positive as an operator, but its spatial realignment need not be.
    # h=+Z tensor Z. In the realigned diagonal subspace, eigenvalues are
    # 2 cosh(beta) and -2 sinh(beta).
    beta=0.7
    spatial=np.array([[np.exp(-beta),np.exp(beta)],[np.exp(beta),np.exp(-beta)]])
    out['realignment_eigenvalues']=np.linalg.eigvalsh(spatial).tolist()
    assert out['realignment_eigenvalues'][0] < 0
    out['seconds']=time.time()-start
    name='checks_cluster.json' if args.cluster else 'checks_small.json'
    (ROOT/name).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'file':name,'seconds':out['seconds'],'xx_error':xx_error,
                      'rings':{k:len(v['rings']) for k,v in out['models'].items()}},indent=2))

if __name__ == '__main__':
    main()
