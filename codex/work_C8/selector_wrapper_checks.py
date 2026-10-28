"""Checks of the microscopic selector/wrapper composition.

No finite-coupling CFT is computed. All parameters and Hamiltonians are
generated here. Uses only numpy/scipy and the local construction modules.
"""
import itertools
import json
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import eigsh
from arrow_wrapper_checks import local_bond, legal_words, arrow_data
from length_filter_checks import filter_data, plus, tensors
from vacancy_parent import charge_hamiltonian
from mobile_selector import overlap, rod_graph


def sparse_open(words, h, d):
    lookup = {s: j for j, s in enumerate(words)}
    tr = [[(a//d, a % d, h[a, c]) for a in np.flatnonzero(h[:, c])]
          for c in range(d*d)]
    rows, cols, vals = [], [], []
    for col, word in enumerate(words):
        for j in range(len(word)-1):
            for x, y, z in tr[word[j]*d+word[j+1]]:
                target = word[:j]+(x, y)+word[j+2:]
                assert target in lookup
                rows.append(lookup[target]); cols.append(col); vals.append(z)
    return sp.coo_matrix((vals, (rows, cols)), shape=(len(words),)*2).tocsr()


def selector_parameters(k, U, v):
    b = k+1
    a = lambda l: 4*np.sin(np.pi/(2*np.asarray(l)))**2
    ec = a(b+1)-a(b)
    muc = 2-((b+1)*a(b)-b*a(b+1))
    eta = min((a(b-1)-2*a(b)+a(b+1))/2,
              (a(b)-2*a(b+1)+a(b+2))/3)
    t = 2*U*(1/eta+b/muc+1/abs(ec))
    e1 = charge_hamiltonian(k, 1, U, v)[3][0]
    e2 = charge_hamiltonian(k, 2, U, v)[3][0]
    mu = t*muc+4*U-(b+1)*e1+b*e2
    e = t*ec+e2-e1
    lengths = np.array([l for l in range(2, 10001) if l not in (b, b+1)])
    dl = a(lengths)-a(b)-(lengths-b)*ec
    ratio = dl/(1+np.abs(lengths-b))
    assert np.min(ratio) >= eta-1e-13 and mu > 0 and e < 0
    return dict(b=b,t=float(t),mu=float(mu),energy_density=float(e),
                eta=float(eta),K_selector=float(4*t+9*U),
                wrong_length_margin=float(t*eta-U),
                e1=float(e1),e2=float(e2),
                sampled_min_convexity_ratio=float(np.min(ratio)))


def physical_bond(wrap, par, u):
    """Track-grouped matrix, then permuted to physical site order."""
    s = int(round(np.sqrt(wrap.shape[0])))
    arrows = [(1,0),(0,0),(0,1),(1,1)]
    legal = lambda x,y: arrows[x][1] == arrows[y][0]
    q = np.diag([not legal(x,y) for x in range(4) for y in range(4)])
    pin = np.diag([legal(x,y) and y != 0 for x in range(4) for y in range(4)])
    m = np.diag([int(x==0)+int(y==0) for x in range(4) for y in range(4)])
    hop = np.zeros((16,16)); hop[1*4+2,2*4+3] = hop[2*4+3,1*4+2] = 1
    base_sel = par['K_selector']*q+par['mu']*m/2-par['t']*hop
    base = sp.kron(base_sel,sp.eye(s*s))+sp.kron(pin,sp.csr_matrix(wrap))
    move = np.zeros((16,16)); move[3*4+0,0*4+1] = 1
    convert = np.zeros((s,s)); convert[2,0] = 1
    term = sp.kron(move,sp.kron(convert,sp.eye(s)))
    total = (base-u*(term+term.getH())).tocsr()
    perm = np.array([((a*4+b)*s+x)*s+y
                     for a in range(4) for x in range(s)
                     for b in range(4) for y in range(s)])
    total = total[perm,:][:,perm]
    err = float(sp.linalg.norm(total-total.getH()))
    assert err < 1e-12
    # Every selector transition preserves all three legality indicators.
    fails = 0
    for old,new in [((1,2),(2,3)),((0,1),(3,0))]:
        for left,right in itertools.product(range(4), repeat=2):
            before = [left,*old,right]; after = [left,*new,right]
            fails += ([legal(*before[j:j+2]) for j in range(3)] !=
                      [legal(*after[j:j+2]) for j in range(3)])
    assert fails == 0
    return dict(site_dimension=4*s,bond_dimension=total.shape[0],
                nonzero_entries=total.nnz,hermiticity_error=err,
                overlapping_selector_legality_failures=fails)


def cell_spectra(k, U, v, par):
    phase,hf,af,bf,epsilon = filter_data(k)
    scale = 2*U/epsilon
    hw,_,_ = local_bond(2,scale*hf,scale*af,scale*bf,U,v,2*U,6*U)
    data=[]
    for length in range(2,7):
        words=legal_words(2,length)
        H=sparse_open(words,hw,9)
        if len(words)<400:
            ev=np.linalg.eigvalsh(H.toarray())[:2]
        else:
            ev=np.sort(eigsh(H,k=2,sigma=-4*U-.1,which='LM',
                            tol=1e-11,return_eigenvectors=False,
                            v0=np.random.default_rng(410+length).normal(size=len(words))))
        assert ev[0]>=-4*U-1e-9
        # An illegal wrapper sector is >=-2U, independently of length.
        lower_full=min(ev[0],-2*U)
        excess_lower=(par['mu']-2*par['t']*np.cos(np.pi/length)
                      +lower_full-par['energy_density']*length)
        expected=0. if length in (k+1,k+2) else par['wrong_length_margin']*(1+abs(length-k-1))
        assert excess_lower >= expected-1e-9
        if length in (k+1,k+2):
            target=-4*U+par['e'+str(length-k)]
            assert abs(ev[0]-target)<1e-9
        data.append(dict(length=length,legal_dimension=len(words),
                         legal_ground_energy=float(ev[0]),legal_gap=float(ev[1]-ev[0]),
                         full_cell_excess_lower_bound=float(excess_lower),
                         proved_excess_lower_bound=float(expected)))
    return hw,data


def quantum_frames(k,p,U,v,phase):
    """Graded qubit-swap core: vacuum and symmetric one-particle state.

    Q=4 combines the two-state payload and two-state helix. No dense
    15**N wavefunction is stored. Coefficients are physical basis entries.
    """
    holes,hc,chi,ev=charge_hamiltonian(k,p,U,v)
    basis = [(0,)*k]+[tuple(int(j==i) for j in range(k)) for i in range(k)]
    frame={}
    for hs,amp in zip(holes,chi):
        occ=[j for j in range(k+p) if j not in hs]
        for bits in basis:
            logical=sum(bits)
            logical_amp=1 if logical==0 else 1/np.sqrt(k)
            for fs in itertools.product(range(2),repeat=k):
                filt=np.prod([plus(j*phase)[f] for j,f in enumerate(fs)])
                occupied=iter([3+4*tag+2*x+f for tag,x,f in
                               zip([0]+[1]*(k-2)+[2],bits,fs)])
                word=tuple((0 if j<occ[0] else 2 if j>occ[-1] else 1)
                           if j in hs else next(occupied) for j in range(k+p))
                vec=frame.setdefault(word,np.zeros(2,complex))
                vec[logical]=amp*filt*logical_amp
    r=1-sum(x*x for hs,x in zip(holes,chi) if 0 in hs)
    return frame,holes,chi,float(r)


def apply_endpoint(frame,op,position):
    out=np.zeros((2,2),complex)
    for word,vector in frame.items():
        for x in np.flatnonzero(op[:,word[position]]):
            target=list(word);target[position]=int(x);target=tuple(target)
            if target in frame:
                out+=np.outer(frame[target].conj(),vector)*op[x,word[position]]
    return out


def quantum_transport(k,U,v,par):
    phase=filter_data(k)[0]
    f1,holes1,chi1,r1=quantum_frames(k,1,U,v,phase)
    f2,holes2,chi2,r2=quantum_frames(k,2,U,v,phase)
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
    porterr=0.;isoerr=0.
    for f,r in [(f1,r1),(f2,r2)]:
        isoerr=max(isoerr,np.linalg.norm(sum(np.outer(z.conj(),z) for z in f.values())-np.eye(2)))
        for op,tag,pos in [(X,0,0),(Y,0,0),(Y,2,-1),(-X,2,-1)]:
            phys=np.zeros((15,15),complex)
            phys[3+4*tag:7+4*tag,3+4*tag:7+4*tag]=np.kron(op,np.eye(2))
            projected=apply_endpoint(f,phys,pos)
            porterr=max(porterr,np.linalg.norm(projected-r*op/np.sqrt(k)))
    # Contract the actual one-site conversion, retaining all four logical
    # two-cell columns. The shifted cut is after k+2 sites.
    converted=np.zeros((4,4),complex);plain=np.zeros((4,4),complex)
    for s,x in f1.items():
        for t,y in f2.items():
            word=s+t;initial=np.kron(x,y)
            aa=f2.get(word[:k+2]);bb=f1.get(word[k+2:])
            if aa is not None and bb is not None:
                plain+=np.outer(np.kron(aa,bb).conj(),initial)
            if word[k+1]==0:
                changed=word[:k+1]+(2,)+word[k+2:]
                aa=f2.get(changed[:k+2]);bb=f1.get(changed[k+2:])
                if aa is not None and bb is not None:
                    converted+=np.outer(np.kron(aa,bb).conj(),initial)
    lookup={h:i for i,h in enumerate(holes2)}
    beta=sum(chi1[i]*chi2[lookup[(h[0],k+1)]] for i,h in enumerate(holes1))**2
    convert_error=np.linalg.norm(converted-beta*np.eye(4))
    hop=overlap(k+1)**2*converted
    hop_error=np.linalg.norm(hop-overlap(k+1)**2*beta*np.eye(4))
    assert max(isoerr,porterr,convert_error,hop_error,np.linalg.norm(plain))<1e-12
    # The selector boundary projector must contain X# as well as R#.
    selector_r_probability=1-(2/(k+1))*np.sin(np.pi*(k)/(k+1))**2
    return dict(payload_dimension=2,wrapper_site_dimension=15,
                composite_site_dimension=60,
                physical_support_sizes=[len(f1),len(f2)],
                isometry_error=float(isoerr),port_error=float(porterr),
                r1=r1,r2=r2,beta=float(beta),
                conversion_error=float(convert_error),
                unconverted_transport_norm=float(np.linalg.norm(plain)),
                marker_hop_coefficient_per_u=float(overlap(k+1)**2*beta),
                combined_marker_hop_error=float(hop_error),
                incorrect_R_only_selector_factor=float(selector_r_probability))


def main():
    k,U,v=3,1.,.1
    par=selector_parameters(k,U,v)
    hw,spectra=cell_spectra(k,U,v,par)
    bond=physical_bond(hw,par,.01)
    quantum=quantum_transport(k,U,v,par)
    graphs=[rod_graph(k+1,M,M//2) for M in (4,6,8)]
    out=dict(k=k,U=U,v=v,parameters=par,physical_bond=bond,
             complete_cell_selection=spectra,quantum_payload=quantum,
             first_order_marker_graphs=graphs,
             limitation='No full-ring diagonalization or finite-coupling CFT estimate.')
    with open('selector_wrapper_checks.json','w') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
