"""Project the actual microscopic hopping columns through explicit cell frames.

This check does not reconstruct the hopping tensor formula. It reads every
column of MicroscopicBond.column(component='hopping') and contracts its actual
entries with independently built endpoint matrix-unit compression maps.
The one-logical-fermion fixture is not a c-d counterexample.
"""
from collections import defaultdict
import copy
import json
import math
from pathlib import Path
import time
import numpy as np

from bond import MicroscopicBond, one_fermion_core, plus
from checks import core_checks, wrapper_states_frame


def wrapper_endpoint_units(bond, core_frame, holes, endpoint):
    """Return V* (|out><in| at one endpoint) V from the full wrapper code.

    The environment is an actual word, so no endpoint-occupation formula,
    port-compression coefficient, or equalizer identity enters this step.
    """
    states, frame = wrapper_states_frame(bond, core_frame, holes)
    end = 0 if endpoint == 'left' else len(states[0])-1
    groups = defaultdict(lambda: np.zeros((bond.wrapper_dimension, frame.shape[1]), complex))
    for state, amplitudes in zip(states, frame):
        environment = state[:end]+state[end+1:]
        groups[environment][state[end]] += amplitudes
    units = np.zeros((bond.wrapper_dimension, bond.wrapper_dimension, 2, 2), complex)
    for slab in groups.values():
        units += np.einsum('am,bn->abmn', slab.conj(), slab)
    identity = np.einsum('aamn->mn', units)
    assert np.linalg.norm(identity-np.eye(2)) < 1e-12
    return units


def selector_endpoint_density(length, endpoint):
    """Partial trace of the explicit sine-wave selector state, independently."""
    words = []
    for head in range(1, length):
        word = (0,) + (1,)*(head-1) + (2,) + (3,)*(length-head-1)
        amplitude = math.sqrt(2/length)*math.sin(math.pi*head/length)
        words.append((word, amplitude))
    end = 0 if endpoint == 'left' else length-1
    groups = defaultdict(lambda: np.zeros(4, complex))
    for word, amplitude in words:
        environment = word[:end]+word[end+1:]
        groups[environment][word[end]] += amplitude
    density = sum(np.outer(slab, slab.conj()) for slab in groups.values())
    assert abs(np.trace(density)-1) < 1e-14
    return density


def endpoint_units(bond, core_frame, holes, endpoint):
    """Full physical-site compression map, indexed by actual bond basis labels."""
    length = bond.k+holes
    wrapper = wrapper_endpoint_units(bond, core_frame, holes, endpoint)
    selector = selector_endpoint_density(length, endpoint)
    position = 0 if endpoint == 'left' else length-1
    helix_a = plus(position*bond.phi)
    helix_b = plus((position-length+1)*bond.phi)
    units = {}
    for out in range(bond.dimension):
        so, wo, ao, bo = bond.decode(out)
        for source in range(bond.dimension):
            si, wi, ai, bi = bond.decode(source)
            # Matrix-unit compression uses conj(output amplitude)*input amplitude.
            scalar = (selector[si,so] * helix_a[ao].conj()*helix_a[ai]
                      * helix_b[bo].conj()*helix_b[bi])
            if scalar == 0 or not np.any(wrapper[wo,wi]):
                continue
            units[out,source] = scalar*wrapper[wo,wi]
    identity = sum(units.get((j,j),np.zeros((2,2))) for j in range(bond.dimension))
    assert np.linalg.norm(identity-np.eye(2)) < 1e-12
    return units


def actual_hopping_entries(bond):
    """Evaluate every physical input column, including all illegal sectors."""
    entries = []
    nonzero_columns = 0
    for column in range(bond.dimension**2):
        values = bond.column(column, 'hopping')
        nonzero_columns += bool(values)
        for row,value in values.items():
            entries.append((row,column,value))
    return entries, nonzero_columns


def project_entries(bond, entries, left_units, right_units, preceding_selector=None):
    result = np.zeros((4,4), complex)
    contributing_entries = 0
    for row,column,value in entries:
        out_left,out_right = divmod(row,bond.dimension)
        in_left,in_right = divmod(column,bond.dimension)
        if preceding_selector is not None and bond.decode(in_left)[0] != preceding_selector:
            continue
        a = left_units.get((out_left,in_left))
        b = right_units.get((out_right,in_right))
        if a is not None and b is not None:
            result += value*np.kron(a,b)
            contributing_entries += 1
    return result,contributing_entries


def main():
    start = time.perf_counter()
    bond = MicroscopicBond(one_fermion_core(),s=1e-7)
    core_frame,_ = core_checks(bond)
    right = {holes:endpoint_units(bond,core_frame,holes,'right') for holes in (1,2)}
    left = {holes:endpoint_units(bond,core_frame,holes,'left') for holes in (1,2)}
    entries, nonzero_columns = actual_hopping_entries(bond)
    raising = np.array([[0,0],[1,0]],complex)
    hop = np.kron(raising,raising.conj().T)+np.kron(raising.conj().T,raising)
    target = -bond.s/(2*math.cos(bond.core.ports[0].fermi_momentum))*hop
    report = dict(status='Direct finite-fixture microscopic projection; no counterexample or infrared claim',
                  local_dimension=bond.dimension,
                  all_input_columns_evaluated=bond.dimension**2,
                  nonzero_hopping_columns=nonzero_columns,
                  nonzero_hopping_entries=len(entries),
                  length_pairs=[])
    # A separate fixture copy without either compensator gives a useful
    # negative control: the raw projected coefficient retains both occupancies.
    raw = copy.copy(bond)
    raw.EL = raw.ER = np.eye(2)
    raw_entries,_ = actual_hopping_entries(raw)
    for lholes in (1,2):
        for rholes in (1,2):
            projected,count = project_entries(bond,entries,right[lholes],left[rholes])
            error = np.linalg.norm(projected-target)/np.linalg.norm(target)
            assert error < 2e-12, (lholes,rholes,projected,target)
            # Test both pieces of P_cut separately, rather than only their sum.
            x_piece,_ = project_entries(bond,entries,right[lholes],left[rholes],2)
            r_piece,_ = project_entries(bond,entries,right[lholes],left[rholes],3)
            density = selector_endpoint_density(bond.k+lholes,'right')
            x_error = np.linalg.norm(x_piece-density[2,2]*target)/np.linalg.norm(target)
            r_error = np.linalg.norm(r_piece-density[3,3]*target)/np.linalg.norm(target)
            assert max(x_error,r_error) < 2e-12
            assert np.linalg.norm(r_piece-target)/np.linalg.norm(target) > .1
            raw_projected,_ = project_entries(raw,raw_entries,right[lholes],left[rholes])
            raw_target = bond.r[lholes-1]*bond.r[rholes-1]*target
            raw_error = np.linalg.norm(raw_projected-raw_target)/np.linalg.norm(target)
            assert raw_error < 2e-12
            report['length_pairs'].append(dict(
                left_length=bond.k+lholes,right_length=bond.k+rholes,
                contributing_actual_matrix_entries=count,
                normalized_hopping_amplitude=float((projected[2,1]/bond.s).real),
                expected_normalized_amplitude=float((target[2,1]/bond.s).real),
                relative_matrix_error=float(error),
                x_cut_weight=float(density[2,2].real),
                r_cut_weight=float(density[3,3].real),
                x_cut_relative_error=float(x_error),r_cut_relative_error=float(r_error),
                without_compensators_amplitude=float((raw_projected[2,1]/bond.s).real),
                raw_occupancy_product=float(bond.r[lholes-1]*bond.r[rholes-1]),
                raw_relative_matrix_error=float(raw_error)))
    report['runtime_seconds']=time.perf_counter()-start
    path=Path(__file__).resolve().with_name('port_projection_results.json')
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
