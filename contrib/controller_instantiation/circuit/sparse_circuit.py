#!/usr/bin/env python3
"""Explicit reversible sparse-port circuit; no history-Hamiltonian claim.

All schedules are functions of the graph alone. Classical-looking registers are
coherent neutral bits. Gates form a unitary circuit on the entire enlarged space;
only its promised initialized subspace is required to encode the Fock input.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import argparse
import json
import math
from typing import Iterator


def width(n: int) -> int:
    return max(1, (n - 1).bit_length())


@dataclass(frozen=True)
class Graph:
    u: int
    v: int
    edges: tuple[tuple[int, int], ...]

    def __post_init__(self):
        assert self.u > 0 and self.v > 0 and self.edges
        assert len(set(self.edges)) == len(self.edges)
        assert all(0 <= a < self.u and 0 <= b < self.v for a, b in self.edges)

    @property
    def q(self): return self.u + self.v
    @property
    def n(self): return len(self.edges)
    @property
    def species_width(self): return width(self.q)
    @property
    def label_width(self): return width(max(self.u, self.n))
    @property
    def m(self): return max(1, self.species_width, 2 * self.label_width)
    @property
    def sigma_dimension(self): return self.q + 3 + 2 * (self.q + 1) * self.m


def projective_plane(p: int) -> Graph:
    """PG(2,p): canonical nonzero vectors, incidence via dot product mod p."""
    if p < 2 or any(p % k == 0 for k in range(2, math.isqrt(p) + 1)):
        raise ValueError('p must be prime')
    vectors = [(1, x, y) for x in range(p) for y in range(p)]
    vectors += [(0, 1, x) for x in range(p)] + [(0, 0, 1)]
    edges = tuple((a, b) for a, x in enumerate(vectors) for b, y in enumerate(vectors)
                  if sum(c * d for c, d in zip(x, y)) % p == 0)
    return Graph(len(vectors), len(vectors), edges)


def flag_word(g: Graph, kind: str, index: int) -> tuple[int, ...]:
    counts = [2] * g.q
    if kind == 'g':
        counts[index] -= 1
    elif kind == 'h':
        a, b = g.edges[index]
        counts[a] += 1
        counts[g.u + b] -= 2
    else:
        raise ValueError(kind)
    return tuple(s for s, count in enumerate(counts) for _ in range(count))


def logical_word(g: Graph, occupations: int) -> tuple[int, ...]:
    return tuple(a if occupations >> j & 1 else g.u + b
                 for j, (a, b) in enumerate(g.edges))


def reference_word(g: Graph) -> tuple[int, ...]:
    return tuple(s for s in range(g.q) for _ in range(2))


def branch_output(g: Graph, occupation: int, kind: int, label: int,
                  right: bool = False) -> tuple[tuple[int, ...], int]:
    """Independent direct Kraus action, returning species word and sign exponent."""
    if not kind:
        return (label,) + flag_word(g, 'g', label) + logical_word(g, occupation), 0
    a, b = g.edges[label]
    occupied = occupation >> label & 1
    sign = (occupation & ((1 << label) - 1)).bit_count() % 2
    if right:
        sign ^= (occupation.bit_count() - occupied) % 2
    word = ((g.u + b,) + flag_word(g, 'h' if occupied else 'g', label if occupied else a)
            + logical_word(g, occupation ^ (1 << label)))
    return word, sign


class Builder:
    """Bennett-style Boolean network, all new wires begin in zero.

    Wires never alias their controls. Temporary gates are reversed explicitly,
    rather than treating an irreversible Python function as a quantum oracle.
    """
    def __init__(self):
        self.wires = 0
        self.ops: list[tuple] = []

    def alloc(self, n=1):
        out = tuple(range(self.wires, self.wires + n))
        self.wires += n
        return out

    def bit(self): return self.alloc()[0]

    def const(self, value, n):
        out = self.alloc(n)
        for j, wire in enumerate(out):
            if value >> j & 1: self.ops.append(('X', wire))
        return out

    def xor(self, *args):
        out = self.bit()
        for a in args: self.ops.append(('CX', a, out))
        return out

    def neg(self, a):
        out = self.xor(a)
        self.ops.append(('X', out))
        return out

    def and_(self, a, b):
        if a == b: return self.xor(a)
        out = self.bit()
        self.ops.append(('CCX', a, b, out))
        return out

    def or_(self, a, b):
        return self.xor(a, b, self.and_(a, b))

    def eq(self, bits, value):
        literals = [a if value >> j & 1 else self.neg(a) for j, a in enumerate(bits)]
        if not literals: return self.const(1, 1)[0]
        result = literals[0]
        for a in literals[1:]: result = self.and_(result, a)
        return result

    def same(self, left, right):
        return self.eq(tuple(self.xor(a, b) for a, b in zip(left, right)), 0)

    def _select_word(self, control, one, zero, n):
        result = []
        for j in range(n):
            r = self.xor(control) if (one ^ zero) >> j & 1 else self.bit()
            if zero >> j & 1: self.ops.append(('X', r))
            result.append(r)
        return tuple(result)

    def conditional_swap(self, control, left, right):
        lo, ro = [], []
        for a, b in zip(left, right):
            delta = self.and_(control, self.xor(a, b))
            lo.append(self.xor(a, delta)); ro.append(self.xor(b, delta))
        return tuple(lo), tuple(ro)

    def reverse_since(self, start):
        block = self.ops[start:].copy()
        if any(op[0] not in {'X', 'CX', 'CCX'} for op in block):
            raise ValueError('Only Boolean compute blocks may use reverse_since')
        self.ops.extend(reversed(block))


def selection_bits(b: Builder, g: Graph, kind, label):
    upper = [b.and_(b.neg(kind), b.eq(label, a)) for a in range(g.u)]
    edge = [b.and_(kind, b.eq(label, e)) for e in range(g.n)]
    return upper, edge


def combine_constants(b: Builder, cases, n):
    return tuple(b.xor(*(control for control, value in cases if value >> j & 1))
                 for j in range(n))


def target_network(b: Builder, g: Graph, x, kind, label, right):
    m = g.species_width
    n = [b.eq(x[j], a) for j, (a, _) in enumerate(g.edges)]
    upper, edge = selection_bits(b, g, kind, label)
    creation = [b.and_(edge[j], b.neg(n[j])) for j in range(g.n)]
    annihilation = [b.and_(edge[j], n[j]) for j in range(g.n)]
    port_cases = [(upper[a], a) for a in range(g.u)]
    port_cases += [(edge[e], g.u + v) for e, (_, v) in enumerate(g.edges)]
    y = [combine_constants(b, port_cases, m)]
    flags_g = [flag_word(g, 'g', a) for a in range(g.u)]
    flags_h = [flag_word(g, 'h', e) for e in range(g.n)]
    for site in range(2 * g.q - 1):
        cases = [(upper[a], flags_g[a][site]) for a in range(g.u)]
        cases += [(creation[e], flags_g[a][site]) for e, (a, _) in enumerate(g.edges)]
        cases += [(annihilation[e], flags_h[e][site]) for e in range(g.n)]
        y.append(combine_constants(b, cases, m))
    for e, (a, v) in enumerate(g.edges):
        y.append(b._select_word(b.xor(n[e], edge[e]), a, g.u + v, m))
    terms = []
    for e in range(g.n):
        parity = b.xor(*n[:e])
        if right: parity = b.xor(parity, *n, n[e])
        terms.append(b.and_(edge[e], parity))
    return y, b.xor(*terms)


def inverse_network(b: Builder, g: Graph, y):
    m = g.species_width
    port, flag, payload = y[0], y[1:2*g.q], y[2*g.q:]
    upper = [b.eq(port, a) for a in range(g.u)]
    flat_flag = tuple(w for item in flag for w in item)
    def flag_value(word):
        return sum(s << (m * j) for j, s in enumerate(word))
    gs = [b.eq(flat_flag, flag_value(flag_word(g, 'g', a))) for a in range(g.u)]
    hs = [b.eq(flat_flag, flag_value(flag_word(g, 'h', e))) for e in range(g.n)]
    edge = [b.and_(b.eq(port, g.u + v), b.xor(gs[a], hs[e]))
            for e, (a, v) in enumerate(g.edges)]
    kind = b.xor(*edge)
    label = combine_constants(b, [(upper[a], a) for a in range(g.u)] +
                              [(edge[e], e) for e in range(g.n)], g.label_width)
    x = []
    for e, (a, v) in enumerate(g.edges):
        n = b.xor(b.eq(payload[e], a), edge[e])
        x.append(b._select_word(n, a, g.u + v, m))
    x += [b.const(s, m) for s in reference_word(g)]
    return x, kind, label


def routing_network(b: Builder, x, target):
    """Fixed bubble schedule; each planned swap is a coherent bit, not a branch."""
    current = list(x)
    plan = []
    for i in range(len(x) - 1):
        found = b.bit()  # zero
        for j in range(len(x) - 1, i, -1):
            found = b.or_(found, b.same(current[j], target[i]))
            current[j-1], current[j] = b.conditional_swap(found, current[j-1], current[j])
            plan.append((found, j-1, j))
    return plan


def rotations(g: Graph):
    result = []
    ell = g.label_width
    for M in (g.u, g.n):
        for depth in range(ell):
            span = 1 << (ell - depth)
            tail = M % span
            if tail:
                A = min(tail, span // 2)
                result.append((Fraction(A, tail), Fraction(tail - A, tail)))
            else:
                result.append((Fraction(1), Fraction(0)))
    return tuple(result)


def prepare_branch(b: Builder, g: Graph, kind, label):
    """An explicit fixed H/CH/indexed-R circuit, with predicates uncomputed."""
    b.ops.append(('H', kind))
    ell = g.label_width
    for branch_kind, M in enumerate((g.u, g.n)):
        for depth in range(ell):
            # Controls are generated for nonempty full prefixes and at most one partial.
            span = 1 << (ell-depth)
            for pref in range((M + span - 1) // span):
                count = min(span, M - pref * span)
                if count <= 0: continue
                start = len(b.ops)
                branch = kind if branch_kind else b.neg(kind)
                # eq expects LSB-first; the selected high register bits already are LSB-first.
                prefix_bits = label[ell-depth:]
                control = b.and_(branch, b.eq(prefix_bits, pref))
                computed = b.ops[start:].copy()
                if count == span:
                    b.ops.append(('CH', control, label[ell-depth-1]))
                else:
                    b.ops.append(('R', branch_kind * ell + depth, control, label[ell-depth-1]))
                b.ops.extend(reversed(computed))


@dataclass
class Stage:
    offset: int
    kind: int
    label: tuple[int, ...]
    prep: tuple[tuple, ...]
    body: tuple[tuple, ...]


@dataclass
class Circuit:
    graph: Graph
    wires: int
    data_sites: int
    stages: tuple[Stage, ...]

    def abstract_ops(self, include_prep=True):
        for stage in self.stages:
            if include_prep: yield from stage.prep
            yield from stage.body

    def manifest(self):
        counts = Counter()
        digest = sha256()
        for op in self.abstract_ops():
            counts[op[0]] += 1
            digest.update((json.dumps(op, separators=(',', ':')) + '\n').encode())
        local_counts = Localizer(self).gate_counts()
        return dict(u=self.graph.u, v=self.graph.v, q=self.graph.q, N=self.graph.n,
                    edges=[list(e) for e in self.graph.edges], neutral_bits=self.wires,
                    charged_sites=self.data_sites, total_sites=self.data_sites+self.wires,
                    local_sigma_dimension=self.graph.sigma_dimension,
                    index_states=self.graph.m, opcode_counts=dict(sorted(counts.items())),
                    abstract_gate_count=sum(counts.values()), abstract_gate_sha256=digest.hexdigest(),
                    local_gate_count=sum(local_counts.values()), local_opcode_counts=dict(sorted(local_counts.items())),
                    exact_rotation_probabilities=[[[x.numerator,x.denominator],[y.numerator,y.denominator]]
                                                  for x,y in rotations(self.graph)],
                    input_contract='All neutral registers are unique BLANK; INIT makes bit0 before the circuit and restores BLANK afterwards.',
                    scope='Externally scheduled circuit. No initialized TI history controller or infrared conclusion.')


def make_circuit(g: Graph, two_ports=True) -> Circuit:
    b = Builder()
    stages = []
    for right in range(2 if two_ports else 1):
        offset = right * 2 * g.q
        length = g.n + 2*g.q
        kind, label = b.bit(), b.alloc(g.label_width)
        start = len(b.ops)
        prepare_branch(b, g, kind, label)
        prep = tuple(b.ops[start:])
        start = len(b.ops)
        x = [b.alloc(g.species_width) for _ in range(length)]
        query_x = [('QUERY', offset+i, bit, wire) for i, bits in enumerate(x) for bit, wire in enumerate(bits)]
        b.ops.extend(query_x)
        compute_start = len(b.ops)
        target, sign = target_network(b, g, x, kind, label, bool(right))
        plan = routing_network(b, x, target)
        compute = b.ops[compute_start:].copy()
        b.ops.extend(('CSWAPD', control, offset+i, offset+j) for control,i,j in plan)
        b.ops.append(('P', sign))
        b.ops.extend(reversed(compute))
        y = [b.alloc(g.species_width) for _ in range(length)]
        query_y = [('QUERY', offset+i, bit, wire) for i,bits in enumerate(y) for bit,wire in enumerate(bits)]
        b.ops.extend(query_y)
        inv_start = len(b.ops)
        recovered_x, recovered_kind, recovered_label = inverse_network(b, g, y)
        inverse_compute = b.ops[inv_start:].copy()
        for src, dst in zip((w for a in recovered_x for w in a), (w for a in x for w in a)):
            b.ops.append(('CX', src, dst))
        b.ops.append(('CX', recovered_kind, kind))
        b.ops.extend(('CX', src, dst) for src,dst in zip(recovered_label,label))
        b.ops.extend(reversed(inverse_compute))
        b.ops.extend(reversed(query_y))
        body = tuple(b.ops[start:])
        stages.append(Stage(offset,kind,label,prep,body))
    return Circuit(g,b.wires,g.n+(4 if two_ports else 2)*g.q,tuple(stages))


def simulate_packed(circuit: Circuit, cases):
    """Execute all Boolean/physical gates with bit-parallel exact arithmetic.

    cases = (input occupation, ((kind,label), ...)); branch amplitudes are
    supplied separately by prepare_branch, which is checked independently.
    """
    g = circuit.graph
    mask = (1 << len(cases))-1
    bits = [0]*circuit.wires
    data = [[0]*g.species_width for _ in range(circuit.data_sites)]
    for case_no,(occupation,branches) in enumerate(cases):
        words = logical_word(g,occupation) + reference_word(g)*len(circuit.stages)
        for site,s in enumerate(words):
            for j in range(g.species_width):
                if s>>j&1: data[site][j] |= 1<<case_no
        for stage,(kind,label) in zip(circuit.stages,branches):
            if kind: bits[stage.kind] |= 1<<case_no
            for j,wire in enumerate(stage.label):
                if label>>j&1: bits[wire] |= 1<<case_no
    phase = 0
    for op in circuit.abstract_ops(include_prep=False):
        name,*a = op
        if name=='X': bits[a[0]] ^= mask
        elif name=='CX': bits[a[1]] ^= bits[a[0]]
        elif name=='CCX': bits[a[2]] ^= bits[a[0]] & bits[a[1]]
        elif name=='P': phase ^= bits[a[0]]
        elif name=='QUERY': bits[a[2]] ^= data[a[0]][a[1]]
        elif name=='CSWAPD':
            control,left,right = a
            for j in range(g.species_width):
                delta=(data[left][j]^data[right][j])&bits[control]
                data[left][j]^=delta;data[right][j]^=delta
        else: raise ValueError(op)
    dirty = 0
    for value in bits: dirty |= value
    outputs=[]
    for i in range(len(cases)):
        word=tuple(sum(((plane>>i)&1)<<j for j,plane in enumerate(site)) for site in data)
        outputs.append((word,(phase>>i)&1,(dirty>>i)&1))
    return outputs


# Fixed numerical/algebraic local opcode families. Every one-site family acts on
# the first site of the pair and trivially on the second.
LOCAL_OPCODES = ('SWAP','X','H','Z','CNOT','CV','CV_DAG','TYPE','SHIFT',
                 'SHIFT_DAG','QUERY','LOAD','EXCHANGE','ROT','ROT_DAG','CH','INIT')
LOCAL_INVERSES = {name: name for name in LOCAL_OPCODES}
for _a, _b in [('CV', 'CV_DAG'), ('SHIFT', 'SHIFT_DAG'), ('ROT', 'ROT_DAG')]:
    LOCAL_INVERSES[_a] = _b
    LOCAL_INVERSES[_b] = _a


class Localizer:
    """Stream nearest-neighbor gates without expanding a large stream in memory."""
    def __init__(self,circuit):
        self.c=circuit
        self.size=circuit.data_sites+circuit.wires

    def pair(self,opcode,left,right):
        assert left!=right
        # Move right operand next to left; routing swaps are undone afterwards.
        route=[]
        if right>left:
            route=list(range(right-1,left,-1))
            location=left
        else:
            # Move the first logical operand left until it is just left of right.
            # For ordered operands, exchange the adjacent pair as well.
            route=list(range(left-1,right-1,-1))
            location=right
        for i in route: yield ('SWAP',i)
        yield (opcode,location)
        for i in reversed(route): yield ('SWAP',i)

    def one(self,opcode,site):
        if site<self.size-1: yield (opcode,site)
        else:
            yield ('SWAP',site-1); yield (opcode,site-1); yield ('SWAP',site-1)

    def indexed_query(self,data,bit,target):
        yield from self.one('TYPE',target)
        for _ in range(bit): yield from self.one('SHIFT',target)
        yield from self.pair('QUERY',data,target)
        for _ in range(bit): yield from self.one('SHIFT_DAG',target)
        yield from self.one('TYPE',target)

    def controlled_data_swap(self,control,left,right):
        # Gather the three operands into an adjacent (data,head,data) block.
        yield from self.one('TYPE',control)
        positions=list(range(self.size))
        route=[]
        for destination,logical in enumerate((left,control,right)):
            current=positions.index(logical)
            for site in range(current-1,destination-1,-1):
                positions[site],positions[site+1]=positions[site+1],positions[site]
                route.append(site)
        for site in route: yield ('SWAP',site)
        yield ('LOAD',0); yield ('EXCHANGE',1); yield ('LOAD',0)
        for site in reversed(route): yield ('SWAP',site)
        yield from self.one('TYPE',control)

    def abstract(self,op):
        name,*a=op
        shift=self.c.data_sites
        if name in {'X','H','P'}:
            yield from self.one({'P':'Z'}.get(name,name),shift+a[0])
        elif name in {'CX','CH'}:
            yield from self.pair({'CX':'CNOT'}.get(name,name),shift+a[0],shift+a[1])
        elif name=='CCX':
            ca,cb,target=(shift+w for w in a)
            for opcode,x,y in [('CV',cb,target),('CNOT',ca,cb),('CV_DAG',cb,target),
                               ('CNOT',ca,cb),('CV',ca,target)]:
                yield from self.pair(opcode,x,y)
        elif name=='QUERY': yield from self.indexed_query(a[0],a[1],shift+a[2])
        elif name=='CSWAPD': yield from self.controlled_data_swap(shift+a[0],a[1],a[2])
        elif name=='R':
            index,control,target=a;control+=shift;target+=shift
            yield from self.one('TYPE',target)
            for _ in range(index): yield from self.one('SHIFT',target)
            yield from self.pair('ROT',target,control)
            for _ in range(index): yield from self.one('SHIFT_DAG',target)
            yield from self.one('TYPE',target)
        else: raise ValueError(op)

    def gate_counts(self):
        """Exact emitted counts without iterating every routing SWAP."""
        counts=Counter()
        def one(opcode,site):
            counts[opcode]+=1
            if site==self.size-1:counts['SWAP']+=2
        def pair(opcode,left,right):
            counts[opcode]+=1
            distance=(right-left-1) if right>left else left-right
            counts['SWAP']+=2*distance
        offset=self.c.data_sites
        for wire in range(self.c.wires):
            one('INIT',offset+wire);one('INIT',offset+wire)
        for name,*a in self.c.abstract_ops():
            if name in {'X','H','P'}:one({'P':'Z'}.get(name,name),offset+a[0])
            elif name in {'CX','CH'}:pair({'CX':'CNOT'}.get(name,name),offset+a[0],offset+a[1])
            elif name=='CCX':
                ca,cb,target=(offset+w for w in a)
                for opcode,x,y in [('CV',cb,target),('CNOT',ca,cb),('CV_DAG',cb,target),
                                   ('CNOT',ca,cb),('CV',ca,target)]:pair(opcode,x,y)
            elif name=='QUERY':
                data,index,target=a;target+=offset
                one('TYPE',target);one('TYPE',target)
                for _ in range(index):one('SHIFT',target);one('SHIFT_DAG',target)
                pair('QUERY',data,target)
            elif name=='CSWAPD':
                control,left,right=a;control+=offset
                one('TYPE',control);one('TYPE',control)
                chosen=(left,control,right)
                swaps=sum(logical-sum(previous<logical for previous in chosen[:i])
                          for i,logical in enumerate(chosen))
                counts['SWAP']+=2*swaps
                counts['LOAD']+=2;counts['EXCHANGE']+=1
            elif name=='R':
                index,control,target=a;control+=offset;target+=offset
                one('TYPE',target);one('TYPE',target)
                for _ in range(index):one('SHIFT',target);one('SHIFT_DAG',target)
                pair('ROT',target,control)
            else:raise ValueError(name)
        return +counts

    def stream(self):
        for wire in range(self.c.wires): yield from self.one('INIT',self.c.data_sites+wire)
        for op in self.c.abstract_ops(): yield from self.abstract(op)
        for wire in reversed(range(self.c.wires)):
            yield from self.one('INIT',self.c.data_sites+wire)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--p',type=int,default=2)
    parser.add_argument('--manifest',default='-')
    parser.add_argument('--abstract-jsonl')
    parser.add_argument('--local-jsonl')
    args=parser.parse_args()
    circuit=make_circuit(projective_plane(args.p))
    manifest=circuit.manifest()
    text=json.dumps(manifest,indent=2)+'\n'
    if args.manifest=='-': print(text,end='')
    else:
        with open(args.manifest,'w') as f:f.write(text)
    if args.abstract_jsonl:
        with open(args.abstract_jsonl,'w') as f:
            for op in circuit.abstract_ops():f.write(json.dumps(op)+'\n')
    if args.local_jsonl:
        with open(args.local_jsonl,'w') as f:
            for op in Localizer(circuit).stream():f.write(json.dumps(op)+'\n')


if __name__=='__main__':main()
