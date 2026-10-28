"""Source in the deliberately small integer/array language compiled by machine.py.

This file is DATA to the compiler. No function in this file is called by Python
when generating a circuit. Gate numbers are documented in README.md.
"""


def clone(a, lo, hi):
    r = []
    for i in range(lo, hi):
        r.append(a[i])
    return r


def extend(a, b):
    for i in range(len(b)):
        a.append(b[i])
    return 0


def reverse_extend(a, b):
    for i in range(len(b) - 1, -1, -1):
        a.append(b[i])
    return 0


def width(n):
    r = 1
    while (1 << r) < n:
        r = r + 1
    return r


def graph(u, v, edges):
    n = len(edges)
    z = u
    if n > z:
        z = n
    return [u, v, edges, width(u + v), width(z), u + v, n]


def plane(p):
    vectors = []
    for x in range(p):
        for y in range(p):
            vectors.append([1, x, y])
    for x in range(p):
        vectors.append([0, 1, x])
    vectors.append([0, 0, 1])
    edges = []
    for a in range(len(vectors)):
        for b in range(len(vectors)):
            product = 0
            for j in range(3):
                product = product + vectors[a][j] * vectors[b][j]
            if product % p == 0:
                edges.append([a, b])
    return graph(len(vectors), len(vectors), edges)


def alloc(b, n):
    out = []
    for i in range(n):
        out.append(b[0])
        b[0] = b[0] + 1
    return out


def bit(b):
    out = b[0]
    b[0] = out + 1
    return out


def const(b, value, n):
    out = alloc(b, n)
    for j in range(n):
        if (value >> j) & 1:
            b[1].append([1, out[j]])
    return out


def xor_bits(b, args):
    out = bit(b)
    for i in range(len(args)):
        b[1].append([2, args[i], out])
    return out


def neg(b, a):
    out = xor_bits(b, [a])
    b[1].append([1, out])
    return out


def and_bits(b, a, c):
    if a == c:
        return xor_bits(b, [a])
    out = bit(b)
    b[1].append([3, a, c, out])
    return out


def or_bits(b, a, c):
    t = and_bits(b, a, c)
    return xor_bits(b, [a, c, t])


def eq_bits(b, bits, value):
    literals = []
    for j in range(len(bits)):
        a = bits[j]
        if ((value >> j) & 1) == 0:
            a = neg(b, a)
        literals.append(a)
    if len(literals) == 0:
        out = const(b, 1, 1)
        return out[0]
    result = literals[0]
    for j in range(1, len(literals)):
        result = and_bits(b, result, literals[j])
    return result


def same(b, left, right):
    bits = []
    for j in range(len(left)):
        bits.append(xor_bits(b, [left[j], right[j]]))
    return eq_bits(b, bits, 0)


def select_word(b, control, one, zero, n):
    out = []
    for j in range(n):
        if ((one ^ zero) >> j) & 1:
            r = xor_bits(b, [control])
        else:
            r = bit(b)
        if (zero >> j) & 1:
            b[1].append([1, r])
        out.append(r)
    return out


def cswap_words(b, control, left, right):
    lo = []
    ro = []
    for j in range(len(left)):
        a = left[j]
        c = right[j]
        delta = and_bits(b, control, xor_bits(b, [a, c]))
        lo.append(xor_bits(b, [a, delta]))
        ro.append(xor_bits(b, [c, delta]))
    return [lo, ro]


def combine(b, cases, n):
    out = []
    for j in range(n):
        controls = []
        for i in range(len(cases)):
            if (cases[i][1] >> j) & 1:
                controls.append(cases[i][0])
        out.append(xor_bits(b, controls))
    return out


def flag(g, kind, index):
    counts = []
    for s in range(g[5]):
        counts.append(2)
    if kind == 0:
        counts[index] = counts[index] - 1
    else:
        a = g[2][index][0]
        v = g[0] + g[2][index][1]
        counts[a] = counts[a] + 1
        counts[v] = counts[v] - 2
    out = []
    for s in range(g[5]):
        for j in range(counts[s]):
            out.append(s)
    return out


def target(b, g, x, kind, label, right):
    n = []
    for j in range(g[6]):
        n.append(eq_bits(b, x[j], g[2][j][0]))
    upper = []
    for a in range(g[0]):
        branch = neg(b, kind)
        upper.append(and_bits(b, branch, eq_bits(b, label, a)))
    edge = []
    for e in range(g[6]):
        edge.append(and_bits(b, kind, eq_bits(b, label, e)))
    creation = []
    for j in range(g[6]):
        creation.append(and_bits(b, edge[j], neg(b, n[j])))
    annihilation = []
    for j in range(g[6]):
        annihilation.append(and_bits(b, edge[j], n[j]))
    cases = []
    for a in range(g[0]):
        cases.append([upper[a], a])
    for e in range(g[6]):
        cases.append([edge[e], g[0] + g[2][e][1]])
    y = [combine(b, cases, g[3])]
    flags_g = []
    flags_h = []
    for a in range(g[0]):
        flags_g.append(flag(g, 0, a))
    for e in range(g[6]):
        flags_h.append(flag(g, 1, e))
    for site in range(2 * g[5] - 1):
        cases = []
        for a in range(g[0]):
            cases.append([upper[a], flags_g[a][site]])
        for e in range(g[6]):
            cases.append([creation[e], flags_g[g[2][e][0]][site]])
        for e in range(g[6]):
            cases.append([annihilation[e], flags_h[e][site]])
        y.append(combine(b, cases, g[3]))
    for e in range(g[6]):
        y.append(select_word(b, xor_bits(b, [n[e], edge[e]]), g[2][e][0], g[0] + g[2][e][1], g[3]))
    terms = []
    for e in range(g[6]):
        parity = xor_bits(b, clone(n, 0, e))
        if right:
            args = [parity]
            extend(args, n)
            args.append(n[e])
            parity = xor_bits(b, args)
        terms.append(and_bits(b, edge[e], parity))
    return [y, xor_bits(b, terms)]


def flag_value(g, word):
    value = 0
    for j in range(len(word)):
        value = value + (word[j] << (g[3] * j))
    return value


def inverse(b, g, y):
    port = y[0]
    flat_flag = []
    for j in range(1, 2 * g[5]):
        extend(flat_flag, y[j])
    upper = []
    for a in range(g[0]):
        upper.append(eq_bits(b, port, a))
    gs = []
    hs = []
    for a in range(g[0]):
        gs.append(eq_bits(b, flat_flag, flag_value(g, flag(g, 0, a))))
    for e in range(g[6]):
        hs.append(eq_bits(b, flat_flag, flag_value(g, flag(g, 1, e))))
    edge = []
    for e in range(g[6]):
        a = g[2][e][0]
        v = g[2][e][1]
        first = eq_bits(b, port, g[0] + v)
        second = xor_bits(b, [gs[a], hs[e]])
        edge.append(and_bits(b, first, second))
    kind = xor_bits(b, edge)
    cases = []
    for a in range(g[0]):
        cases.append([upper[a], a])
    for e in range(g[6]):
        cases.append([edge[e], e])
    label = combine(b, cases, g[4])
    x = []
    for e in range(g[6]):
        a = g[2][e][0]
        v = g[2][e][1]
        n = xor_bits(b, [eq_bits(b, y[2 * g[5] + e], a), edge[e]])
        x.append(select_word(b, n, a, g[0] + v, g[3]))
    for s in range(g[5]):
        for j in range(2):
            x.append(const(b, s, g[3]))
    return [x, kind, label]


def route(b, x, target_words):
    current = clone(x, 0, len(x))
    plan = []
    for i in range(len(x) - 1):
        found = bit(b)
        for j in range(len(x) - 1, i, -1):
            found = or_bits(b, found, same(b, current[j], target_words[i]))
            pair = cswap_words(b, found, current[j - 1], current[j])
            current[j - 1] = pair[0]
            current[j] = pair[1]
            plan.append([found, j - 1, j])
    return plan


def prepare(b, g, kind, label):
    b[1].append([4, kind])
    for branch_kind in range(2):
        M = g[0]
        if branch_kind:
            M = g[6]
        for depth in range(g[4]):
            span = 1 << (g[4] - depth)
            for pref in range((M + span - 1) // span):
                count = M - pref * span
                if count > span:
                    count = span
                start = len(b[1])
                branch = kind
                if branch_kind == 0:
                    branch = neg(b, kind)
                prefix = clone(label, g[4] - depth, len(label))
                control = and_bits(b, branch, eq_bits(b, prefix, pref))
                computed = clone(b[1], start, len(b[1]))
                if count == span:
                    b[1].append([5, control, label[g[4] - depth - 1]])
                else:
                    b[1].append([6, branch_kind * g[4] + depth, control, label[g[4] - depth - 1]])
                reverse_extend(b[1], computed)
    return 0


def queries(b, g, offset, words):
    result = []
    for i in range(len(words)):
        for j in range(g[3]):
            result.append([7, offset + i, j, words[i][j]])
    extend(b[1], result)
    return result


def circuit(g):
    b = [0, []]
    for right in range(2):
        offset = right * 2 * g[5]
        length = g[6] + 2 * g[5]
        kind = bit(b)
        label = alloc(b, g[4])
        prepare(b, g, kind, label)
        x = []
        for i in range(length):
            x.append(alloc(b, g[3]))
        queries(b, g, offset, x)
        start = len(b[1])
        target_result = target(b, g, x, kind, label, right)
        plan = route(b, x, target_result[0])
        computed = clone(b[1], start, len(b[1]))
        for i in range(len(plan)):
            b[1].append([8, plan[i][0], offset + plan[i][1], offset + plan[i][2]])
        b[1].append([9, target_result[1]])
        reverse_extend(b[1], computed)
        y = []
        for i in range(length):
            y.append(alloc(b, g[3]))
        query_y = queries(b, g, offset, y)
        start = len(b[1])
        recovered = inverse(b, g, y)
        inverse_computed = clone(b[1], start, len(b[1]))
        for i in range(length):
            for j in range(g[3]):
                b[1].append([2, recovered[0][i][j], x[i][j]])
        b[1].append([2, recovered[1], kind])
        for j in range(len(label)):
            b[1].append([2, recovered[2][j], label[j]])
        reverse_extend(b[1], inverse_computed)
        reverse_extend(b[1], query_y)
    return [b[0], g[6] + 4 * g[5], b[1]]


def gate(state, opcode, address):
    if state[1] == 2:
        emit(opcode, address)
    else:
        while state[0] < address:
            emit(18, 0)
            state[0] = state[0] + 1
        while state[0] > address:
            emit(19, 0)
            state[0] = state[0] - 1
        emit(opcode, 0)
    return 0


def one(state, opcode, site, size):
    if site < size - 1:
        gate(state, opcode, site)
    else:
        gate(state, 1, site - 1)
        gate(state, opcode, site - 1)
        gate(state, 1, site - 1)
    return 0


def pair(state, opcode, left, right):
    if right > left:
        for i in range(right - 1, left, -1):
            gate(state, 1, i)
        gate(state, opcode, left)
        for i in range(left + 1, right):
            gate(state, 1, i)
    else:
        for i in range(left - 1, right - 1, -1):
            gate(state, 1, i)
        gate(state, opcode, right)
        for i in range(right, left):
            gate(state, 1, i)
    return 0


def controlled_swap(state, control, left, right, size):
    one(state, 8, control, size)
    positions = []
    route_sites = []
    for i in range(size):
        positions.append(i)
    chosen = [left, control, right]
    for destination in range(3):
        current = 0
        while positions[current] != chosen[destination]:
            current = current + 1
        for site in range(current - 1, destination - 1, -1):
            tmp = positions[site]
            positions[site] = positions[site + 1]
            positions[site + 1] = tmp
            route_sites.append(site)
    for i in range(len(route_sites)):
        gate(state, 1, route_sites[i])
    gate(state, 12, 0)
    gate(state, 13, 1)
    gate(state, 12, 0)
    for i in range(len(route_sites) - 1, -1, -1):
        gate(state, 1, route_sites[i])
    one(state, 8, control, size)
    return 0


def lower(state, op, shift, size):
    code = op[0]
    if code == 1:
        one(state, 2, shift + op[1], size)
    elif code == 4:
        one(state, 3, shift + op[1], size)
    elif code == 9:
        one(state, 4, shift + op[1], size)
    elif code == 2:
        pair(state, 5, shift + op[1], shift + op[2])
    elif code == 5:
        pair(state, 16, shift + op[1], shift + op[2])
    elif code == 3:
        a = shift + op[1]
        b = shift + op[2]
        c = shift + op[3]
        pair(state, 6, b, c)
        pair(state, 5, a, b)
        pair(state, 7, b, c)
        pair(state, 5, a, b)
        pair(state, 6, a, c)
    elif code == 7:
        t = shift + op[3]
        one(state, 8, t, size)
        for i in range(op[2]):
            one(state, 9, t, size)
        pair(state, 11, op[1], t)
        for i in range(op[2]):
            one(state, 10, t, size)
        one(state, 8, t, size)
    elif code == 8:
        controlled_swap(state, shift + op[1], op[2], op[3], size)
    elif code == 6:
        c = shift + op[2]
        t = shift + op[3]
        one(state, 8, t, size)
        for i in range(op[1]):
            one(state, 9, t, size)
        pair(state, 14, t, c)
        for i in range(op[1]):
            one(state, 10, t, size)
        one(state, 8, t, size)
    return 0


def generate(g, mode):
    c = circuit(g)
    if mode < 2:
        return c
    state = [0, mode]
    size = c[0] + c[1]
    for i in range(c[0]):
        one(state, 17, c[1] + i, size)
    for i in range(len(c[2])):
        lower(state, c[2][i], c[1], size)
    for i in range(c[0] - 1, -1, -1):
        one(state, 17, c[1] + i, size)
    return [c[0], c[1], len(c[2])]


def main(p, mode):
    return generate(plane(p), mode)


def tiny(mode):
    return generate(graph(1, 1, [[0, 0]]), mode)
