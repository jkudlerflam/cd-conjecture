#!/usr/bin/env python3
"""A symbolic instruction bound, obtained from the restricted source AST.

Every loop is pessimistically allowed B iterations, including loops that are
actually constant or logarithmic. Both arms of each conditional are charged.
No empirical fitting enters the polynomial. README gives the input-specific
proof that B=128 A suffices for every loop of this particular program.
"""
import ast
from functools import lru_cache
import json
from pathlib import Path


def add(*polynomials):
    out = [0] * max(map(len, polynomials), default=1)
    for poly in polynomials:
        for i, coefficient in enumerate(poly):
            out[i] += coefficient
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def times_B(poly):
    return (0,) + tuple(poly)


def evaluate(poly, B):
    value = 0
    for coefficient in reversed(poly):
        value = value * B + coefficient
    return value


def analyze():
    module = ast.parse(Path(__file__).with_name("program.py").read_text())
    functions = {node.name: node for node in module.body if isinstance(node, ast.FunctionDef)}

    def expression(node):
        if isinstance(node, (ast.Constant, ast.Name)):
            return (1,)
        if isinstance(node, ast.List):
            return add((1,), *(expression(x) for x in node.elts))
        if isinstance(node, ast.Subscript):
            return add((1,), expression(node.value), expression(node.slice))
        if isinstance(node, ast.BinOp):
            return add((1,), expression(node.left), expression(node.right))
        if isinstance(node, ast.Compare):
            return add((1,), expression(node.left), expression(node.comparators[0]))
        if isinstance(node, ast.UnaryOp):
            return add((2,), expression(node.operand))
        if isinstance(node, ast.Call):
            arguments = add(*(expression(x) for x in node.args))
            if isinstance(node.func, ast.Attribute):
                return add((1,), expression(node.func.value), arguments)
            name = node.func.id
            return add((1,), arguments, function(name) if name in functions else (0,))
        raise ValueError(ast.dump(node))

    def statements(nodes):
        return add(*(statement(node) for node in nodes))

    def statement(node):
        if isinstance(node, ast.Assign):
            target = node.targets[0]
            result = add((1,), expression(node.value))
            if isinstance(target, ast.Subscript):
                result = add(result, expression(target.value), expression(target.slice))
            return result
        if isinstance(node, (ast.Expr, ast.Return)):
            return add((1,), expression(node.value))
        if isinstance(node, ast.If):
            return add((2,), expression(node.test), statements(node.body), statements(node.orelse))
        if isinstance(node, ast.While):
            test = add((1,), expression(node.test))
            return add(test, times_B(add(test, statements(node.body), (1,))))
        if isinstance(node, ast.For):
            args = node.iter.args
            start = ast.Constant(0) if len(args) == 1 else args[0]
            stop = args[0] if len(args) == 1 else args[1]
            setup = add((6,), expression(start), expression(stop))
            return add(setup, times_B(add((9,), statements(node.body))))
        raise ValueError(ast.dump(node))

    @lru_cache(None)
    def function(name):
        return statements(functions[name].body)

    graph = {
        name: sorted({node.func.id for node in ast.walk(body)
                      if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                      and node.func.id in functions})
        for name, body in functions.items()
    }
    @lru_cache(None)
    def depth(name):
        return 1 + max((depth(child) for child in graph[name]), default=0)
    poly = function("main")
    assert all(c >= 0 for c in poly)
    return dict(scope="Source-derived word-instruction bound, not local-TM transitions",
                loop_trip_bound="B = 128 A, A = (N+4q)^2 max(1,ceil(log2 q),ceil(log2 max(u,N)))",
                T_polynomial_coefficients_in_ascending_powers_of_B=list(poly),
                polynomial_degree=len(poly)-1,
                coefficient_sum=sum(poly),
                simpler_bound=f"T <= {sum(poly)} B^{len(poly)-1} for B >= 1",
                maximum_call_depth=max(map(depth, functions)),
                function_polynomials={name:list(function(name)) for name in sorted(functions)},
                largest_literal_array_arity=max(len(n.elts) for n in ast.walk(module) if isinstance(n, ast.List)),
                largest_function_argument_count=max(len(n.args.args) for n in functions.values()),
                journal_words_per_step_upper_bound=64,
                working_and_history_words_upper_bound="100 T + 40000",
                serialized_word_bit_width_upper_bound="b = 2 A + ceil(log2(T+1)) + 64",
                serialized_data_and_journal_bits_upper_bound="(100 T + 40000) (2 b + 4)")


if __name__ == "__main__":
    result = analyze()
    Path(__file__).with_name("resource_bounds.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k:v for k,v in result.items() if k != "function_polynomials"}, indent=2))
