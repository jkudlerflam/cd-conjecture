#!/usr/bin/env python3
"""Closed integer/array bytecode, its compiler, and a journal-reversible VM.

This is an explicit word RAM, NOT a claim to be a local finite-alphabet TM.
All runtime values are integers. Arrays are addressed heap blocks of integers;
no Python function from program.py is executed and no library algorithm is an
instruction. Integer arithmetic and heap addressing still require tape-level
compilation. The exact bytecode is independent of the input prime.
"""
from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path


BINOPS = {
    ast.Add: "ADD", ast.Sub: "SUB", ast.Mult: "MUL",
    ast.FloorDiv: "DIV", ast.Mod: "MOD", ast.LShift: "SHL",
    ast.RShift: "SHR", ast.BitAnd: "AND", ast.BitXor: "XOR",
    ast.Eq: "EQ", ast.NotEq: "NE", ast.Lt: "LT", ast.LtE: "LE",
    ast.Gt: "GT", ast.GtE: "GE",
}


class Compiler:
    def __init__(self, source: str):
        module = ast.parse(source)
        self.functions = {}
        self.code = []
        for node in module.body:
            if isinstance(node, ast.FunctionDef):
                self.functions[node.name] = {"arguments": [a.arg for a in node.args.args]}
            elif not (isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant)):
                raise ValueError("Only function declarations and a docstring are allowed")
        for node in module.body:
            if not isinstance(node, ast.FunctionDef):
                continue
            self.regs = {name: i for i, name in enumerate(self.functions[node.name]["arguments"])}
            self.functions[node.name]["entry"] = len(self.code)
            self.temp = 0
            for statement in node.body:
                self.statement(statement)
            self.functions[node.name]["registers"] = len(self.regs)
        self.payload = {"functions": self.functions, "code": self.code}

    def emit(self, *op):
        self.code.append(list(op))
        return len(self.code) - 1

    def patch(self, pc):
        self.code[pc][1] = len(self.code)

    def register(self, name):
        if name not in self.regs:
            self.regs[name] = len(self.regs)
        return self.regs[name]

    def temporary(self):
        self.temp += 1
        return self.register("$" + str(self.temp))

    def expression(self, node):
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            self.emit("PUSH", int(node.value))
        elif isinstance(node, ast.Name):
            self.emit("GETREG", self.register(node.id))
        elif isinstance(node, ast.List):
            for item in node.elts:
                self.expression(item)
            self.emit("NEW", len(node.elts))
        elif isinstance(node, ast.Subscript):
            self.expression(node.value)
            self.expression(node.slice)
            self.emit("GET")
        elif isinstance(node, ast.BinOp):
            self.expression(node.left)
            self.expression(node.right)
            self.emit(BINOPS[type(node.op)])
        elif isinstance(node, ast.Compare) and len(node.ops) == 1:
            self.expression(node.left)
            self.expression(node.comparators[0])
            self.emit(BINOPS[type(node.ops[0])])
        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            self.emit("PUSH", 0)
            self.expression(node.operand)
            self.emit("SUB")
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute) and node.func.attr == "append":
                assert len(node.args) == 1 and not node.keywords
                self.expression(node.func.value)
                self.expression(node.args[0])
                self.emit("APPEND")
            elif isinstance(node.func, ast.Name):
                assert not node.keywords
                for argument in node.args:
                    self.expression(argument)
                if node.func.id == "len":
                    assert len(node.args) == 1
                    self.emit("LEN")
                elif node.func.id == "emit":
                    assert len(node.args) == 2
                    self.emit("EMIT")
                elif node.func.id in self.functions:
                    assert len(node.args) == len(self.functions[node.func.id]["arguments"])
                    self.emit("CALL", node.func.id)
                else:
                    raise ValueError("No external function: " + ast.unparse(node))
            else:
                raise ValueError(ast.unparse(node))
        else:
            raise ValueError("Unsupported expression: " + ast.dump(node))

    def statement(self, node):
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name):
                self.expression(node.value)
                self.emit("SETREG", self.register(target.id))
            elif isinstance(target, ast.Subscript):
                self.expression(target.value)
                self.expression(target.slice)
                self.expression(node.value)
                self.emit("SET")
            else:
                raise ValueError(ast.dump(node))
        elif isinstance(node, ast.Expr):
            self.expression(node.value)
            self.emit("POP")
        elif isinstance(node, ast.Return):
            self.expression(node.value)
            self.emit("RET")
        elif isinstance(node, ast.If):
            self.expression(node.test)
            skip = self.emit("JZ", None)
            for statement in node.body:
                self.statement(statement)
            jump = self.emit("JMP", None)
            self.patch(skip)
            for statement in node.orelse:
                self.statement(statement)
            self.patch(jump)
        elif isinstance(node, ast.While):
            assert not node.orelse
            start = len(self.code)
            self.expression(node.test)
            stop = self.emit("JZ", None)
            for statement in node.body:
                self.statement(statement)
            self.emit("JMP", start)
            self.patch(stop)
        elif isinstance(node, ast.For):
            assert isinstance(node.target, ast.Name) and not node.orelse
            assert isinstance(node.iter, ast.Call) and isinstance(node.iter.func, ast.Name)
            assert node.iter.func.id == "range" and not node.iter.keywords
            args = node.iter.args
            start = ast.Constant(0) if len(args) == 1 else args[0]
            stop = args[0] if len(args) == 1 else args[1]
            step = 1 if len(args) < 3 else ast.literal_eval(args[2])
            assert step in (-1, 1)
            register = self.register(node.target.id)
            end_register = self.temporary()
            self.expression(start)
            self.emit("SETREG", register)
            self.expression(stop)
            self.emit("SETREG", end_register)
            check = len(self.code)
            self.emit("GETREG", register)
            self.emit("GETREG", end_register)
            self.emit("LT" if step == 1 else "GT")
            jump = self.emit("JZ", None)
            for statement in node.body:
                self.statement(statement)
            self.emit("GETREG", register)
            self.emit("PUSH", step)
            self.emit("ADD")
            self.emit("SETREG", register)
            self.emit("JMP", check)
            self.patch(jump)
        else:
            raise ValueError("Unsupported statement: " + ast.dump(node))


@dataclass
class Frame:
    registers: list[int]
    return_pc: int


class Machine:
    def __init__(self, payload, entry, arguments, *, journal=False, sink=None):
        self.program = payload
        function = payload["functions"][entry]
        assert len(arguments) == len(function["arguments"])
        self.frames = [Frame(list(arguments) + [0] * (function["registers"] - len(arguments)), -1)]
        self.pc = function["entry"]
        self.stack = []
        self.heap = {}
        self.next_address = 1
        self.output = []
        self.sink = sink
        assert not (journal and sink), "A journalled run retains its reversible output track"
        self.journal = [] if journal else None
        self.halted = False
        self.steps = 0
        self.max_stack = 0
        self.max_frames = 1
        self.max_integer_bits = 0
        self.heap_words = 0
        self.emissions = 0
        self.op_counts = Counter()

    def step(self):
        assert not self.halted
        old_pc = self.pc
        op, *args = self.program["code"][self.pc]
        self.pc += 1
        consumed = []
        produced = []
        mutation = None
        def pop(n=1):
            values = self.stack[-n:] if n else []
            if n:
                del self.stack[-n:]
            consumed.extend(values)
            return values
        if op == "PUSH":
            produced = args
        elif op == "GETREG":
            produced = [self.frames[-1].registers[args[0]]]
        elif op == "SETREG":
            value, = pop()
            regs = self.frames[-1].registers
            mutation = ("REG", args[0], regs[args[0]])
            regs[args[0]] = value
        elif op == "NEW":
            values = pop(args[0])
            address = self.next_address
            self.next_address += 1
            self.heap[address] = values.copy()
            self.heap_words += len(values)
            mutation = ("NEW", address)
            produced = [address]
        elif op == "GET":
            address, index = pop(2)
            assert 0 <= index < len(self.heap[address])
            produced = [self.heap[address][index]]
        elif op == "SET":
            address, index, value = pop(3)
            assert 0 <= index < len(self.heap[address])
            mutation = ("SET", address, index, self.heap[address][index])
            self.heap[address][index] = value
        elif op == "APPEND":
            address, value = pop(2)
            self.heap[address].append(value)
            self.heap_words += 1
            mutation = ("APPEND", address)
            produced = [0]
        elif op == "LEN":
            address, = pop()
            produced = [len(self.heap[address])]
        elif op == "JMP":
            self.pc = args[0]
        elif op == "JZ":
            value, = pop()
            if value == 0:
                self.pc = args[0]
        elif op == "CALL":
            function = self.program["functions"][args[0]]
            values = pop(len(function["arguments"]))
            self.frames.append(Frame(values + [0] * (function["registers"] - len(values)), self.pc))
            mutation = ("CALL",)
            self.pc = function["entry"]
        elif op == "RET":
            value, = pop()
            frame = self.frames.pop()
            mutation = ("RET", frame)
            self.pc = frame.return_pc
            self.halted = self.pc == -1
            produced = [value]
        elif op == "POP":
            pop()
        elif op == "EMIT":
            a, b = pop(2)
            if self.sink is None:
                self.output.append((a, b))
            else:
                self.sink(a, b)
            self.emissions += 1
            mutation = ("EMIT",)
            produced = [0]
        elif op in set(BINOPS.values()):
            a, b = pop(2)
            if op == "ADD": value = a + b
            elif op == "SUB": value = a - b
            elif op == "MUL": value = a * b
            elif op == "DIV": value = a // b
            elif op == "MOD": value = a % b
            elif op == "SHL": value = a << b
            elif op == "SHR": value = a >> b
            elif op == "AND": value = a & b
            elif op == "XOR": value = a ^ b
            elif op == "EQ": value = int(a == b)
            elif op == "NE": value = int(a != b)
            elif op == "LT": value = int(a < b)
            elif op == "LE": value = int(a <= b)
            elif op == "GT": value = int(a > b)
            elif op == "GE": value = int(a >= b)
            produced = [value]
        else:
            raise ValueError(op)
        self.stack.extend(produced)
        if self.journal is not None:
            self.journal.append((old_pc, consumed, len(produced), mutation))
        self.steps += 1
        self.op_counts[op] += 1
        self.max_stack = max(self.max_stack, len(self.stack))
        self.max_frames = max(self.max_frames, len(self.frames))
        for value in produced:
            self.max_integer_bits = max(self.max_integer_bits, abs(value).bit_length())

    def undo(self):
        """Literal inverse on every state reached by a journalled forward run."""
        assert self.journal
        old_pc, consumed, produced_count, mutation = self.journal.pop()
        op = self.program["code"][old_pc][0]
        if produced_count:
            del self.stack[-produced_count:]
        if mutation:
            kind, *args = mutation
            if kind == "REG":
                index, old = args
                self.frames[-1].registers[index] = old
            elif kind == "NEW":
                address, = args
                self.heap_words -= len(self.heap[address])
                del self.heap[address]
                self.next_address -= 1
                assert address == self.next_address
            elif kind == "SET":
                address, index, old = args
                self.heap[address][index] = old
            elif kind == "APPEND":
                self.heap[args[0]].pop()
                self.heap_words -= 1
            elif kind == "CALL":
                self.frames.pop()
            elif kind == "RET":
                self.frames.append(args[0])
            elif kind == "EMIT":
                self.output.pop()
                self.emissions -= 1
        self.stack.extend(consumed)
        self.pc = old_pc
        self.halted = False
        self.steps -= 1
        self.op_counts[op] -= 1

    def run(self, limit=None):
        while not self.halted and (limit is None or self.steps < limit):
            self.step()
        return self.stack[-1] if self.halted else None

    def array(self, address):
        return self.heap[address]

    def statistics(self):
        return dict(steps=self.steps, heap_blocks=len(self.heap), heap_words=self.heap_words,
                    max_integer_bits=self.max_integer_bits, max_stack=self.max_stack,
                    max_call_depth=self.max_frames, emissions=self.emissions,
                    bytecode_counts=dict(sorted(self.op_counts.items())))

    def state(self):
        """Semantic state; monotone diagnostic maxima are deliberately excluded."""
        return (self.pc, self.halted, self.next_address, self.stack.copy(),
                [(f.registers.copy(), f.return_pc) for f in self.frames],
                {a: v.copy() for a, v in self.heap.items()}, self.output.copy(),
                self.steps, self.heap_words, self.emissions)


def compile_program():
    return Compiler(Path(__file__).with_name("program.py").read_text()).payload


def program_manifest(payload):
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return dict(instruction_count=len(payload["code"]), function_count=len(payload["functions"]),
                opcode_families=sorted({instruction[0] for instruction in payload["code"]}),
                maximum_registers_per_frame=max(f["registers"] for f in payload["functions"].values()),
                sha256=hashlib.sha256(encoded).hexdigest(),
                scope="Fixed word-RAM program; arithmetic and addressing are not yet local tape transitions")


if __name__ == "__main__":
    payload = compile_program()
    Path(__file__).with_name("bytecode.json").write_text(json.dumps(payload, indent=2) + "\n")
    manifest = program_manifest(payload)
    Path(__file__).with_name("program_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
