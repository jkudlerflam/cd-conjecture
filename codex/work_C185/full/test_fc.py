"""Validate the C185 copy of fc.py with the DM term against exact diagonalization of the open chain (L = 10)."""
import json, subprocess, numpy as np, sys
sys.path.insert(0, '.')
import band as B
L = 10; lam = 0.3; jxpy = -0.045; dhz = -0.045
subprocess.run(['python3', 'fc.py', 'T', str(L), '--r', '8', '--chis', '16', '32', '64', '64', '--dmy', str(lam), '--jxpy', str(jxpy), '--dhz', str(dhz), '--outdir', 'test', '--nocorr'], check=True, capture_output=True)
d = json.load(open('test/T_10.json'))
cp = d['couplings']
h = B.bond(cp['Jx'], cp['Jy'], cp['D'], cp['Jz'], cp['hz']) + lam*(np.kron(B.SX, B.SZ) - np.kron(B.SZ, B.SX))
I2 = np.eye(2)
def emb(op2, j):
    return np.kron(np.kron(np.eye(2**j), op2), np.eye(2**(L-j-2)))
H = sum(emb(h, j) for j in range(L-1))
H += np.kron(0.5*cp['hz']*B.SZ, np.eye(2**(L-1))) + np.kron(np.eye(2**(L-1)), 0.5*cp['hz']*B.SZ)
E = np.linalg.eigvalsh(H)[0]
print(f"DMRG E = {d['E']:.12f}, ED E = {E:.12f}, diff = {d['E']-E:.2e}, trunc = {d['trunc']:.1e}")
