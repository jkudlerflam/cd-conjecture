"""Calibration of degen_opt.W at the two known points: Fibonacci (ZZ + 2Z field) and XZ."""
import numpy as np, degen_opt as d
fib = np.zeros(12); fib[8] = 1.0; fib[11] = 2.0     # J_zz = 1, B_z = 2 (n = (1-Z)/2 blockade: cost of adjacent pair)
fib2 = np.zeros(12); fib2[8] = 1.0; fib2[11] = -2.0
xz = np.zeros(12); xz[2] = 1.0                       # J_xz
heis = np.zeros(12); heis[0] = heis[4] = heis[8] = -1.0
for name, p in (("Fib +", fib), ("Fib -", fib2), ("XZ", xz), ("FM Heis", heis)):
    for (K8, K10) in ((47, 123), (8, 16), (9, 11)):
        print(name, (K8, K10), f"W = {d.W(p, {8: K8, 10: K10}):.3e}")
