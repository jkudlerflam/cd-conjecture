"""Local subset of the stage-E scan (chi 64), output xz_local/. Usage: python3 local_xz.py i1 i2 ..."""
import os, sys
os.environ['XZ_OUT'] = 'xz_local'
from multiprocessing import Pool
import xz_scan
if __name__ == '__main__':
    P = xz_scan.points()
    idx = [int(a) for a in sys.argv[1:]]
    with Pool(len(idx)) as p:
        for fn in p.imap_unordered(xz_scan.work, [(i, P[i], 64) for i in idx]):
            print('done', fn, flush=True)
