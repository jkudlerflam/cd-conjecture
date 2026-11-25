"""Publish only fetched small outputs; check snapshot/source provenance."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shutil
import tarfile
from common import ROOT,HERE,sha256,write_json


def collect(run):
    src=ROOT/'runs/typhon'/run;dest=HERE/'results'/run;dest.mkdir(parents=True,exist_ok=True)
    assert sha256(src/'snapshot.tar.gz')==(src/'snapshot.sha256').read_text().split()[0]
    hashes={}
    with tarfile.open(src/'snapshot.tar.gz') as tar:
        for m in tar.getmembers():
            if m.isfile():hashes[m.name]=hashlib.sha256(tar.extractfile(m).read()).hexdigest()
    for name in ('manifest.json','snapshot.tar.gz','snapshot.sha256','job-id','accounting.psv','logs','results'):
        p=src/name
        if not p.exists():continue
        files=[p] if p.is_file() else list(p.rglob('*'))
        for f in files:
            if not f.is_file():continue
            assert f.stat().st_size<45_000_000,(f,'oversize')
            assert f.suffix not in ('.h5','.hdf5'),(f,'large state should remain remote')
            target=dest/f.relative_to(src);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,target)
    write_json(dest/'source_sha256.json',hashes)
    print('collected',run,(src/'accounting.psv').read_text().splitlines()[-1])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run');a=p.parse_args();collect(a.run)
