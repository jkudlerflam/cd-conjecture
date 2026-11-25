"""Shared provenance and baseline imports; Python >=3.11, numpy/scipy/TeNPy."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BASE = ROOT / 'contrib/d2_xy_yx_chirality_imps_20261118'
sys.path.insert(0, str(BASE / 'scripts'))


def sha256(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def source_hashes():
    files = sorted(HERE.glob('*.py')) + [BASE / 'scripts' / n for n in
             ('ims.py', 'ofd.py', 'meas_quartet.py', 'ana_quartet.py', 'quartet_table.py')]
    return {str(p.relative_to(ROOT)): sha256(p) for p in files if p.exists()}
