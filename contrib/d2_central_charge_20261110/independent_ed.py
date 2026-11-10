"""Load the immutable, independently implemented Pauli-sector ED reference."""
import hashlib,importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent
SHA='f2724018ef05d393cf4476d7689be8e91b9122f5a6c1b5402d651ce8a10ec1d1'
SOURCE=HERE/'source_snapshots'/f'{SHA}_ed_screen.py'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SHA
spec=importlib.util.spec_from_file_location('_cc_independent_ed',SOURCE)
ed=importlib.util.module_from_spec(spec);spec.loader.exec_module(ed)
