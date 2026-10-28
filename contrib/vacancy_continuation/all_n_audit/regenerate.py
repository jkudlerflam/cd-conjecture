"""Replay C39 certificates in this audit's directory, extending through n=4."""
import importlib.util
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
here = Path(__file__).resolve().parent
source = here.parents[2] / "codex/work_C39/certify.py"
spec = importlib.util.spec_from_file_location("c39_certificate", source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
os.chdir(here)
for n in range(1, 5):
    module.certificate(n)
