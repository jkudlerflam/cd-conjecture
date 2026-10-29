"""Local terms: controls from contrib/d2_followup_220_1282/models.py and survey points (9-vectors
in candidates.json, convention of strata.py).  Returns (h, r, c_exact, note, field)."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'd2_followup_220_1282'))
import numpy as np
import strata
from vumps import bond
import models as M0

CAND = os.path.join(HERE, 'candidates.json')


def model(name):
    if os.path.exists(CAND):
        c = json.load(open(CAND))
        if name in c:
            J, f = strata.J_f(c[name]['v'])
            return bond(J, f), 2, None, c[name].get('note', name), f
    h, r, c, note = M0.model(name)
    return h, r, c, note, None
