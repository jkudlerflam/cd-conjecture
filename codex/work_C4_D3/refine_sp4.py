#!/usr/bin/env python3
import json, numpy as np
from d3 import Sector,measure,append
# Both sides of the short-size high-entropy region; representative bulk points.
angles=[0,.015,.03,.06544984695,.10,.1308996939,.16,.1963495408,.2617993878,.3926990817,.5235987756,.7853981634,1.0471975512,1.1780972451,1.308996939,1.374446786,1.5707963268,2.0943951024,3.1415926536,3.4033920414,3.6651914292,3.926990817,4.7123889804]
for L in [8,10,12]:
 sec=Sector('sp4',L,cartan=True)
 for t in angles:
  append('refine_sp4.jsonl',measure(sec,[np.cos(t),np.sin(t)],f'sp4_angle{t:.11f}',8))
