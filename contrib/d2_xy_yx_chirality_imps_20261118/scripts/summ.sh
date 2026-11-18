#!/bin/bash
# compact table: run, chi, E, S, xi, m(kappa), B, c_loc(xi), dlnm/dlnxi
cd ~/cc_spinchain/v3
PYTHONPATH= ~/cc_spinchain/venv/bin/python -s an2.py "$@" 2>/dev/null | sed -E 's/Sc=nan //; s/c_conf= *nan //; s/c_conf= *//; s/E=-1\.0([0-9]{9})[0-9]*/E=-1.0\1/; s/xi0= *[0-9.]+ xi1= *[0-9.]+ //; s/c_xi0= *[-0-9.]* //; s/Z=[0-9.-]+ //; s/ +/ /g'
