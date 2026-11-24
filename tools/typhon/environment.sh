#!/bin/bash
# Sourced by batch.sh. The local setup command also feeds this file to IAS bash -se.
unset PYTHONHOME PYTHONPATH
export PYTHONNOUSERSITE=1
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2025.06

prefix="$HOME/cd-conjecture/envs/numerics-v1"
if [[ ${1:-} == setup ]]; then
    umask 077
    mkdir -p "$(dirname "$prefix")"
    (
        flock -x 9
        if [[ ! -f "$prefix/ready" ]]; then
            python -m venv --system-site-packages "$prefix"
            "$prefix/bin/python" -m pip --disable-pip-version-check install \
                --only-binary=:all: physics-tenpy==1.1.1 cvxpy==1.9.2
            "$prefix/bin/python" -c 'import tenpy, cvxpy, numpy, scipy, sympy, mpmath, matplotlib, threadpoolctl; assert {"CLARABEL", "SCS"} <= set(cvxpy.installed_solvers())'
            "$prefix/bin/python" -m pip freeze --local > "$prefix/overlay.txt"
            touch "$prefix/ready"
        fi
    ) 9>"$prefix.lock"
fi
if [[ ! -f "$prefix/ready" ]]; then
    echo "Run the Typhon setup command to initialize $prefix" >&2
    return 1
fi
source "$prefix/bin/activate"
