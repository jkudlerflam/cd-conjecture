"""C197: is the collaborator's c ~ 1.45 chain a crossover from the nearby frustration-free point?
Usage: python3 run197.py <D> <hz> <L> [P].  Fixed Jx = 1/2, Jy = -1/2, Jz = 5/8, parity P (default 0, even), chi = 128, seed 17.
Uses the collaborator's unchanged wrapper (contrib/d2_order_one_20261110/dmrg_screen.py) around the C190 solver,
with outputs here and checkpoints outside the repo."""
import sys, os
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / 'contrib/d2_order_one_20261110'))
import dmrg_screen as S
D, hz, L = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
name = f'D{sys.argv[1]}_h{sys.argv[2]}'
point = dict(name=name, couplings=dict(Jx=0.5, Jy=-0.5, Jz=0.625, Dz=D, hz=hz))
P = int(sys.argv[4]) if len(sys.argv) > 4 else 0
S.run_point(point, L, P, chi=128, seed=17, output_dir=HERE / 'out',
            checkpoint_dir=os.environ.get('C197_CKPT', '/tmp/c197_ckpt'))
