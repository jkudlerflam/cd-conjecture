import subprocess, math, os
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
for b, w in [(6, 0.001), (6, 0.003), (6, 0.008), (8, 0.0005), (8, 0.0015), (8, 0.004)]:
    u = math.sqrt(1 + b * b / 4) - w
    log = open(f'fes_local/log_b{b}_w{w}.txt', 'w')
    subprocess.Popen(['timeout', '280', 'python3', 'fes.py', 'stag', '--u', f'{u:.8f}', '--b', str(b), '--chis', '8', '12', '16', '24', '32',
                      '--maxit', '400', '--outdir', 'fes_local'], stdout=log, stderr=log, env=env)
