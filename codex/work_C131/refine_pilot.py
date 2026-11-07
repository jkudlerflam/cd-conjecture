from optimize_q import run
for L in (2,3,4):
    for seed in range(12):
        run(L,seed,4000,f'refined_L{L}_s{seed}.json')
