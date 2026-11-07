# heavy first
tasks = []
for n in (18, 17):
    for m in ['Heisenberg', 'XX', 'XXZ_0.5', 'XXZ_2', 'XXZ_-0.9', 'TFIM_1']:
        tasks.append(f'survey.py {m} {n} {n}')
for n in (13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3):
    for i in range(30):
        tasks.append(f'survey.py random_{i:02d} {n} {n}')
for L in (16, 15, 14, 13):
    for i in range(30):
        tasks.append(f'rings.py random_{i:02d} {L} {L}')
for L in (18, 17):
    for m in ['Heisenberg', 'XX', 'XXZ_0.5', 'XXZ_2', 'XXZ_-0.9', 'TFIM_1']:
        tasks.append(f'rings.py {m} {L} {L}')
print('\n'.join(tasks))
