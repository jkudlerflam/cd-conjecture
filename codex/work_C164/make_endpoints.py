"""C164: collect qtune endpoints for the DMRG stage. Keeps endpoints whose final state passes qtune.quality_ok,
sorted by the final common gap relative to the bandwidth; writes the best `top` (default 48) to endpoints.jsonl.
usage: make_endpoints.py out.jsonl top qtune_files..."""
import sys, json
import qtune
out, top, files = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
E = []
for fn in files:
    for l in open(fn):
        r = json.loads(l); inf = r['hist'][-1]['info']
        if qtune.quality_ok(inf):
            E.append((0.5 * (inf['EA'] + inf['EB']) / inf['bw'], r))
E.sort(key=lambda t: t[0])
with open(out, 'w') as fo:
    for g, r in E[:top]:
        fo.write(json.dumps(r) + '\n')
print(len(E), 'quality endpoints;', min(top, len(E)), 'written; gaps/bw', [round(g, 4) for g, _ in E[:top]])
