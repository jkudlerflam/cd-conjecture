from core import *
out=[]
for name,h,ns in [('XX_field',xxz(0,.5),[2,3,4,5]),('Heisenberg',xxz(1),[4,5,6,7,8])]:
    for n in ns:
        rec,r,x=LTI(n,True,'u1').solve(h)
        out.append(dict(name=name,**rec)); print(json.dumps(out[-1]),flush=True)
        Path('matrices').mkdir(exist_ok=True);np.savez(f'matrices/{name}_{n}.npz',rho=r,X=x,h=h)
        write_json('benchmark.json',out)
