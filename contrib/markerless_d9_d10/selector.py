"""Cell-length selector of the markerless chain at zero port.  A cell with n data sites (length b = n + 1) has
energy nu - 2cos(pi/(b+1)); the midpoint nu = 2 - (delta_b + delta_{b-1})/2 selects b uniquely."""
import mpmath as mp, json
mp.mp.dps = 40
a = lambda x: 2*(1 - mp.cos(mp.pi/(x+1)))
delta = lambda x: (x+1)*a(x) - x*a(x+1)
out = []
for n in (9, 10, 11, 12, 14, 16, 20):
    b = n + 1
    nu = 2 - (delta(b) + delta(b-1))/2
    e = lambda x: (nu - 2 + a(x))/x
    comp = min(range(1, 60*b), key=lambda x: mp.inf if x == b else e(x))
    rec = dict(n=n, cell_length=b, nu=mp.nstr(nu, 20), e_b=mp.nstr(e(b), 12), density_gap=mp.nstr(e(comp) - e(b), 8), nearest_competitor=comp,
               clock_gap=mp.nstr(2*mp.cos(mp.pi/(b+1)) - 2*mp.cos(2*mp.pi/(b+1)), 8))
    assert e(b) < 0 and e(comp) > e(b)
    out.append(rec); print(rec)
json.dump(out, open('data/selector.json', 'w'), indent=1)
