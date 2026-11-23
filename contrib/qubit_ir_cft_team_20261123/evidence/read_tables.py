"""Re-fit rounded target ladders and re-evaluate the printed estimator algebra."""
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import mean


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "d2_xy_yx_chirality_imps_20261118" / "tables"
TAG = "tg_d0.18750_h1.78300_rot"
NUMBER = r"[-+\d.eE]+"


def fit(x, y):
    xbar, ybar = mean(x), mean(y)
    norm = sum((v - xbar) ** 2 for v in x)
    slope = sum((a - xbar) * (b - ybar) for a, b in zip(x, y)) / norm
    intercept = ybar - slope * xbar
    residual = sum((b - slope * a - intercept) ** 2 for a, b in zip(x, y))
    return dict(slope=slope, slope_regression_se=math.sqrt(residual / (len(x) - 2) / norm))


def field(line, label):
    return float(re.search(re.escape(label) + r"\s*(" + NUMBER + r")", line)[1])


def target_block(lines, label):
    header = "== " + TAG + " " + label
    i = lines.index(header) + 1
    out = []
    for line in lines[i:]:
        if line.startswith("== "):
            break
        out.append(line)
    return out


def infinite(lines):
    rows = []
    for line in target_block(lines, "_pin"):
        if line.lstrip().startswith("chi="):
            rows.append(dict(chi=int(field(line, "chi=")), S=field(line, " S="),
                             xi=field(line, " xi="), xi_even=field(line, " xi0="),
                             xi_odd=field(line, " xi1="), kappa=field(line, " m=")))
    fits = []
    for lower in (96, 128, 160):
        window = [r for r in rows if r["chi"] >= lower]
        for channel in ("xi", "xi_even"):
            x = [math.log(r[channel]) for r in window]
            entropy = fit([v / 6 for v in x], [r["S"] for r in window])
            chirality = fit(x, [-math.log(abs(r["kappa"])) for r in window])
            fits.append(dict(chi_min=lower, channel=channel, n=len(x), entropy_c=entropy,
                             kappa_pseudoorder_x=chirality))
    assert [r["chi"] for r in rows] == [32, 48, 64, 96, 128, 160, 200, 256]
    return dict(rows=rows, fits=fits)


def open_entropy(lines):
    rows = []
    for line in target_block(lines, "par=0"):
        if not line.lstrip().startswith("L="):
            continue
        length = int(field(line, "L="))
        cuts = [l for l in range(1, length) if abs(l - length / 2) <= length / 8]
        x = mean(math.log(2 * length / math.pi * math.sin(math.pi * l / length)) / 6 for l in cuts)
        rows.append(dict(L=length, Xbar=x, Sbar=field(line, " Sbar="),
                         dSbar=field(line, " dSbar(chi)=")))
    good = [r for r in rows if r["L"] >= 96 and abs(r["dSbar"]) <= .02]
    assert [r["L"] for r in good] == [96, 128, 192, 256, 384]
    estimates = []
    for lower in (96, 128, 192):
        window = [r for r in good if r["L"] >= lower]
        estimates.append(dict(L_min=lower, lengths=[r["L"] for r in window],
                              **fit([r["Xbar"] for r in window], [r["Sbar"] for r in window])))
    local_slopes = [(b["Sbar"] - a["Sbar"]) / (b["Xbar"] - a["Xbar"]) for a, b in zip(good, good[1:])]
    return dict(rows=rows, fits=estimates, local_slopes_96_through_384=local_slopes)


def correlators(lines):
    rows, odd = [], None
    for line in lines:
        if line.startswith(TAG + "_L"):
            length = int(re.search(r"_L(\d+)_", line)[1])
            eta, pair, sf = [field(line, x) for x in (" eta_b=", " pair3=", " Kc=")]
            err_b, err_p, err_sf = [float(re.search(re.escape(x) + NUMBER + r"\((" + NUMBER + r")\)", line)[1])
                                    for x in (" eta_b=", " pair3=", " Kc=")]
            rows.append(dict(L=length, eta_b=eta, estimator_err_b=err_b, eta_P3=pair,
                             estimator_err_P3=err_p, eta_P2=field(line, " pair2="),
                             Kc_sf=sf, estimator_err_Kc_sf=err_sf,
                             pair_inferred_Kc=1/pair, Ks_pair_formula=1/(4*eta - pair),
                             conjugate_variance_product_if_matched=sf*pair,
                             quartet_dimension_if_matched=2*pair))
        if line.startswith("d=0.18750 h=1.78300 "):
            odd = [dict(L=int(length), chi=int(chi), odd_minus_even=float(diff))
                   for length, chi, diff in re.findall(r"L=(\d+) chi=(\d+) dE=(" + NUMBER + r")", line)]
    assert sorted(r["L"] for r in rows) == [96, 128, 192, 256, 384]
    assert [r["L"] for r in odd] == [64, 96, 128]
    return dict(rows=sorted(rows, key=lambda r: r["L"]), OBC_odd_even=odd)


def main():
    files = ("target_ladders.txt", "open_chain_ladders.txt", "correlators_odd_even.txt")
    data = [(SOURCE / file).read_bytes() for file in files]
    lines = [b.decode().splitlines() for b in data]
    result = dict(source_sha256={f: hashlib.sha256(b).hexdigest() for f, b in zip(files, data)},
                  iMPS=infinite(lines[0]), OBC=open_entropy(lines[1]), exponents=correlators(lines[2]),
                  scope="Rounded target table regressions and estimator arithmetic. Raw target correlators unavailable.")
    (HERE / "table_readings.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
