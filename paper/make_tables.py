#!/usr/bin/env python3
"""Regenerate paper/tables/*.tex from results/spline_analysis.json
(written by analyse_spline.py).  Run from the repository root:
    python paper/make_tables.py
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
A = json.loads((ROOT / "results" / "spline_analysis.json").read_text())
OUT = ROOT / "paper" / "tables"
OUT.mkdir(exist_ok=True)

r = np.array(A["r"]); V = np.array(A["V"]); sV = np.array(A["sigma"])
par = {k: np.array(v) for k, v in A["parametric"].items()}
C = json.loads((ROOT / "results" / "centroid_analysis.json").read_text())
Vc = np.array(C["V"]); sVc = np.array(C["sigma"])           # centroid fit, already χ²-scaled
REF = json.loads((ROOT / "results" / f"spline_{A['ref']}.json").read_text())
SCALE = np.sqrt(max(REF["chi2_data"] / REF["dof"], 1.0))    # PDG-style scaling of the J-resolved fit errors
sV = sV * SCALE
print(f"J-resolved error scale factor {SCALE:.2f}; centroid {C['scale']:.2f}")


def at(arr, x):
    return float(np.interp(x, r, arr))


# ---- Table 1: V(r) ± σ vs shifted parametric potentials
rows = []
for x in (0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
    rows.append(f"{x:g} & ${at(V, x):+.3f}\\pm{at(sV, x):.3f}$ & ${at(Vc, x):+.3f}\\pm{at(sVc, x):.3f}$ & "
                + " & ".join(f"${at(par[k], x):+.3f}$" for k in ("Martin (running)", "Log (running)", "Cornell")) + r" \\")
(OUT / "potential.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")

# ---- Table 2: stability scan
def key_order(k):
    K = int(k.split("_")[0][1:]); lam = float(k.split("_")[1][3:]); tag = k.split("__")[1] if "__" in k else ""
    return (K, lam, tag)

rows = []
for k in sorted(A["stability"], key=key_order):
    s = A["stability"][k]
    tag = k.split("__")[1] if "__" in k else "--"
    tag = {"noD": "no $D$", "no3P": "no $3P$", "floor5": "floor 5"}.get(tag, tag.replace("kappa", r"$\kappa=") + ("$" if tag.startswith("kappa") else ""))
    ndata = 27 - {"noD": 4, "no3P": 2}.get(k.split("__")[1] if "__" in k else "", 0)
    sc = np.sqrt(max(s["chi2"] / (ndata - s["K"] - 3), 1.0))
    rows.append(f"{s['K']} & {s['lam']:g} & {tag} & {s['chi2']:.1f} & {s['prior']:.1f} & "
                + " & ".join(f"{1000*sc*s['sig'][q]:.0f}" for q in ("0.3", "1", "3", "6")) + " & "
                + " & ".join(f"${s['V'][q]:+.3f}$" for q in ("0.3", "3", "6")) + r" \\")
(OUT / "stability.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")

# ---- Table 3: B_c with fit error and variant half-range
LAB = {"11S0": r"$1\,^1S_0$", "13S1": r"$1\,^3S_1$", "21S0": r"$2\,^1S_0$", "23S1": r"$2\,^3S_1$",
       "13P0": r"$1\,^3P_0$", "1P1": r"$1^{+}$", "1P1'": r"$1^{+'}$", "13P2": r"$1\,^3P_2$",
       "13D1": r"$1\,^3D_1$", "13D3": r"$1\,^3D_3$", "33S1": r"$3\,^3S_1$"}
EXP = {"11S0": (6274.5, 0.3), "13S1": (6339.0, 1.7), "21S0": (6871.2, 1.0)}
core = {k: v for k, v in A["bc_variants"].items() if "lam100" not in k}


def half(f):
    vals = [f(v) for v in core.values()]
    return (max(vals) - min(vals)) / 2


rows = []
for l, (m, e) in A["bc"].items():
    ex = EXP.get(l)
    rows.append(f"{LAB[l]} & ${m:.0f}\\pm{SCALE*e:.0f}\\pm{half(lambda v: v[l]):.0f}$ & "
                + (f"${ex[0]:.1f}\\pm{ex[1]}$" if ex else "--") + r" \\")
h = A["hyp"]
rows.append(f"$\\Delta_{{1S}}$ & ${h['d1S'][0]:.1f}\\pm{SCALE*h['d1S'][1]:.1f}\\pm{half(lambda v: v['13S1']-v['11S0']):.1f}$ & $64.5\\pm1.7$ \\\\")
rows.append(f"$\\Delta_{{2S}}$ & ${h['d2S'][0]:.1f}\\pm{SCALE*h['d2S'][1]:.1f}\\pm{half(lambda v: v['23S1']-v['21S0']):.1f}$ & $35.4\\pm2.4$ \\\\")
rows.append(f"$M(2\\,^1S_0)-M(1\\,^1S_0)$ & ${h['s21'][0]:.0f}\\pm{SCALE*h['s21'][1]:.0f}\\pm{half(lambda v: v['21S0']-v['11S0']):.0f}$ & $596.7\\pm1.0$ \\\\")
(OUT / "bc.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")
print("wrote", [p.name for p in OUT.glob("*.tex")])

# ---- Table 4: centroid-fit residuals
NAME = {"cc": r"$c\bar c$", "bb": r"$b\bar b$"}
rows = [f"{NAME[q['pair']]} ${q['label'][0]}{q['label'][1]}$ & {q['exp']:.1f} & {q['model']:.1f} & ${q['model']-q['exp']:+.1f}$ \\\\" for q in C["residuals"]]
(OUT / "centroid_residuals.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")

# ---- Table 5: centroid-fit stability
rows = []
for k in sorted(C["stability"], key=key_order):
    s = C["stability"][k]
    tag = k.split("__")[1] if "__" in k else "--"
    tag = {"noD": "no $D$", "no3S3P": "no $3S$, $3P$"}.get(tag, tag.replace("kappa", r"$\kappa=") + ("$" if tag.startswith("kappa") else ""))
    rows.append(f"{s['K']} & {s['lam']:g} & {tag} & {s['chi2']:.1f} & {s['prior']:.1f} & "
                + " & ".join(f"{1000*s['sig'][q]:.0f}" for q in ("0.3", "3", "6")) + " & "
                + " & ".join(f"${s['V'][q]:+.3f}$" for q in ("0.3", "3", "6")) + f" & {s['m_c']:.3f} & {s['m_b']:.3f} \\\\")
(OUT / "centroid_stability.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")

# ---- Table 6: B_c centroids from the centroid fit
CL = {"1S": "$1S$", "2S": "$2S$", "1P": "$1P$", "1D": "$1D$", "3S": "$3S$"}
rows = []
for l, (m, e, v) in C["bc"].items():
    ex = C["bc_exp"].get(l)
    rows.append(f"{CL[l]} & ${m:.0f}\\pm{e:.0f}\\pm{v:.0f}$ & " + (f"${ex[0]:.1f}\\pm{ex[1]:.1f}$" if ex else "--") + " \\\\")
s21, e21, v21, x21, ex21, sig = C["s21"]
rows.append(f"$2S-1S$ & ${s21:.1f}\\pm{e21:.1f}\\pm{v21:.1f}$ & ${x21:.1f}\\pm{ex21:.1f}$ \\\\")
(OUT / "centroid_bc.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")
print("wrote", sorted(p.name for p in OUT.glob("*.tex")))
