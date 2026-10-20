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


def at(arr, x):
    return float(np.interp(x, r, arr))


# ---- Table 1: V(r) ± σ vs shifted parametric potentials
rows = []
for x in (0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
    rows.append(f"{x:g} & ${at(V, x):+.3f}\\pm{at(sV, x):.3f}$ & "
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
    rows.append(f"{s['K']} & {s['lam']:g} & {tag} & {s['chi2']:.1f} & {s['prior']:.1f} & "
                + " & ".join(f"{1000*s['sig'][q]:.0f}" for q in ("0.3", "1", "3", "6")) + " & "
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
    rows.append(f"{LAB[l]} & ${m:.0f}\\pm{e:.0f}\\pm{half(lambda v: v[l]):.0f}$ & "
                + (f"${ex[0]:.1f}\\pm{ex[1]}$" if ex else "--") + r" \\")
h = A["hyp"]
rows.append(f"$\\Delta_{{1S}}$ & ${h['d1S'][0]:.1f}\\pm{h['d1S'][1]:.1f}\\pm{half(lambda v: v['13S1']-v['11S0']):.1f}$ & $64.5\\pm1.7$ \\\\")
rows.append(f"$\\Delta_{{2S}}$ & ${h['d2S'][0]:.1f}\\pm{h['d2S'][1]:.1f}\\pm{half(lambda v: v['23S1']-v['21S0']):.1f}$ & $35.4\\pm2.4$ \\\\")
rows.append(f"$M(2\\,^1S_0)-M(1\\,^1S_0)$ & ${h['s21'][0]:.0f}\\pm{h['s21'][1]:.0f}\\pm{half(lambda v: v['21S0']-v['11S0']):.0f}$ & $596.7\\pm1.0$ \\\\")
(OUT / "bc.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")
print("wrote", [p.name for p in OUT.glob("*.tex")])
