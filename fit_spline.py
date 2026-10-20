#!/usr/bin/env python3
"""Fit the spline potential:  python fit_spline.py K lam [--exclude D 3P] [--floor 10]
Writes results/spline_K<K>_lam<lam>[tag].json"""
import argparse
from pathlib import Path
from quarkonium.flexible import fit_spline

ap = argparse.ArgumentParser()
ap.add_argument("K", type=int); ap.add_argument("lam", type=float)
ap.add_argument("--exclude", nargs="*", default=[]); ap.add_argument("--floor", type=float, default=10.0); ap.add_argument("--kappa", type=float, default=None)
a = ap.parse_args()
tag = ("__no" + "".join(a.exclude) if a.exclude else "") + (f"__floor{a.floor:g}" if a.floor != 10.0 else "") + (f"__kappa{a.kappa:g}" if a.kappa is not None else "")
r = fit_spline(K=a.K, lam=a.lam, sigma_th=a.floor, exclude=tuple(a.exclude), verbose=False, kappa_fixed=a.kappa)
out = Path("results") / f"spline_K{a.K}_lam{a.lam:g}{tag}.json"
r.to_json(out)
print(f"K={a.K} lam={a.lam:g}{tag}: chi2_data={r.chi2_data:.1f} chi2_prior={r.chi2_prior:.1f} dof={r.dof} ({r.seconds:.0f}s)")
for n, x, e in zip(r.names, r.x, r.cov.diagonal() ** 0.5):
    print(f"   {n:10s} = {x:8.4f} ± {e:.4f}")
