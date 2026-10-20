#!/usr/bin/env python3
"""Fit the spline potential to spin-averaged centroids only:
    python fit_centroid.py K lam [--kappa 0] [--floor 10] [--exclude D 3S 3P]
Writes results/centroid_K<K>_lam<lam>[tag].json"""
import argparse
from pathlib import Path
from quarkonium.flexible import fit_spline_centroid

ap = argparse.ArgumentParser()
ap.add_argument("K", type=int); ap.add_argument("lam", type=float)
ap.add_argument("--exclude", nargs="*", default=[]); ap.add_argument("--floor", type=float, default=10.0); ap.add_argument("--kappa", type=float, default=0.0)
a = ap.parse_args()
tag = ("__no" + "".join(a.exclude) if a.exclude else "") + (f"__floor{a.floor:g}" if a.floor != 10.0 else "") + (f"__kappa{a.kappa:g}" if a.kappa != 0.0 else "")
r = fit_spline_centroid(K=a.K, lam=a.lam, sigma_th=a.floor, kappa=a.kappa, exclude=tuple(a.exclude), verbose=False)
out = Path("results") / f"centroid_K{a.K}_lam{a.lam:g}{tag}.json"
r.to_json(out)
print(f"K={a.K} lam={a.lam:g}{tag}: chi2_data={r.chi2_data:.1f} chi2_prior={r.chi2_prior:.1f} dof={r.dof} ({r.seconds:.0f}s)")
for n, x, e in zip(r.names, r.x, r.cov.diagonal() ** 0.5):
    print(f"   {n:10s} = {x:8.4f} ± {e:.4f}")
for p, l, m, mod, pull in r.residuals:
    print(f"   {p} {l}: exp {1000*m:.1f} model {1000*mod:.1f} pull {pull:+.2f}")
