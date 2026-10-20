#!/usr/bin/env python3
"""Analysis of the spin-averaged (centroid) spline fits.

  1. V(r) with band for the reference centroid fit, vs the J-resolved fit and
     the shifted parametric potentials.
  2. Stability over K, λ, κ, floor, data subset.
  3. B_c spin-averaged centroids with propagated covariance, and the clean
     2S−1S test against the measured B_c centroids.
  4. Residuals of the fitted cc̄ / bb̄ centroids (what the model gets wrong).

Writes results/centroid_analysis.md/.json and paper/figures/figB*.pdf
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np

from quarkonium.fit import FAMILIES
from quarkonium.flexible import SplineFit, make_spline_potential, predict_centroids

RES = Path("results"); FIG = Path("paper/figures"); FIG.mkdir(parents=True, exist_ok=True)
BC_LAB = ["1S", "2S", "1P", "1D", "3S"]
# measured B_c centroids: 1S from B_c(6274.5±0.3) and B_c*(ATLAS 6339.0±1.7);
# 2S from B_c(2S)(6871.2±1.0) + (3/4)·Δ_hyp(2S) (ATLAS+CMS 35.4±2.4)
BC_EXP = {"1S": ((6274.5 + 3 * 6339.0) / 4, np.sqrt(0.3**2 + (3 * 1.8)**2) / 4),      # ATLAS 64.5 ± 1.4 +1.0/−1.4 symmetrised
          "2S": (6871.2 + 0.75 * 35.4, np.hypot(1.0, 0.75 * 2.5))}


def load_all():
    return {Path(p).stem.replace("centroid_", ""): SplineFit.from_json(p)
            for p in sorted(glob.glob(str(RES / "centroid_K*_lam*.json")))}


def bc_with_band(f: SplineFit, labels=BC_LAB, h=1e-3, rel_step=2e-3):
    K = f.K
    x = f.x.copy()

    def pred(xx):
        return predict_centroids(make_spline_potential(f.knots, xx[:K - 1], *xx[K - 1:]), "cb", labels, h=h)
    c = pred(x)
    J = np.zeros((len(labels), len(x)))
    for j in range(len(x)):
        if np.sqrt(f.cov[j, j]) < 1e-6:          # fixed parameter
            continue
        dx = rel_step * max(abs(x[j]), 1e-2)
        xp, xm = x.copy(), x.copy(); xp[j] += dx; xm[j] -= dx
        pp, pm = pred(xp), pred(xm)
        J[:, j] = [(pp[l] - pm[l]) / (2 * dx) for l in labels]
    cv = J @ f.cov @ J.T * 1e6
    return {l: 1000 * c[l] for l in labels}, cv


def main(ref_key="K6_lam1"):
    fits = load_all()
    print("fits:", list(fits))
    ref = fits[ref_key]
    def scale_of(f):
        """PDG-style error scaling √(χ²/dof) when χ²/dof > 1; dof counts the prior
        residuals as (K−2) soft data points so that it stays positive for K=10."""
        dof = max(f.dof + f.K - 2, 1)
        return np.sqrt(max((f.chi2_data + f.chi2_prior) / dof, 1.0))
    scale = scale_of(ref)
    r = np.geomspace(0.1, 12.0, 300)
    md = [f"# Model-independent potential from the spin-averaged cc̄+bb̄ spectrum\n",
          f"Reference fit: {ref_key} (χ²_data = {ref.chi2_data:.1f}, χ²_prior = {ref.chi2_prior:.1f}, dof = {ref.dof}, "
          f"error scale factor √(χ²/dof) = {scale:.2f} applied to all propagated errors below).\n"]

    V, sV = ref.band(r)
    sV = sV * scale
    par = {}
    for fam, lab in (("martin_run_geometric", "Martin"), ("log_run_geometric", "Log"), ("cornell", "Cornell")):
        d = json.loads((RES / f"fit_{fam}.json").read_text())
        Pp = FAMILIES[fam]["make"](*d["x"])
        par[lab] = Pp(r) - Pp(np.array([1.0]))[0] + (-1.166)
    jres = SplineFit.from_json(RES / "spline_K8_lam1.json")
    Vj, sVj = jres.band(r)
    md.append("## Fitted parameters\n")
    for n, x, e in zip(ref.names, ref.x, np.sqrt(ref.cov.diagonal()) * scale):
        md.append(f"* {n} = {x:.4f} ± {e:.4f}")
    md.append("")
    md.append("## Residuals of the fitted centroids (MeV)\n\n| level | exp | model | model − exp | pull |\n|---|---|---|---|---|")
    for p, l, m, mod, pull in ref.residuals:
        md.append(f"| {p} {l} | {1000*m:.1f} | {1000*mod:.1f} | {1000*(mod-m):+.1f} | {pull:+.2f} |")
    md.append("")
    md.append(f"σ_V < 50 MeV for {r[sV<0.05].min():.2f} ≤ r ≤ {r[sV<0.05].max():.2f} GeV⁻¹; < 100 MeV for {r[sV<0.1].min():.2f} ≤ r ≤ {r[sV<0.1].max():.2f}.\n")
    md.append("| r | V ± σ (centroid fit) | J-resolved fit | Martin | Log | Cornell |\n|---|---|---|---|---|---|")
    pick = lambda rr: int(np.argmin(np.abs(r - rr)))
    for rr in (0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        i = pick(rr)
        md.append(f"| {rr} | {V[i]:+.3f} ± {sV[i]:.3f} | {Vj[i]:+.3f} ± {sVj[i]:.3f} | " + " | ".join(f"{par[k][i]:+.3f}" for k in par) + " |")
    md.append("")

    md.append("## Stability\n\n| fit | K | λ | χ²_data | χ²_prior | dof | σ_V(0.3) | σ_V(3) | σ_V(6) | V(0.3) | V(3) | V(6) | m_c | m_b |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    stab = {}
    for k, f in fits.items():
        Vk, sk = f.band(r); sk = sk * scale_of(f)
        stab[k] = dict(K=f.K, lam=f.lam, chi2=f.chi2_data, prior=f.chi2_prior, dof=f.dof,
                       sig={rr: float(sk[pick(rr)]) for rr in (0.3, 3, 6)}, V={rr: float(Vk[pick(rr)]) for rr in (0.3, 3, 6)},
                       m_c=float(f.x[f.K]), m_b=float(f.x[f.K + 1]))
        md.append(f"| {k} | {f.K} | {f.lam:g} | {f.chi2_data:.1f} | {f.chi2_prior:.1f} | {f.dof} | "
                  + " | ".join(f"{sk[pick(rr)]*1000:.0f}" for rr in (0.3, 3, 6)) + " | "
                  + " | ".join(f"{Vk[pick(rr)]:+.3f}" for rr in (0.3, 3, 6)) + f" | {f.x[f.K]:.3f} | {f.x[f.K+1]:.3f} |")
    md.append("")

    md.append("## B_c spin-averaged centroids (MeV)\n")
    bc, cv = bc_with_band(ref)
    cv = cv * scale**2
    idx = {l: i for i, l in enumerate(BC_LAB)}
    err = {l: float(np.sqrt(cv[idx[l], idx[l]])) for l in BC_LAB}
    allbc = {}
    for k, f in sorted(fits.items()):
        P = f.potential()
        allbc[k] = {l: 1000 * v for l, v in predict_centroids(P, "cb", BC_LAB, h=1e-3).items()}
    core = [k for k in allbc if "lam100" not in k]
    half = lambda g: (max(g(allbc[k]) for k in core) - min(g(allbc[k]) for k in core)) / 2
    md.append("| level | prediction ± fit ± variants | experiment |\n|---|---|---|")
    for l in BC_LAB:
        e = BC_EXP.get(l)
        md.append(f"| {l} | {bc[l]:.0f} ± {err[l]:.0f} ± {half(lambda v: v[l]):.0f} | " + (f"{e[0]:.1f} ± {e[1]:.1f}" if e else "–") + " |")
    s21 = bc["2S"] - bc["1S"]
    e21 = float(np.sqrt(cv[idx['2S'], idx['2S']] + cv[idx['1S'], idx['1S']] - 2 * cv[idx['2S'], idx['1S']]))
    v21 = half(lambda v: v["2S"] - v["1S"])
    x21 = BC_EXP["2S"][0] - BC_EXP["1S"][0]; ex21 = np.hypot(BC_EXP["2S"][1], BC_EXP["1S"][1])
    md.append(f"| 2S − 1S | {s21:.1f} ± {e21:.1f} ± {v21:.1f} | {x21:.1f} ± {ex21:.1f} |")
    sig = (s21 - x21) / np.sqrt(e21**2 + v21**2 + ex21**2)
    md.append(f"\n2S−1S excess: {s21-x21:+.1f} MeV, {sig:.1f}σ (fit ⊕ variants ⊕ experiment).\n")
    # the same spacing in the fitted systems
    res = {(p, l): (m, mod) for p, l, m, mod, _ in ref.residuals}
    for pair in ("cc", "bb"):
        m1, mo1 = res[(pair, "1S")]; m2, mo2 = res[(pair, "2S")]
        md.append(f"* {pair} 2S−1S: exp {1000*(m2-m1):.1f}, model {1000*(mo2-mo1):.1f} ({1000*((mo2-mo1)-(m2-m1)):+.1f})")
    md.append("")
    md.append("## B_c centroids across variants\n\n| fit | " + " | ".join(BC_LAB) + " | 2S−1S |\n|---|" + "---|" * (len(BC_LAB) + 1))
    for k in sorted(allbc):
        md.append(f"| {k} | " + " | ".join(f"{allbc[k][l]:.0f}" for l in BC_LAB) + f" | {allbc[k]['2S']-allbc[k]['1S']:.1f} |")

    text = "\n".join(md)
    (RES / "centroid_analysis.md").write_text(text)
    (RES / "centroid_analysis.json").write_text(json.dumps(dict(
        ref=ref_key, scale=float(scale), r=r.tolist(), V=V.tolist(), sigma=sV.tolist(), Vj=Vj.tolist(), sigmaj=sVj.tolist(),
        parametric={k: v.tolist() for k, v in par.items()}, stability=stab,
        residuals=[dict(pair=p, label=l, exp=1000 * m, model=1000 * mod, pull=pull) for p, l, m, mod, pull in ref.residuals],
        bc={l: [bc[l], err[l], half(lambda v, l=l: v[l])] for l in BC_LAB},
        s21=[s21, e21, v21, x21, ex21, float(sig)], bc_exp=BC_EXP, bc_variants=allbc), indent=1))
    print(text)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm", "legend.fontsize": 7.5})
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.fill_between(r, V - sV, V + sV, color="C0", alpha=0.3, label=r"centroid fit $\pm1\sigma$")
    ax.plot(r, V, color="C0", lw=1.2)
    ax.fill_between(r, Vj - sVj, Vj + sVj, color="C3", alpha=0.15, label=r"$J$-resolved fit $\pm1\sigma$")
    ax.plot(r, Vj, color="C3", lw=0.8)
    for (lab, Vp), ls in zip(par.items(), ("--", "-.", ":")):
        ax.plot(r, Vp, ls, color="k", lw=0.9, label=lab)
    for rk in ref.knots:
        ax.axvline(rk, color="0.8", lw=0.5, zorder=0)
    ax.set_xscale("log"); ax.set_xlim(0.1, 12); ax.set_ylim(-3.5, 2.0)
    ax.set_xlabel(r"$r$ [GeV$^{-1}$]"); ax.set_ylabel(r"$V(r)$ [GeV]"); ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "figB1_centroid_band.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 3.4))
    for k, f in sorted(fits.items(), key=lambda kv: (kv[1].K, kv[1].lam, kv[0])):
        Vk, sk = f.band(r); sk = sk * scale_of(f)
        tag = k.split("__")[1] if "__" in k else ""
        ax.plot(r, 1000 * sk, lw=1, label=f"K={f.K}, λ={f.lam:g}" + (f" ({tag})" if tag else ""))
    ax.plot(r, 1000 * sVj, "k--", lw=1, label="J-resolved K=8, λ=1")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 12); ax.set_ylim(2, 2000)
    ax.set_xlabel(r"$r$ [GeV$^{-1}$]"); ax.set_ylabel(r"$\sigma_V$ [MeV]")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False, columnspacing=0.8)
    fig.tight_layout(); fig.savefig(FIG / "figB2_centroid_width.pdf"); plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "K6_lam1")
