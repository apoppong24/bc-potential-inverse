#!/usr/bin/env python3
"""Analysis of the spline-potential fits.

  1. V(r) with 1σ band for the reference fit (K, λ chosen), overlaid with the
     refitted Martin / log / Cornell potentials from the B_c paper.
  2. Stability: bands for all available (K, λ) fits; χ²_data vs λ.
  3. Sensitivity map: dE_i/ds_k for every fitted level (which r each level probes).
  4. Model-independent B_c predictions with propagated band.
  5. Where is V(r) determined?  σ_V(r) and the r-range where σ_V < 20 MeV.

Writes results/spline_analysis.md/.json and paper/figures/*.pdf
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np

from quarkonium import pdg
from quarkonium.fit import FAMILIES, predict
from quarkonium.flexible import SplineFit, make_spline_potential
from quarkonium.uncertainty import predict_with_errors as _pwe

RES = Path("results"); FIG = Path("paper/figures"); FIG.mkdir(parents=True, exist_ok=True)
BC_LAB = ["11S0", "13S1", "21S0", "23S1", "13P0", "1P1", "1P1'", "13P2", "13D1", "13D3", "33S1"]


def load_all():
    fits = {}
    for p in sorted(glob.glob(str(RES / "spline_K*_lam*.json"))):
        f = SplineFit.from_json(p)
        fits[Path(p).stem.replace("spline_", "")] = f
    return fits


def bc_with_band(f: SplineFit, labels=BC_LAB, h=1e-3, rel_step=2e-3):
    """Linear propagation of the full spline-fit covariance to B_c masses."""
    K = f.K
    x = f.x.copy()
    def pred(xx):
        P = make_spline_potential(f.knots, xx[:K - 1], *xx[K - 1:])
        return predict(P, "cb", labels, h=h)
    c = pred(x)
    J = np.zeros((len(labels), len(x)))
    for j in range(len(x)):
        dx = rel_step * max(abs(x[j]), 1e-2)
        xp, xm = x.copy(), x.copy(); xp[j] += dx; xm[j] -= dx
        pp, pm = pred(xp), pred(xm)
        J[:, j] = [(pp[l] - pm[l]) / (2 * dx) for l in labels]
    var = np.einsum("ij,jk,ik->i", J, f.cov, J)
    return {l: (1000 * c[l], 1000 * np.sqrt(max(v, 0))) for l, v in zip(labels, var)}, J


def main(ref_key="K8_lam1"):
    fits = load_all()
    print("fits:", list(fits))
    ref = fits[ref_key]
    r = np.geomspace(0.1, 12.0, 300)
    md = [f"# Model-independent potential from the cc̄+bb̄ spectrum\n", f"Reference fit: {ref_key} "
          f"(χ²_data = {ref.chi2_data:.1f}, χ²_prior = {ref.chi2_prior:.1f}, dof = {ref.dof}).\n"]

    # ---- 1. band + parametric comparison
    V, sV = ref.band(r)
    par = {}
    for fam, lab in (("martin_run_geometric", "Martin (running)"), ("log_run_geometric", "Log (running)"), ("cornell", "Cornell")):
        d = json.loads((RES / f"fit_{fam}.json").read_text())
        Pp = FAMILIES[fam]["make"](*d["x"])
        par[lab] = Pp(r) - Pp(np.array([1.0]))[0] + (-1.166)      # shift to V(1) = V0: compare shapes
    md.append("## Fitted parameters\n")
    for n, x, e in zip(ref.names, ref.x, np.sqrt(ref.cov.diagonal())):
        md.append(f"* {n} = {x:.4f} ± {e:.4f}")
    md.append("")
    # where is V determined
    good = sV < 0.020
    if good.any():
        md.append(f"## Where V(r) is determined\n\nσ_V < 20 MeV for {r[good].min():.2f} ≤ r ≤ {r[good].max():.2f} GeV⁻¹; "
                  f"σ_V < 50 MeV for {r[sV<0.05].min():.2f} ≤ r ≤ {r[sV<0.05].max():.2f} GeV⁻¹.\n")
    md.append("Parametric potentials are shifted by a constant to V(1) = V₀ (the constant is absorbed by the quark masses), so this compares shapes.\n")
    md.append("| r [GeV⁻¹] | V ± σ [GeV] | Martin | Log | Cornell |")
    md.append("|---|---|---|---|---|")
    for rr in (0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        i = int(np.argmin(np.abs(r - rr)))
        md.append(f"| {rr} | {V[i]:+.3f} ± {sV[i]:.3f} | " + " | ".join(f"{par[k][i]:+.3f}" for k in par) + " |")
    md.append("")

    # ---- 2. stability table
    md.append("## Stability over K and λ\n")
    md.append("| fit | K | λ | χ²_data | χ²_prior | σ_V(0.3) | σ_V(1) | σ_V(3) | σ_V(6) | V(0.3) | V(3) | V(6) |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    stab = {}
    for k, f in fits.items():
        Vk, sk = f.band(r)
        pick = lambda rr: int(np.argmin(np.abs(r - rr)))
        stab[k] = dict(K=f.K, lam=f.lam, chi2=f.chi2_data, prior=f.chi2_prior,
                       sig={rr: float(sk[pick(rr)]) for rr in (0.3, 1, 3, 6)}, V={rr: float(Vk[pick(rr)]) for rr in (0.3, 3, 6)})
        md.append(f"| {k} | {f.K} | {f.lam:g} | {f.chi2_data:.1f} | {f.chi2_prior:.1f} | "
                  + " | ".join(f"{sk[pick(rr)]*1000:.0f}" for rr in (0.3, 1, 3, 6)) + " | "
                  + " | ".join(f"{Vk[pick(rr)]:+.3f}" for rr in (0.3, 3, 6)) + " |")
    md.append("\n(σ_V in MeV, V in GeV)\n")

    # ---- 3. sensitivity map from the fit Jacobian: refit-free estimate via finite differences on knots
    K = ref.K
    labels = [(p, l) for p, l, *_ in ref.residuals]
    base = {}
    P0 = ref.potential()
    for pair in ("cc", "bb"):
        labs = [l for p, l in labels if p == pair]
        base[pair] = predict(P0, pair, labs, h=2e-3)
    S = np.zeros((len(labels), K - 1))
    for k in range(K - 1):
        xx = ref.x.copy(); dx = 0.02; xx[k] += dx
        Pk = make_spline_potential(ref.knots, xx[:K - 1], *xx[K - 1:])
        for pair in ("cc", "bb"):
            labs = [l for p, l in labels if p == pair]
            pk = predict(Pk, pair, labs, h=2e-3)
            for i, (p, l) in enumerate(labels):
                if p == pair:
                    S[i, k] = (pk[l] - base[pair][l]) / dx
    k1 = int(np.argmin(np.abs(np.log(ref.knots))))
    free_knots = np.delete(ref.knots, k1)
    md.append("## Sensitivity dE/ds_k (dimensionless): which knot each level feels\n")
    md.append("| level | " + " | ".join(f"r={rk:.2f}" for rk in free_knots) + " | peak r |")
    md.append("|---|" + "---|" * (K - 1) + "---|")
    sens = {}
    for i, (p, l) in enumerate(labels):
        peak = free_knots[int(np.argmax(np.abs(S[i])))]
        sens[f"{p} {l}"] = dict(row=S[i].tolist(), peak=float(peak))
        md.append(f"| {p} {l} | " + " | ".join(f"{v:+.2f}" for v in S[i]) + f" | {peak:.2f} |")
    md.append("")

    # ---- 4. B_c predictions with band
    md.append("## B_c predictions (MeV) with the model-independent band\n")
    bc, Jbc = bc_with_band(ref)
    md.append("| state | mass ± σ | PDG / ATLAS |")
    md.append("|---|---|---|")
    exp = {"11S0": pdg.BC["11S0"], "21S0": pdg.BC["21S0"], "13S1": pdg.BC["13S1"]}
    for l in BC_LAB:
        e = exp.get(l)
        md.append(f"| {l} | {bc[l][0]:.0f} ± {bc[l][1]:.0f} | " + (f"{e[0]:.1f} ± {e[1]}" if e else "–") + " |")
    d1 = bc["13S1"][0] - bc["11S0"][0]; d2 = bc["23S1"][0] - bc["21S0"][0]
    # errors on differences need the covariance; do it properly
    idx = {l: i for i, l in enumerate(BC_LAB)}
    cv = Jbc @ ref.cov @ Jbc.T * 1e6
    def err_diff(a, b):
        return np.sqrt(max(cv[idx[a], idx[a]] + cv[idx[b], idx[b]] - 2 * cv[idx[a], idx[b]], 0))
    md.append(f"\nΔ_1S = {d1:.1f} ± {err_diff('13S1','11S0'):.1f} MeV (ATLAS 64.5 ± 1.7);  Δ_2S = {d2:.1f} ± {err_diff('23S1','21S0'):.1f} (35.4 ± 2.4);  "
              f"Δ_1S−Δ_2S = {d1-d2:.1f} (CMS 29.1 ± 1.7)\n")

    # ---- 5. spread of B_c predictions across all spline variants (λ, K, data, κ, floor)
    md.append("## B_c predictions across spline variants (MeV) — the model-independent systematic\n")
    md.append("| fit | " + " | ".join(BC_LAB) + " | Δ1S | Δ2S |")
    md.append("|---|" + "---|" * (len(BC_LAB) + 2))
    allbc = {}
    for k, f in sorted(fits.items()):
        P = f.potential()
        pr = {l: 1000 * v for l, v in predict(P, "cb", BC_LAB, h=1e-3).items()}
        allbc[k] = pr
        md.append(f"| {k} | " + " | ".join(f"{pr[l]:.0f}" for l in BC_LAB)
                  + f" | {pr['13S1']-pr['11S0']:.1f} | {pr['23S1']-pr['21S0']:.1f} |")
    core = [k for k in allbc if "lam100" not in k]          # λ=100 is over-regularised (χ²_prior ≫ 1)
    md.append("| **half-range (excl. λ=100)** | " + " | ".join(f"{(max(allbc[k][l] for k in core)-min(allbc[k][l] for k in core))/2:.0f}" for l in BC_LAB)
              + f" | {(max(allbc[k]['13S1']-allbc[k]['11S0'] for k in core)-min(allbc[k]['13S1']-allbc[k]['11S0'] for k in core))/2:.1f}"
              + f" | {(max(allbc[k]['23S1']-allbc[k]['21S0'] for k in core)-min(allbc[k]['23S1']-allbc[k]['21S0'] for k in core))/2:.1f} |")
    md.append("")

    text = "\n".join(md)
    (RES / "spline_analysis.md").write_text(text)
    (RES / "spline_analysis.json").write_text(json.dumps(dict(
        ref=ref_key, r=r.tolist(), V=V.tolist(), sigma=sV.tolist(), parametric={k: v.tolist() for k, v in par.items()},
        stability=stab, sensitivity=sens, knots=free_knots.tolist(),
        bc={l: list(v) for l, v in bc.items()},
        hyp=dict(d1S=[d1, err_diff('13S1', '11S0')], d2S=[d2, err_diff('23S1', '21S0')],
                 s21=[bc["21S0"][0] - bc["11S0"][0], err_diff('21S0', '11S0')]),
        bc_variants=allbc), indent=1))
    print(text)

    # ---- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm", "legend.fontsize": 7.5})

    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.fill_between(r, V - sV, V + sV, color="C0", alpha=0.25, label=r"spline $\pm1\sigma$")
    ax.fill_between(r, V - 2 * sV, V + 2 * sV, color="C0", alpha=0.10)
    ax.plot(r, V, color="C0", lw=1.2)
    for (lab, Vp), ls in zip(par.items(), ("--", "-.", ":")):
        ax.plot(r, Vp, ls, color="k", lw=0.9, label=lab)
    for rk in ref.knots:
        ax.axvline(rk, color="0.8", lw=0.5, zorder=0)
    ax.set_xscale("log"); ax.set_xlim(0.1, 12); ax.set_ylim(-3.5, 2.0)
    ax.set_xlabel(r"$r$ [GeV$^{-1}$]"); ax.set_ylabel(r"$V(r)$ [GeV]"); ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "figA1_potential_band.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 3.4))
    for k, f in sorted(fits.items(), key=lambda kv: (kv[1].K, kv[1].lam, kv[0])):
        if "__" in k and not k.endswith("noD"):
            continue                                        # floor / kappa / no3P variants: in the table only
        Vk, sk = f.band(r)
        ax.plot(r, 1000 * sk, lw=1, label=f"K={f.K}, λ={f.lam:g}" + (" (no D)" if "noD" in k else ""))
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 12); ax.set_ylim(2, 2000)
    ax.set_xlabel(r"$r$ [GeV$^{-1}$]"); ax.set_ylabel(r"$\sigma_V$ [MeV]")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False, columnspacing=0.8)
    fig.tight_layout(); fig.savefig(FIG / "figA2_band_width.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    im = ax.imshow(np.abs(S), aspect="auto", cmap="Blues", origin="lower")
    ax.set_xticks(range(K - 1)); ax.set_xticklabels([f"{rk:.2f}" for rk in free_knots])
    ax.set_yticks(range(len(labels))); ax.set_yticklabels([f"{p} {l}" for p, l in labels], fontsize=6)
    ax.set_xlabel(r"knot position $r_k$ [GeV$^{-1}$]"); plt.colorbar(im, ax=ax, label=r"$|dE/ds_k|$")
    fig.tight_layout(); fig.savefig(FIG / "figA3_sensitivity.pdf"); plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "K8_lam1")
