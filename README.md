# bc-potential-inverse — what does the quarkonium spectrum know about V(r)?

[![ORCID](https://img.shields.io/badge/ORCID-0009--0009--4678--0034-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0009-4678-0034)

Code, fits and manuscript for

> P. A. Oppong, *What does the quarkonium spectrum know about the heavy-quark
> potential? A model-independent extraction of V(r) from charmonium and
> bottomonium, and its consequences for B_c* (2026).

Instead of assuming a Cornell, logarithmic or power-law form, the
quark–antiquark potential is represented by a cubic spline (6–12 knots,
0.15 ≤ r ≤ 8 GeV⁻¹) plus a Coulomb term and fitted, together with m_c, m_b
and a running α_s, to 27 J-resolved cc̄ + bb̄ levels (PDG 2026) with a
smoothness prior.  The output is V(r) **with a 1σ band**, a map of which
level constrains which distance, and B_c predictions whose uncertainty is
derived rather than estimated from a spread of models.

Two extractions are made: from 27 J-resolved levels with first-order
spin-dependent corrections (`fit_spline.py`), and from the 11 spin-averaged
multiplet centroids alone, with no spin-dependent input (`fit_centroid.py`).
All propagated errors are scaled by √(χ²/dof) (1.61 and 1.29).

Headline results (`results/spline_analysis.md`, `results/centroid_analysis.md`):

* The two extractions agree at every r. V(r) is determined (relative to
  V(1 GeV⁻¹)) to 40–70 MeV for 0.7 ≤ r ≤ 7 GeV⁻¹; blind below 0.2 GeV⁻¹.
* χ² = 41 for 27 levels — better than any parametric form (log 59, Martin 70,
  Cornell 134 → Δχ² = 93 for 6 extra parameters).
* Even with a free shape, a flavour-independent static potential fits the
  spin-averaged cc̄ + bb̄ centroids no better than χ²/dof ≈ 2.6 with a 10 MeV
  floor (≈ 17 MeV floor needed): its accuracy is 15–20 MeV, with a robust
  residual pattern (cc̄ 2S–1S 12 MeV too small, bb̄ 16 MeV too large).
* B_c: 1S 6264 ± 7 ± 3, 1³S₁ 6336 ± 8 ± 4, 2S 6874 ± 6 ± 2 MeV (fit ± variants);
  Δ_hyp(1S) = 71 ± 9, Δ_hyp(2S) = 38 ± 6 MeV (ATLAS 64.5, 35.4).
  Spin-averaged centroids 1S 6315 ± 8, 2S 6900 ± 9 vs measured 6322.9, 6897.8.
* The B_c 2S–1S excess that every parametric fit produces (610–616 vs 597)
  is reproduced at 0.8σ in the clean centroid test: it is a property of the
  static-potential description (which misses the cc̄/bb̄ spacings by the same
  amount), not of the cb̄ system.

Companion paper (parametric potentials, same data and spin-dependent
treatment): https://github.com/apoppong24/bc-spectrum-paper.  Solver package
and thesis reproduction: https://github.com/apoppong24/bc-quarkonium.

## Run

```bash
pip install -r requirements.txt
python fit_spline.py 8 1                       # reference fit (~10 min) -> results/spline_K8_lam1.json
python fit_spline.py 8 1 --exclude D           # variants: --exclude D | 3P, --floor 5, --kappa 0.3
for K in 6 10 12; do python fit_spline.py $K 1; done
for lam in 0.1 10 100; do python fit_spline.py 8 $lam; done
python analyse_spline.py                       # band, comparison, sensitivity, B_c -> results/spline_analysis.{md,json}, paper/figures/figA*
python fit_centroid.py 6 1                     # centroid-only fit (~3 min); variants: K 8/10, lam 0.1/10, --kappa 0.4, --exclude D | 3S 3P
python analyse_centroid.py                     # -> results/centroid_analysis.{md,json}, paper/figures/figB*
python paper/make_tables.py                    # -> paper/tables/*.tex  (needs both analysis JSONs)
cd paper && latexmk -pdf main.tex
pytest                                         # solver regression tests
```

All twenty fits used in the paper (12 J-resolved, 8 centroid) are committed under `results/` (with their
logs), so `analyse_spline.py` and the manuscript can be rebuilt without
refitting.

## Layout

```
quarkonium/flexible.py   spline potential, anchor V(1 GeV⁻¹) = −1.166 GeV, prior, J-resolved and centroid fits, band
quarkonium/              Numerov solver, spin-dependent shifts, PDG 2026 data, parametric fits
fit_spline.py            CLI (J-resolved): K lam [--exclude ...] [--floor F] [--kappa k]
fit_centroid.py          CLI (centroids):  K lam [--exclude ...] [--floor F] [--kappa k]
analyse_spline.py        band vs parametric forms, stability table, sensitivity map, B_c with covariance
analyse_centroid.py      centroid band vs J-resolved band, residuals, B_c centroids and the 2S−1S test
paper/                   main.tex, refs.bib, make_tables.py, figures/, tables/
results/                 spline_K*_lam*.json, centroid_K*_lam*.json fits (+ logs), fit_*.json parametric refits, *_analysis.{md,json}
tests/                   solver / spin-structure regression tests
```

## Method in one paragraph

V(r) = −κ/r + S(r), S a natural cubic spline in ln r.  The knot at r = 1 GeV⁻¹
is fixed (an additive constant is degenerate with m_c + m_b).  Residuals are
the 27 level pulls (10 MeV theory floor) followed by √λ × second differences
of the knot values; the Gauss–Newton covariance of that full residual vector
gives the band on V(r) by linear propagation and, through finite-difference
Jacobians, the errors on B_c masses and hyperfine splittings.  K ∈ {6,8,10,12}
and λ ∈ {0.1,1,10,100} are scanned; λ = 100 is over-regularised and excluded
from the systematic.  The centroid fit uses the (2J+1)-weighted multiplet
centroids (three incomplete bb̄ multiplets completed from the pattern of the
lower ones, with an assigned error) and fixes κ, which nothing in the
spin-averaged spectrum constrains.

## Acknowledgement

Extends the author's MSc project at the University of Glasgow (supervisor
Christine Davies).  An AI assistant (Claude, Anthropic) was used to implement
the spline fit and analysis scripts and to help draft the manuscript.

## License

MIT (see `LICENSE`).
