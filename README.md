# bc-potential-inverse — what does the quarkonium spectrum know about V(r)?

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

Headline results (`results/spline_analysis.md`):

* V(r) determined to 30–50 MeV for 0.7 ≤ r ≤ 7 GeV⁻¹; blind below 0.2 GeV⁻¹.
* χ² = 40 for 27 levels — better than any parametric form (log 59, Martin 70).
* Log / power-law potentials lie within the band; Cornell's linear rise is
  excluded at 3σ at r ≈ 3 GeV⁻¹.
* B_c: 1S 6264 ± 4 ± 3, 1³S₁ 6336 ± 5 ± 4, 2S 6874 ± 4 ± 2 MeV (fit ± variants);
  Δ_hyp(1S) = 71 ± 8, Δ_hyp(2S) = 38 ± 5 MeV (ATLAS 64.5, 35.4).
* The B_c 2S–1S spacing comes out 12 ± 7 MeV above experiment and cannot be
  fixed by any change of V(r): the size expected for the flavour-dependent
  1/m correction to a static potential.

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
python analyse_spline.py                       # band, comparison, sensitivity, B_c -> results/spline_analysis.{md,json}, paper/figures/
python paper/make_tables.py                    # -> paper/tables/*.tex
cd paper && latexmk -pdf main.tex
pytest                                         # solver regression tests
```

All twelve fits used in the paper are committed under `results/` (with their
logs), so `analyse_spline.py` and the manuscript can be rebuilt without
refitting.

## Layout

```
quarkonium/flexible.py   spline potential, anchor V(1 GeV⁻¹) = −1.166 GeV, prior, fit, band
quarkonium/              Numerov solver, spin-dependent shifts, PDG 2026 data, parametric fits
fit_spline.py            CLI: K lam [--exclude ...] [--floor F] [--kappa k]
analyse_spline.py        band vs parametric forms, stability table, sensitivity map, B_c with covariance
paper/                   main.tex, refs.bib, make_tables.py, figures/, tables/
results/                 spline_K*_lam*.json fits, fit_*.json parametric refits, spline_analysis.{md,json}
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
from the systematic.

## Acknowledgement

Extends the author's MSc project at the University of Glasgow (supervisor
Christine Davies).  An AI assistant (Claude, Anthropic) was used to implement
the spline fit and analysis scripts and to help draft the manuscript.

## License

MIT (see `LICENSE`).
