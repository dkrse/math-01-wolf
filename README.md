# Empirical subleading correction to Wolf's formula for consecutive prime gaps

Reproducibility repository for the paper *"An empirical subleading correction
to Wolf's formula for consecutive prime gaps"* (Sestak, 2026). Contains the
C generator, Python analysis pipeline, 11 reproducible Jupyter notebooks,
generated figures and tables, and the LaTeX source of the manuscript.

## Headline result

On a dataset of *n* = 140 points covering 25 values of *N* ∈ \[10⁵, 10¹³\]
(primary filter ρ = g·C_g/ln N ∈ [0.05, 1.10]) the residual of Wolf's
formula `W(g, N) = C_g · Li_2(N) · exp(-g·C_g/ln N)` is well described by
the two-parameter, no-intercept model M₁\*

```
R(g, N) = log( N_g_emp(N) / W(g, N) ) ≈ 0.855·√ω(g) − 0.202·log log N
```

with median hold-out-one-*N* `R²_CV = 0.92` (plateau ≈ 0.94 for
*N* ≥ 10⁹, 0.874 at the held-out *N* = 10¹³) and `ΔAIC ≈ −39` against the
reparametrised Wolf alternative M₂. A five-parameter joint model M₃
(reparametrised Wolf + linear residual) fits better still
(`ΔAIC = −110.8` vs M₁\*); M₁\* is kept as the minimally parametrised
descriptor, not the best fit. Under the primary filter ω(g) ∈ {1, 2}, so
√ω, log ω, ω, ω−1, √(ω−1) are statistically indistinguishable; a
relaxed-filter robustness check (ρ_max = 4, n = 553, ω ∈ {1, 2, 3},
outside Wolf's regime of validity) prefers a linear ω.

## Repository layout

```
.
├── README.md                   ← this file
├── .gitignore
├── requirements.txt            ← pinned Python dependencies
│
├── app-c/                      ← C generator (primesieve-based)
│   ├── main-csv.c              data + histogram producer
│   ├── main.c                  minimal demo
│   ├── Makefile
│   └── README.md
│
├── app-python/                 ← descriptive analysis (upstream)
│   ├── analyze_gaps.py         histograms, Hardy-Littlewood, Wolf, χ²
│   ├── analyze_sequence.py     champions, GPY ratio, autocorrelation
│   └── README.md
│
├── 01..11_*.ipynb              ← 11 reproducible notebooks (paper layer)
│   ├── 01_validate_histogram      C-sieve vs numpy at N = 10⁸
│   ├── 02_validate_theory         independent C_g, Li_2, pair counts
│   ├── 03_validate_sequence       champions, GPY, autocorrelation
│   ├── 04_fixed_theory            corrected C_g, χ² recompute
│   ├── 05_ml_bug_impact           bug impact on residual model
│   ├── 06_bootstrap_ci            B = 1000 bootstrap on audit range
│   ├── 07_extended_to_1e12        extension to N = 10¹²
│   ├── 08_missing_figures         remaining paper figures
│   ├── 09_rev01_alt_forms         alternative ω-features + LR test on c
│   ├── 10_rev01_robustness        filter / weight / block-bootstrap scan
│   └── 11_rev01_m1_vs_m3          formal M₁* vs M₃ comparison
│
├── ml_data/                    ← cached gap histograms (CSV per N)
│   └── gaps_N{1e5 ... 1e13}.csv      25 values of N (gitignored above 1e10)
│
├── main-csv-static             ← compiled C binary used by notebooks
│
├── outputs/                    ← generated artefacts
│
├── paper/                      ← manuscript
│
└── docs/                       ← documentation, audit trail, project notes
```

## Quickstart

### Build the C generator

```bash
cd app-c
make                            # builds dynamic + static binaries into build/
cp build/main-csv-static ../    # used by notebooks
```

Requires `libprimesieve` headers (`apt install libprimesieve-dev` on Debian /
Ubuntu).

### Set up Python

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Run the notebooks

```bash
jupyter lab                     # then open 01_validate_histogram.ipynb
```

Notebooks 01 → 11 may be run in order. Each is self-contained, runs end-to-
end, and uses caches in `ml_data/` for histograms above N = 10⁹. The full
chain takes ≈30 minutes on a modern laptop (sieves up to N = 10¹³ included;
the 10¹² and 10¹³ caches are reused if present in `ml_data/`).

## Data flow

```
   primesieve  ──>  C generator (app-c)  ──>  ml_data/*.csv
                                                  │
                       ┌──────────────────────────┤
                       ▼                          ▼
         analyze_gaps.py            01..11_*.ipynb (model selection,
         analyze_sequence.py         AIC, bootstrap, CV, M₁* vs M₃)
                       │                          │
                       └──────────┬───────────────┘
                                  ▼
                outputs/{figures, results, tables}
                                  │
                                  ▼
                     paper/paper-en.{tex,pdf}
```

## Citation

If you use this code, please cite

> Sestak, K. (2026). *An empirical subleading correction to Wolf's formula
> for consecutive prime gaps.* Manuscript / preprint.

Author: Kristian Sestak (kristian.sestak@gmail.com).
Repository: https://github.com/dkrse/math-01-wolf
