# gstar-extremal-effects

Code and data for the numerical example (Section 5) of the paper

> A. Nurwahid, *Disjoint cycles and extremal traces of resolvents of stochastic matrices* (submitted).

Author: Arif Nurwahid, Graduate School of Science and Engineering, Ehime University [![ORCID](https://img.shields.io/badge/ORCID-0009--0005--8802--9629-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0005-8802-9629)

The example uses monthly rainfall in the 38 provinces of Indonesia. It fits a space-time autoregressive model (STAR and GSTAR(1;1)) under two spatial weight matrices, a geographic one ($W_{\mathrm{geo}}$, four nearest provinces) and a correlation-based one ($W_{\mathrm{cor}}$), and computes the average direct and spillover effects of the long-run multiplier $M=(I-\Phi_0-\Phi_1W)^{-1}$. It then computes how far these effects can move when $W$ changes while the estimated coefficients stay fixed, using the bounds of the paper. The script reproduces Table 1 and every number in Section 5.

## Contents

```
data/
  province_rainfall_monthly.csv   monthly rainfall, 38 provinces, 1985-01 to 2025-04 (mm/month)
  provinces.csv                   ISO 3166-2 code, province code, name, centroid lat/lon
  README.md                       sources and processing of the data
scripts/
  pattern_tools.py                helper functions (data, weight matrices, estimation, cycles, integer programs)
  example_provinces.py            Section 5 (Table 1 and text)
results/
  example_provinces.txt           output of example_provinces.py --all-cycles
  example_provinces.json          all numbers of Section 5 at full precision
```

## Requirements

Python 3.10 or later and

```
pip install -r requirements.txt
```

The integer programs use the HiGHS solver through `scipy.optimize.milp`, so no separate solver is needed. The results in `results/` were produced with Python 3.13, NumPy 2.5, SciPy 1.18, pandas 3.0 and NetworkX 3.7.

## How to run

From the repository root:

```
python scripts/example_provinces.py
```

This takes under one minute on a laptop and writes `results/example_provinces.json`. The option `--all-cycles` also enumerates all 576,397 simple cycles of the geographic digraph (a few seconds more). The printed output should match `results/example_provinces.txt`.

## What the script computes

1. Standardised monthly anomalies of each series (by calendar month).
2. The weight matrices $W_{\mathrm{geo}}$ and $W_{\mathrm{cor}}$, and the structure of their digraphs: strong connectivity, the longest cycle (by integer programming), and the largest number of disjoint 2-cycles.
3. Pooled least-squares estimates of the homogeneous STAR model and the resulting average direct and spillover effects (Table 1).
4. The smallest and largest average effects over all weight matrices with the same pattern, and over all weight matrices, with the coefficients held fixed.
5. Location-wise least-squares estimates of the GSTAR(1;1) model under $W_{\mathrm{geo}}$, and the range of the average direct effect over all weight matrices with the pattern of $W_{\mathrm{geo}}$.

Where an extreme is found by an integer program over cycles up to a fixed length, the script also reports a rigorous bound on what longer cycles could add.

## Data

The rainfall series are areal means of CHIRPS v2.0 (Funk et al., 2015) over province polygons. See `data/README.md` for sources, processing and terms of use.

## Citation

If you use this code or data, please cite the version archived at Zenodo together with the paper (details will be added after publication) and the data sources listed in `data/README.md`. A `CITATION.cff` file is included.

## License

The code is released under the MIT License (see `LICENSE`). Terms for the data files are given in `data/README.md`.
