# HIF and PACAP gene expression in SPS-induced rat PTSD model

Statistical analysis of HIF-1α, HIF-2α, HIF-3α, PACAP and PAI-1
gene expression in hippocampus and cortex of rats with SPS-induced
PTSD model. Includes normality testing, ANOVA/Kruskal-Wallis with
post-hoc comparisons, and publication-ready boxplot visualizations.

## Setup
```bash
pip install numpy scipy matplotlib openpyxl
```

## Data files

Place in `data/` folder:
- `PTSD_Hypoxia_Cortex.xlsx`
- `PTSD_Hypoxia_Hipo.xlsx`

Each file contains 5 sheets (HIF1, HIF2, HIF3, PACAP, PAI1).
Data starts from row 7. Groups are arranged in columns (6 columns per group):
- Control (offset 0)
- PTSD (offset 6)

## Usage
```bash
# With IQR outlier removal (default)
python qpcr_art1.py

# Without outlier removal
python qpcr_art1.py --no-outliers
```

## Statistical pipeline

1. **Outlier removal** — IQR 1.5× per group
2. **Normality** — Shapiro-Wilk test per group
3. **Overall test** — one-way ANOVA (all normal) or Kruskal-Wallis
4. **Post-hoc** — Student's t-test (after ANOVA) or Mann-Whitney U (after Kruskal)

## Output
```
results/YYYY-MM-DD_HH-MM-SS/
    Art1_qPCR_Cortex.png / .pdf
    Art1_qPCR_Hippocampus.png / .pdf
    stats_summary.csv
```

`stats_summary.csv` contains per-gene results:
normality p-values, overall test, post-hoc p-value, significance,
mean ± SD ± SEM, n, outliers removed.

## Citation

Porkhalo D. et al. (2025) — manuscript in preparation