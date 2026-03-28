# HIF and PACAP gene expression in SPS-induced rat PTSD model

Statistical analysis of HIF-1α, HIF-2α, HIF-3α, PACAP and PAI-1
gene expression in hippocampus and medial prefrontal cortex of rats with SPS-induced
PTSD model. Includes normality testing, ANOVA/Kruskal-Wallis with
post-hoc comparisons, and publication-ready boxplot visualizations.

## Setup
```bash
pip install numpy scipy matplotlib openpyxl
```

## Data files

Place in `data/` folder:

**Molecular data (qPCR):**
- `PCR_Results_Cortex.xlsx`
- `PCR_Results_Hippocampus.xlsx`

Each file contains 5 sheets (HIF-1, HIF-2, HIF-3, PACAP, PAI-1).
Data starts from row 10. Groups are arranged in columns (6 columns per group):
- Control (offset 0)
- PTSD (offset 6)

**Behavioural data:**
- `Open_Field_Test_Results.xlsx`
- `Elevated_Plus_Maze_Test_Results.xlsx`
- `Dark-Light_Box_Test_Results.xlsx`

Each file contains columns: `animal №`, `group`, `treatment`, followed by measured parameters.
Groups: `control` / `ptsd`. Treatment: `non`.

## Usage
```bash
# Molecular analysis
python molecular.py
python molecular.py --no-outliers

# Behavioural analysis (all tests)
python behavioural.py
python behavioural.py --test EPM
python behavioural.py --no-outliers
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
    Art1_OF.png / .pdf
    Art1_OF_stats.csv
    Art1_EPM.png / .pdf
    Art1_EPM_stats.csv
    Art1_DLB.png / .pdf
    Art1_DLB_stats.csv
```

CSV files contain per-parameter results:
normality p-values, overall test, post-hoc p-value, significance,
mean ± SD ± SEM, n, outliers removed.

## Citation

Porkhalo D. et al. (2025) — manuscript in preparation
