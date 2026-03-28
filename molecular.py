"""
qPCR Statistical Analysis (Control vs PTSD)
Boxplot with jitter, ggplot2 style.
Statistical pipeline:
    1. Outlier removal — IQR 1.5×
    2. Normality — Shapiro-Wilk per group
    3. Overall test — one-way ANOVA (all normal) or Kruskal-Wallis
    4. Post-hoc — t-test (after ANOVA) or Mann-Whitney U (after Kruskal)
Usage:
    python molecular.py
    python molecular.py --no-outliers
Data files in data/ folder:
    data/PCR_Results_Cortex.xlsx
    data/PCR_Results_Hippocampus.xlsx
"""
import os, csv, argparse
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import openpyxl
import warnings
warnings.filterwarnings('ignore')

# ── Config ──────────────────────────────────────────────────────────────────────
DATA_FILES = {
    'Cortex':      'data/PCR_Results_Cortex.xlsx',
    'Hippocampus': 'data/PCR_Results_Hippocampus.xlsx',
}

SHEETS      = ['HIF-1', 'HIF-2', 'HIF-3', 'PACAP', 'PAI-1']
GENE_LABELS = {
    'HIF-1': 'HIF-1α', 'HIF-2': 'HIF-2α', 'HIF-3': 'HIF-3α',
    'PACAP': 'PACAP',   'PAI-1': 'PAI-1',
}

CTRL_OFF = 0
PTSD_OFF = 6
DATA_ROW = 9

CTRL_COLOR = '#388E3C'
PTSD_COLOR = '#7B1FA2'
plt.rcParams.update({
    'font.family':  'DejaVu Sans',
    'pdf.fonttype': 42,
    'ps.fonttype':  42,
})

# ── Statistics ──────────────────────────────────────────────────────────────────
def remove_iqr(vals, factor=1.5):
    if len(vals) < 4: return vals
    q1, q3 = np.percentile(vals, 25), np.percentile(vals, 75)
    iqr = q3 - q1
    return [v for v in vals if q1 - factor*iqr <= v <= q3 + factor*iqr]

def shapiro_wilk(v):
    if len(v) < 3: return False, None
    _, p = stats.shapiro(v)
    return p >= 0.05, round(p, 4)

def run_stats(groups: dict) -> dict:
    group_names = list(groups.keys())
    values_list = list(groups.values())
    normality = {}
    for name, vals in groups.items():
        is_norm, p_sw = shapiro_wilk(vals)
        normality[name] = {'normal': is_norm, 'p': p_sw}
    all_normal = all(v['normal'] for v in normality.values())
    if all_normal:
        stat, p_overall = stats.f_oneway(*values_list)
        overall = {'test': 'One-way ANOVA', 'statistic': round(stat, 4),
                   'p': round(p_overall, 6)}
    else:
        stat, p_overall = stats.kruskal(*values_list)
        overall = {'test': 'Kruskal-Wallis', 'statistic': round(stat, 4),
                   'p': round(p_overall, 6)}
    posthoc = {}
    for i in range(len(group_names)):
        for j in range(i + 1, len(group_names)):
            g1, g2 = group_names[i], group_names[j]
            v1, v2 = groups[g1], groups[g2]
            if len(v1) < 3 or len(v2) < 3: continue
            if all_normal:
                _, p = stats.ttest_ind(v1, v2)
                test = 't-test'
            else:
                _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
                test = 'Mann-Whitney U'
            p = round(p, 6)
            sl = sig_label(p)
            posthoc[(g1, g2)] = {
                'test': test, 'p': p, 'sig': sl,
                'direction': (f'{g2} ↑' if np.mean(v2) > np.mean(v1)
                              else f'{g2} ↓') if sl not in ('ns', '') else '—'
            }
    return {'normality': normality, 'all_normal': all_normal,
            'overall': overall, 'posthoc': posthoc}

def sig_label(p):
    if p is None: return ''
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

def load_vals(filepath, sheet, offset):
    wb   = openpyxl.load_workbook(filepath, data_only=True)
    ws   = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    col  = offset + 4
    return [float(r[col]) for r in rows[DATA_ROW:]
            if col < len(r) and r[col] is not None
            and isinstance(r[col], (int, float))]

def print_stats(gene, sr, n_rm_ctrl=0, n_rm_ptsd=0):
    print(f"\n  {gene}:")
    if n_rm_ctrl or n_rm_ptsd:
        print(f"    ⚠ Outliers removed: Control={n_rm_ctrl}, PTSD={n_rm_ptsd}")
    print(f"    Normality (Shapiro-Wilk):")
    for grp, v in sr['normality'].items():
        status = '✅ normal' if v['normal'] else '❌ non-normal'
        p_str  = f"p={v['p']}" if v['p'] is not None else 'n/a'
        print(f"      {grp:<20} {p_str:<12} {status}")
    ov = sr['overall']
    print(f"    Overall: {ov['test']}  stat={ov['statistic']}  p={ov['p']}")
    print(f"    Post-hoc ({('t-test' if sr['all_normal'] else 'Mann-Whitney U')}):")
    for (g1, g2), v in sr['posthoc'].items():
        marker = ' ←' if v['sig'] not in ('ns', '') else ''
        print(f"      {g2} vs {g1:<18} p={v['p']:<10} "
              f"{v['sig']:<4} {v['direction']}{marker}")

# ── Plot ─────────────────────────────────────────────────────────────────────────
def plot_gene(ax, v1, v2, gene_label, posthoc_result, show_ylabel=False):
    bp = ax.boxplot(
        [v1, v2],
        positions=[0, 1],
        widths=0.45,
        patch_artist=True,
        notch=False,
        medianprops=dict(color='#212121', linewidth=2.0),
        whiskerprops=dict(color='#424242', linewidth=1.2),
        capprops=dict(color='#424242', linewidth=1.2),
        flierprops=dict(marker='', markersize=0),
        boxprops=dict(linewidth=1.2),
        zorder=3,
    )
    for patch, color in zip(bp['boxes'], [CTRL_COLOR, PTSD_COLOR]):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
        patch.set_edgecolor('#424242')
    np.random.seed(42)
    for xi, vals in [(0, v1), (1, v2)]:
        jitter = np.random.uniform(-0.13, 0.13, len(vals))
        ax.scatter(xi + jitter, vals,
                   color='#212121', s=28, alpha=0.75,
                   linewidths=0, zorder=5)
    ax.text(0, 0, f'n={len(v1)}', ha='center', va='top', fontsize=10,
            color='#757575', transform=ax.get_xaxis_transform())
    ax.text(1, 0, f'n={len(v2)}', ha='center', va='top', fontsize=10,
            color='#757575', transform=ax.get_xaxis_transform())
    data_max = max(v1 + v2)
    top = data_max * 1.35
    bottom = -1 if min(v1 + v2) >= 0 else min(v1 + v2) * 1.1
    ax.set_ylim(bottom, top)
    ax.set_xlim(-0.6, 1.6)
    ph = posthoc_result.get(('Control', 'PTSD'), {})
    sl = ph.get('sig', '')
    if sl not in ('ns', ''):
        y = top * 0.84
        h = top * 0.04
        ax.plot([0, 0, 1, 1], [y, y+h, y+h, y],
                lw=1.2, color='#212121', clip_on=False)
        ax.text(0.5, y + h*1.1, sl, ha='center', va='bottom',
                fontsize=14, color='#212121', fontweight='bold', clip_on=False)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Control', 'PTSD'],
                       fontsize=13, fontweight='bold', color='#212121')
    ax.set_title(gene_label, fontsize=13, fontweight='bold',
                 color='#212121', pad=7)
    if show_ylabel:
        ax.set_ylabel('[C]rel', fontsize=12, fontweight='bold',
                      color='#424242', labelpad=18)
    ax.set_facecolor('#EBEBEB')
    ax.yaxis.grid(True, color='white', linewidth=1.0, zorder=0)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color('#BDBDBD')
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(axis='x', length=0, pad=18)
    ax.tick_params(axis='y', length=0, labelsize=12)

def make_figure(tissue, data_dict, stats_dict, output_prefix):
    fig, axes = plt.subplots(1, 5, figsize=(20, 6))
    fig.patch.set_facecolor('white')
    for gi, sheet in enumerate(SHEETS):
        v1 = data_dict[sheet]['ctrl']
        v2 = data_dict[sheet]['ptsd']
        ph = stats_dict[sheet]['posthoc']
        plot_gene(axes[gi], v1, v2, GENE_LABELS[sheet], ph,
                  show_ylabel=(gi == 0))
    patches = [
        mpatches.Patch(facecolor=CTRL_COLOR, edgecolor='white', label='Control'),
        mpatches.Patch(facecolor=PTSD_COLOR, edgecolor='white', label='PTSD'),
    ]
    fig.legend(handles=patches, loc='lower center', ncol=2,
               fontsize=13, frameon=True, fancybox=False, edgecolor='#E0E0E0',
               bbox_to_anchor=(0.5, -0.04), handlelength=1.5, handleheight=1.0)
    fig.suptitle(f'Gene expression — {tissue} — Control vs PTSD',
                 fontsize=13, fontweight='bold', color='#212121', y=1.02)
    fig.text(0.5, -0.09,
             'Boxplot: median, IQR, whiskers = 1.5×IQR. '
             'Dots = individual animals (outliers removed, IQR 1.5×). '
             'Normality: Shapiro–Wilk. '
             'Overall: one-way ANOVA or Kruskal–Wallis. '
             'Post-hoc: t-test or Mann–Whitney U. '
             '* p<0.05, ** p<0.01, *** p<0.001.',
             ha='center', fontsize=8, color='#9E9E9E', style='italic')
    plt.tight_layout(rect=[0, 0.06, 1, 1])
    plt.subplots_adjust(wspace=0.40)
    fig.savefig(f'{output_prefix}.png', dpi=300,
                bbox_inches='tight', facecolor='white')
    fig.savefig(f'{output_prefix}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"\n  Figure → {output_prefix}.png / .pdf")

def save_csv(rows, filepath):
    if not rows: return
    fields = list(rows[0].keys())
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  CSV    → {filepath}")

# ── Main ─────────────────────────────────────────────────────────────────────────
def run(remove_outliers=True):
    ts      = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    out_dir = os.path.join('results', ts)
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"qPCR Analysis — Article 1 (Control vs PTSD)")
    print(f"Outlier removal: {'ON (IQR 1.5×)' if remove_outliers else 'OFF'}")
    print(f"Output → {out_dir}/")
    print('='*60)
    all_csv = []
    for tissue, filepath in DATA_FILES.items():
        print(f"\n{'─'*60}")
        print(f"  {tissue}")
        print(f"{'─'*60}")
        data_dict  = {}
        stats_dict = {}
        for sheet in SHEETS:
            gene  = GENE_LABELS[sheet]
            raw_c = load_vals(filepath, sheet, CTRL_OFF)
            raw_p = load_vals(filepath, sheet, PTSD_OFF)
            ctrl = remove_iqr(raw_c) if remove_outliers else raw_c
            ptsd = remove_iqr(raw_p) if remove_outliers else raw_p
            n_rm_ctrl = len(raw_c) - len(ctrl)
            n_rm_ptsd = len(raw_p) - len(ptsd)
            data_dict[sheet] = {'ctrl': ctrl, 'ptsd': ptsd}
            if not ctrl or not ptsd: continue
            sr = run_stats({'Control': ctrl, 'PTSD': ptsd})
            stats_dict[sheet] = sr
            print_stats(gene, sr, n_rm_ctrl, n_rm_ptsd)
            ph = sr['posthoc'].get(('Control', 'PTSD'), {})
            m1, m2 = np.mean(ctrl), np.mean(ptsd)
            s1, s2 = np.std(ctrl, ddof=1), np.std(ptsd, ddof=1)
            all_csv.append({
                'tissue':                tissue,
                'gene':                  gene,
                'normality_control_p':   sr['normality']['Control']['p'],
                'normality_control':     'normal' if sr['normality']['Control']['normal'] else 'non-normal',
                'normality_ptsd_p':      sr['normality']['PTSD']['p'],
                'normality_ptsd':        'normal' if sr['normality']['PTSD']['normal'] else 'non-normal',
                'overall_test':          sr['overall']['test'],
                'overall_p':             sr['overall']['p'],
                'posthoc_test':          ph.get('test', ''),
                'posthoc_p':             ph.get('p', ''),
                'significance':          ph.get('sig', ''),
                'direction':             ph.get('direction', ''),
                'mean_ctrl':             round(m1, 4),
                'sd_ctrl':               round(s1, 4),
                'sem_ctrl':              round(s1/np.sqrt(len(ctrl)), 4),
                'n_ctrl':                len(ctrl),
                'mean_ptsd':             round(m2, 4),
                'sd_ptsd':               round(s2, 4),
                'sem_ptsd':              round(s2/np.sqrt(len(ptsd)), 4),
                'n_ptsd':                len(ptsd),
                'outliers_removed_ctrl': n_rm_ctrl,
                'outliers_removed_ptsd': n_rm_ptsd,
            })
        make_figure(tissue, data_dict, stats_dict,
                    os.path.join(out_dir, f'Art1_qPCR_{tissue}'))
    save_csv(all_csv, os.path.join(out_dir, 'stats_summary.csv'))
    print(f"\nDone! All results → {out_dir}/")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-outliers', action='store_true',
                        help='Skip IQR outlier removal')
    args = parser.parse_args()
    run(remove_outliers=not args.no_outliers)