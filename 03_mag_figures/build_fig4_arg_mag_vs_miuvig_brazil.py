#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# USER: set your base path here

BASE_OUT = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Article2_Brazil')
TSV_OUT = BASE_OUT / 'TSV'
FIG_OUT = BASE_OUT / 'FIG'

# USER: set your base path here

MIUVIG_ARG = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/ARG_PHAGES_MIUVIG_80_80_DB5/ALL_SAMPLES_UNIFIED_CLASSIFIED.tsv')
# USER: set your base path here
MAG_ARG = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/ARGs/ALL_SAMPLES_UNIFIED_with_categoria_final.tsv')
# USER: set your base path here
GTDB_BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/GTDBK/Downloaded_ubu_GTDBK')

OUT_PNG = FIG_OUT / 'Fig4_ARG_MAG_vs_MIUViG_Brazil.png'
OUT_SVG = FIG_OUT / 'Fig4_ARG_MAG_vs_MIUViG_Brazil.svg'
OUT_SOURCE = TSV_OUT / 'Fig4_source_data.tsv'


def parse_genus(classification: str) -> str:
    if not isinstance(classification, str):
        return 'Unclassified_bacteria'
    genus = 'Unclassified_bacteria'
    for p in classification.split(';'):
        p = p.strip()
        if p.startswith('g__'):
            genus = p[3:] or 'Unclassified_bacteria'
    return genus


def load_mag_bin_genus(samples: list[str]) -> pd.DataFrame:
    rows = []
    for s in samples:
        fp = GTDB_BASE / f'P19109_{s}' / 'classify' / 'gtdbtk.bac120.summary.tsv'
        if not fp.exists() or fp.stat().st_size == 0:
            continue
        df = pd.read_csv(fp, sep='\t', low_memory=False)
        if 'user_genome' not in df.columns or 'classification' not in df.columns:
            continue
        for _, r in df.iterrows():
            rows.append({'sample': s, 'bin_norm': str(r['user_genome']), 'bact_genus': parse_genus(r['classification'])})
    return pd.DataFrame(rows)


def main() -> None:
    TSV_OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    mag = pd.read_csv(MAG_ARG, sep='\t', low_memory=False)
    miu = pd.read_csv(MIUVIG_ARG, sep='\t', low_memory=False)

    mag['sample'] = mag['Sample'].astype(str).str.replace('P19109_', '', regex=False)
    mag['bin_norm'] = mag['MAG_norm'].astype(str)
    mag['category'] = mag['CATEGORIA_FINAL'].fillna('Unclassified')

    samples = sorted(mag['sample'].unique(), key=lambda x: int(x))
    bg = load_mag_bin_genus(samples)
    if not bg.empty:
        mag = mag.merge(bg, on=['sample', 'bin_norm'], how='left')
    mag['bact_genus'] = mag['bact_genus'].fillna('Unclassified_bacteria')

    miu['sample'] = miu['Sample'].astype(str)
    miu['bin'] = miu['Bin'].astype(str)
    miu['category'] = miu['Category_final'].fillna('Unclassified')

    # Panel A metrics
    metrics = pd.DataFrame([
        {'dataset': 'MAG', 'metric': 'ARG_hits', 'value': int(len(mag))},
        {'dataset': 'MAG', 'metric': 'ARG_positive_bins', 'value': int(mag[['sample', 'bin_norm']].drop_duplicates().shape[0])},
        {'dataset': 'MIUViG', 'metric': 'ARG_hits', 'value': int(len(miu))},
        {'dataset': 'MIUViG', 'metric': 'ARG_positive_bins', 'value': int(miu[['sample', 'bin']].drop_duplicates().shape[0])},
    ])

    # Panel B heatmap matrix (genus x category)
    hm = (
        mag.groupby(['bact_genus', 'category'], dropna=False)
        .size()
        .reset_index(name='n_hits')
    )
    top_g = hm.groupby('bact_genus')['n_hits'].sum().sort_values(ascending=False).head(20).index.tolist()
    top_c = hm.groupby('category')['n_hits'].sum().sort_values(ascending=False).head(10).index.tolist()
    hm2 = hm[hm['bact_genus'].isin(top_g) & hm['category'].isin(top_c)].copy()
    mat = hm2.pivot(index='bact_genus', columns='category', values='n_hits').fillna(0)

    # Panel C MIUViG presence
    miu_presence = miu.groupby(['sample', 'bin', 'category'], dropna=False).size().reset_index(name='n_hits')
    miu_presence['sample_bin'] = miu_presence['sample'].astype(str) + '|' + miu_presence['bin'].astype(str)

    # source data
    src_rows = []
    for _, r in metrics.iterrows():
        src_rows.append({'panel': 'A', 'label': f"{r['dataset']}|{r['metric']}", 'value': int(r['value'])})
    for idx in mat.index:
        for col in mat.columns:
            src_rows.append({'panel': 'B', 'label': f'{idx}|{col}', 'value': float(mat.loc[idx, col])})
    for _, r in miu_presence.iterrows():
        src_rows.append({'panel': 'C', 'label': f"{r['sample_bin']}|{r['category']}", 'value': int(r['n_hits'])})
    pd.DataFrame(src_rows).to_csv(OUT_SOURCE, sep='\t', index=False)

    # plot
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'svg.fonttype': 'none',
    })

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    axA, axB, axC = axes

    # A
    aplot = metrics.pivot(index='metric', columns='dataset', values='value').fillna(0)
    x = range(len(aplot.index))
    w = 0.35
    axA.bar([i - w / 2 for i in x], aplot['MAG'].values, width=w, label='MAG', color='#4C78A8')
    axA.bar([i + w / 2 for i in x], aplot['MIUViG'].values, width=w, label='MIUViG', color='#F58518')
    axA.set_xticks(list(x))
    axA.set_xticklabels(aplot.index)
    axA.set_ylabel('count')
    axA.set_title('A) Overall ARG burden: MAG vs MIUViG')
    axA.legend(frameon=False)

    # B
    im = axB.imshow(mat.values, aspect='auto')
    axB.set_yticks(range(len(mat.index)))
    axB.set_yticklabels(mat.index)
    axB.set_xticks(range(len(mat.columns)))
    axB.set_xticklabels(mat.columns, rotation=45, ha='right')
    axB.set_title('B) MAG association: bacterial genus x ARG category')
    fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04)

    # C
    if miu_presence.empty:
        axC.text(0.5, 0.5, 'No ARG in MIUViG', ha='center', va='center')
        axC.set_axis_off()
    else:
        p = miu_presence.groupby('sample_bin')['n_hits'].sum().sort_values(ascending=False)
        axC.barh(p.index[::-1], p.values[::-1], color='#F58518')
        axC.set_xlabel('n ARG hits')
        axC.set_title('C) ARG presence in MIUViGs (sample|bin)')

    fig.suptitle('Figure 4 - ARG in MAGs vs MIUViGs (Brazil-only)', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_PNG, dpi=450, bbox_inches='tight')
    fig.savefig(OUT_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    print(f'[OK] figure: {OUT_PNG}')
    print(f'[OK] figure: {OUT_SVG}')
    print(f'[OK] source: {OUT_SOURCE}')


if __name__ == '__main__':
    main()
