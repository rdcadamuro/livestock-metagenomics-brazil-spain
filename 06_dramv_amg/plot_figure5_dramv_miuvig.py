#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# USER: set your base path here

BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Article2_Brazil')
TSV = BASE / 'TSV'
FIG = BASE / 'FIG'

IN_POTENTIAL = TSV / 'AMG_by_MIUViG_main_potential_transfer.tsv'

OUT_PNG = FIG / 'Fig5_DRAMv_AMG_Funnel_MIUViG_Brazil.png'
OUT_SVG = FIG / 'Fig5_DRAMv_AMG_Funnel_MIUViG_Brazil.svg'
OUT_SOURCE = TSV / 'Fig5_source_data.tsv'

GROUP_ORDER = [
    'Broiler litter W/A',
    'Broiler litter No/A',
    'Biodigestor Swine Entrance No/A',
    'Biodigestor Swine In',
    'Biodigestor Swine Out',
]
GROUP_COL = {
    'Broiler litter W/A': '#E69F00',
    'Broiler litter No/A': '#C44E00',
    'Biodigestor Swine Entrance No/A': '#009E73',
    'Biodigestor Swine In': '#56B4E9',
    'Biodigestor Swine Out': '#CC79A7',
}


def _is_na(x: object) -> bool:
    s = str(x).strip()
    return s == '' or s.lower() in {'na', 'nan', 'none', 'null'}


def _pick_function_label(row: pd.Series) -> str:
    # Priority: curated pathway/module annotation first, then KO/gene fallback.
    for col in ['module', 'kegg_hit', 'gene_description', 'ko_id', 'gene']:
        if col in row and not _is_na(row[col]):
            val = str(row[col]).strip()
            val = val.replace(';', ' | ')
            if len(val) > 90:
                return val[:87] + '...'
            return val
    return 'Unknown_function'


def _simplify_function_label(label: str) -> str:
    s = str(label).strip().lower()
    if ('pyrimidine' in s) or ('m00053' in s) or ('k01520' in s):
        return 'DNA nucleotide synthesis (pyrimidines)'
    if ('methionine' in s and 'cysteine' in s) or ('k00789' in s):
        return 'Sulfur amino-acid metabolism (Met/Cys)'
    if ('methionine degradation' in s) or ('k00525' in s) or ('k00558' in s):
        return 'Methionine metabolism'
    if s == '' or s in {'na', 'nan', 'none', 'null'}:
        return 'Unknown function'
    return label


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'svg.fonttype': 'none',
    })

    pot = pd.read_csv(IN_POTENTIAL, sep='\t', low_memory=False)

    # source data export (long format)
    long_rows = []
    if not pot.empty:
        # New Panel A source: potential transfer counts by context.
        pctx = pot.groupby('PlotGroup', dropna=False).size().reset_index(name='count')
        for _, r in pctx.iterrows():
            long_rows.append({
                'panel': 'A_context_potential',
                'sample': 'NA',
                'PlotGroup': str(r['PlotGroup']),
                'stage': 'potential',
                'count': int(r['count']),
            })

        # New Panel B source: functional action heatmap (function x context).
        p2 = pot.copy()
        p2['function_action'] = p2.apply(_pick_function_label, axis=1)
        p2['function_action_simple'] = p2['function_action'].map(_simplify_function_label)
        pfun = p2.groupby(['function_action', 'PlotGroup'], dropna=False).size().reset_index(name='count')
        for _, r in pfun.iterrows():
            long_rows.append({
                'panel': 'B_function_heatmap',
                'sample': 'NA',
                'PlotGroup': str(r['PlotGroup']),
                'stage': str(_simplify_function_label(r['function_action'])),
                'count': int(r['count']),
            })

    pd.DataFrame(long_rows).to_csv(OUT_SOURCE, sep='\t', index=False)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    axA, axB = axes

    # New A: potential by context (old C)
    if pot.empty:
        axA.text(0.5, 0.5, 'No mapped potential hits', ha='center', va='center')
        axA.set_axis_off()
    else:
        pctx = pot.groupby('PlotGroup', dropna=False).size().reindex(GROUP_ORDER).fillna(0).astype(int)
        colors = [GROUP_COL[g] for g in GROUP_ORDER]
        axA.bar(range(len(GROUP_ORDER)), pctx.values, color=colors)
        axA.set_xticks(range(len(GROUP_ORDER)))
        axA.set_xticklabels(GROUP_ORDER, rotation=20, ha='right')
        axA.set_ylabel('n potential AMG genes')
        axA.set_title('A) Potential transfer by context')
        for i, v in enumerate(pctx.values):
            axA.text(i, v + 0.05, str(int(v)), ha='center', va='bottom', fontsize=8)

    # New B: functional action promoted by potential AMG (heatmap)
    if pot.empty:
        axB.text(0.5, 0.5, 'No mapped potential hits', ha='center', va='center')
        axB.set_axis_off()
    else:
        p2 = pot.copy()
        p2['function_action'] = p2.apply(_pick_function_label, axis=1)
        p2['function_action_simple'] = p2['function_action'].map(_simplify_function_label)

        top_funcs = (
            p2['function_action_simple']
            .value_counts()
            .head(12)
            .index
            .tolist()
        )
        mat = (
            p2[p2['function_action_simple'].isin(top_funcs)]
            .groupby(['function_action_simple', 'PlotGroup'], dropna=False)
            .size()
            .reset_index(name='count')
            .pivot(index='function_action_simple', columns='PlotGroup', values='count')
            .reindex(index=top_funcs, columns=GROUP_ORDER)
            .fillna(0)
        )

        im = axB.imshow(mat.values, aspect='auto', interpolation='nearest')
        axB.set_yticks(np.arange(len(mat.index)))
        axB.set_yticklabels(mat.index)
        axB.set_xticks(np.arange(len(mat.columns)))
        axB.set_xticklabels(mat.columns, rotation=20, ha='right')
        axB.set_title('B) AMG functional action in MIUViGs (potential transfer)')
        axB.set_xlabel('Context')
        axB.set_ylabel('AMG function/action')

        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = int(mat.iloc[i, j])
                if v > 0:
                    axB.text(j, i, str(v), ha='center', va='center', fontsize=8, color='white')
        cbar = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04)
        cbar.set_label('n potential AMG genes')

    fig.suptitle('Figure 5 - DRAM-v potential AMG: context distribution and functional action (Brazil-only)', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_PNG, dpi=450, bbox_inches='tight')
    fig.savefig(OUT_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    print(f'[OK] figure: {OUT_PNG}')
    print(f'[OK] figure: {OUT_SVG}')
    print(f'[OK] source: {OUT_SOURCE}')


if __name__ == '__main__':
    main()
