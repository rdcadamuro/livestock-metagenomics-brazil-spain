#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# USER: set your base path here

OUT_BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Article2_Brazil')
OUT_TSV = OUT_BASE / 'TSV'
OUT_FIG = OUT_BASE / 'FIG'

# USER: set your base path here

BASE_ANTI = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/Results/Spacers_AntiPhage_GTDB_Final_1035_1057')
# USER: set your base path here
BASE_SP = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/OneDrive_1_21-02-2026/Spacers_All/Results_Taxon_Spacers')

IMG_A = BASE_ANTI / 'IMG' / 'Fig_Context_Load_Broiler_vs_Biodigestor_Strict.png'
IMG_B = BASE_ANTI / 'IMG' / 'Fig_Biodigestor_InOut_Paired_Medians.png'
IMG_C = BASE_ANTI / 'IMG' / 'Fig_Genus_vs_AntiPhageSystem_Heatmap.png'
IMG_C_SVG = BASE_ANTI / 'IMG' / 'Fig_Genus_vs_AntiPhageSystem_Heatmap.svg'
IMG_D = BASE_SP / 'IMG_CLEAN' / 'Fig_Spacer_Genus_Consensus_CCTK.png'
HEATMAP_C_TSV = BASE_ANTI / 'TSV' / 'anti_phage_genus_system_top_heatmap.tsv'

OUT_PNG = OUT_FIG / 'Fig3_AntiPhage_MAGs_plus_SpacerTaxon_Brazil.png'
OUT_SVG = OUT_FIG / 'Fig3_AntiPhage_MAGs_plus_SpacerTaxon_Brazil.svg'
OUT_SOURCE = OUT_TSV / 'Fig3_source_data.tsv'


def load_img(path: Path):
    if path.exists() and path.stat().st_size > 0:
        return mpimg.imread(path)
    return None


def _is_formal_genus(name: str) -> bool:
    s = str(name).strip()
    return bool(re.match(r'^[A-Z][a-z]+(?:-[a-z]+)?$', s))


def _is_gtdb_placeholder(name: str) -> bool:
    s = str(name).strip()
    if s == '' or s.lower().startswith('unclassified'):
        return True
    if re.search(r'\d', s):
        return True
    if s.startswith('UBA') or s.startswith('JAA'):
        return True
    return not _is_formal_genus(s)


def rebuild_panel_c_heatmap() -> None:
    if not HEATMAP_C_TSV.exists() or HEATMAP_C_TSV.stat().st_size == 0:
        return
    df = pd.read_csv(HEATMAP_C_TSV, sep='\t', low_memory=False)
    if 'bact_genus' not in df.columns:
        return

    mat = df.set_index('bact_genus')
    for c in mat.columns:
        mat[c] = pd.to_numeric(mat[c], errors='coerce').fillna(0.0)
    if mat.empty:
        return

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'svg.fonttype': 'none',
    })

    fig_w = max(10.0, 0.45 * mat.shape[1] + 4.2)
    fig_h = max(8.0, 0.30 * mat.shape[0] + 3.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(mat.values, aspect='auto', cmap='viridis', vmin=0, vmax=100)
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=65, ha='right', fontsize=8)
    ax.set_yticks(np.arange(mat.shape[0]))

    ylabels = []
    placeholder_flags = []
    for g in mat.index.tolist():
        is_placeholder = _is_gtdb_placeholder(g)
        placeholder_flags.append(is_placeholder)
        ylabels.append(f"{g}*" if is_placeholder else g)
    ax.set_yticklabels(ylabels, fontsize=8)
    for tick, is_placeholder in zip(ax.get_yticklabels(), placeholder_flags):
        if not is_placeholder:
            tick.set_fontstyle('italic')

    ax.set_xlabel('Anti-phage system')
    ax.set_ylabel('Bacterial genus')
    ax.set_title('Top genera x anti-phage systems (% MAGs in genus with system)')
    cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.02)
    cbar.set_label('% of MAGs in genus')

    # Keep this very short and clear as requested.
    fig.text(0.5, 0.006, '* GTDB placeholder/provisional label (not formal genus).',
             ha='center', va='bottom', fontsize=8)

    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(IMG_C, dpi=450, bbox_inches='tight')
    fig.savefig(IMG_C_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    OUT_TSV.mkdir(parents=True, exist_ok=True)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    rebuild_panel_c_heatmap()

    # Build source table
    src_rows = []
    src_map = {
        'A_context_boxplot': BASE_ANTI / 'TSV' / 'anti_phage_context_bin_table.tsv',
        'A_context_stats': BASE_ANTI / 'TSV' / 'anti_phage_context_stats.tsv',
        'B_paired': BASE_ANTI / 'TSV' / 'anti_phage_biodigester_paired_medians.tsv',
        'C_heatmap': BASE_ANTI / 'TSV' / 'anti_phage_genus_system_top_heatmap.tsv',
        'D_spacer_taxon': BASE_SP / 'spacer_taxon_abundance_genus_pct.tsv',
        'D_spacer_unified': BASE_SP / 'spacer_taxonomy_unified_by_spacer.tsv',
    }

    for panel, fp in src_map.items():
        if not fp.exists() or fp.stat().st_size == 0:
            src_rows.append({'panel': panel, 'source_file': str(fp), 'n_rows': 0, 'status': 'missing'})
            continue
        try:
            n = max(sum(1 for _ in fp.open()) - 1, 0)
        except Exception:
            n = 0
        src_rows.append({'panel': panel, 'source_file': str(fp), 'n_rows': int(n), 'status': 'ok'})

    pd.DataFrame(src_rows).to_csv(OUT_SOURCE, sep='\t', index=False)

    # Compose panel image
    imgs = [load_img(IMG_A), load_img(IMG_B), load_img(IMG_C), load_img(IMG_D)]
    titles = [
        'A) Anti-phage load by context (strict)',
        'B) Biodigestor paired in/out',
        'C) Anti-phage systems x bacterial genus',
        'D) Spacer taxonomy (consensus genus)',
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    for i, ax in enumerate(axes.flatten()):
        if imgs[i] is None:
            ax.text(0.5, 0.5, f'Missing image:\n{[IMG_A, IMG_B, IMG_C, IMG_D][i]}', ha='center', va='center')
            ax.set_axis_off()
            continue
        ax.imshow(imgs[i])
        ax.set_title(titles[i], fontsize=11)
        ax.axis('off')

    fig.suptitle('Figure 3 - Anti-phage systems in MAGs + spacer taxonomy (Brazil-only)', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_PNG, dpi=450, bbox_inches='tight')
    fig.savefig(OUT_SVG, format='svg', bbox_inches='tight')
    plt.close(fig)

    print(f'[OK] figure: {OUT_PNG}')
    print(f'[OK] figure: {OUT_SVG}')
    print(f'[OK] source: {OUT_SOURCE}')


if __name__ == '__main__':
    main()
