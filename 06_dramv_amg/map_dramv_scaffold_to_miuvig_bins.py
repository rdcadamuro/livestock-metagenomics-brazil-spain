#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import re
from pathlib import Path
import pandas as pd

SAMPLES = [str(x) for x in range(1035, 1058)]
# USER: set your base path here
DRAM_BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/DRamV')
# USER: set your base path here
ANNO_BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/vRhyme_test_P19109_1035/Anotation')

SRC_FILES = {
    'supporting': 'AMG_supporting_evidence.tsv',
    'strict': 'AMG_strict_high_confidence.tsv',
    'potential': 'AMG_potential_transfer.tsv',
}

OUT_ALL = DRAM_BASE / 'DRAMv_MIUViG_scaffold_map_all_hits_1035_1057.tsv'
OUT_SUM = DRAM_BASE / 'DRAMv_MIUViG_scaffold_map_summary_1035_1057.tsv'


def normalize_scaffold(scaffold: object) -> str:
    s = str(scaffold).strip()
    if not s or s.lower() == 'nan':
        return ''
    return re.sub(r'__full-cat_\d+$', '', s)


def load_miuvig_index(sample: str) -> tuple[dict[str, list[tuple[str, str]]], str]:
    fasta = ANNO_BASE / sample / '_sample_level' / f'P19109_{sample}.vMAGs.clean.merged.fasta'
    if not fasta.exists() or fasta.stat().st_size == 0:
        return {}, f'missing_miuvig_fasta:{fasta}'

    idx: dict[str, list[tuple[str, str]]] = {}
    with fasta.open() as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            header = line[1:].strip().split()[0]
            m = re.match(r'^vRhyme_(\d+)__(.+)$', header)
            if not m:
                continue
            bin_id = f"vRhyme_bin_{m.group(1)}"
            contig_base = m.group(2)
            idx.setdefault(contig_base, []).append((header, bin_id))
    return idx, ''


def placeholder_row(sample: str, source_table: str, status: str, note: str) -> dict:
    return {
        'sample': sample,
        'source_table': source_table,
        'gene': 'NA',
        'scaffold_raw': 'NA',
        'scaffold_base': 'NA',
        'miuvig_contig_header': 'NA',
        'miuvig_bin': 'NA',
        'map_status': status,
        'map_note': note,
    }


def main() -> None:
    rows: list[dict] = []
    summary_rows: list[dict] = []

    for sample in SAMPLES:
        dram_dir = DRAM_BASE / f'P19109_{sample}'
        idx, idx_note = load_miuvig_index(sample)

        if not dram_dir.exists():
            for source in SRC_FILES:
                rows.append(placeholder_row(sample, source, 'missing_dramv_sample', 'missing DRAMv sample dir'))
                summary_rows.append(
                    {
                        'sample': sample,
                        'source_table': source,
                        'n_hits_total': 0,
                        'n_mapped': 0,
                        'n_unmapped': 0,
                        'n_ambiguous': 0,
                        'map_rate_pct': 0.0,
                        'status': 'missing_dramv_sample',
                        'note': 'missing DRAMv sample dir',
                    }
                )
            continue

        for source, fname in SRC_FILES.items():
            src = dram_dir / fname
            if not src.exists() or src.stat().st_size == 0:
                rows.append(placeholder_row(sample, source, 'missing_source_table', f'missing file:{src.name}'))
                summary_rows.append(
                    {
                        'sample': sample,
                        'source_table': source,
                        'n_hits_total': 0,
                        'n_mapped': 0,
                        'n_unmapped': 0,
                        'n_ambiguous': 0,
                        'map_rate_pct': 0.0,
                        'status': 'missing_source_table',
                        'note': f'missing file:{src.name}',
                    }
                )
                continue

            df = pd.read_csv(src, sep='\t', low_memory=False)
            if df.empty:
                rows.append(placeholder_row(sample, source, 'empty_source_table', f'empty file:{src.name}'))
                summary_rows.append(
                    {
                        'sample': sample,
                        'source_table': source,
                        'n_hits_total': 0,
                        'n_mapped': 0,
                        'n_unmapped': 0,
                        'n_ambiguous': 0,
                        'map_rate_pct': 0.0,
                        'status': 'empty_source_table',
                        'note': f'empty file:{src.name}',
                    }
                )
                continue

            n_total = len(df)
            n_mapped = 0
            n_unmapped = 0
            n_amb = 0

            for _, r in df.iterrows():
                scaffold_raw = r.get('scaffold', 'NA')
                scaffold_base = normalize_scaffold(scaffold_raw)

                out = {k: r.get(k, 'NA') for k in df.columns}
                out['sample'] = sample
                out['source_table'] = source
                out['scaffold_raw'] = scaffold_raw
                out['scaffold_base'] = scaffold_base

                if idx_note:
                    out['miuvig_contig_header'] = 'NA'
                    out['miuvig_bin'] = 'NA'
                    out['map_status'] = 'missing_miuvig_sample'
                    out['map_note'] = idx_note
                    n_unmapped += 1
                else:
                    matches = idx.get(scaffold_base, [])
                    if len(matches) == 1:
                        out['miuvig_contig_header'] = matches[0][0]
                        out['miuvig_bin'] = matches[0][1]
                        out['map_status'] = 'mapped'
                        out['map_note'] = ''
                        n_mapped += 1
                    elif len(matches) > 1:
                        out['miuvig_contig_header'] = ';'.join([x[0] for x in matches])
                        out['miuvig_bin'] = ';'.join(sorted(set([x[1] for x in matches])))
                        out['map_status'] = 'ambiguous_multi_match'
                        out['map_note'] = 'multiple merged headers matched same base contig'
                        n_amb += 1
                    else:
                        out['miuvig_contig_header'] = 'NA'
                        out['miuvig_bin'] = 'NA'
                        out['map_status'] = 'unmapped'
                        out['map_note'] = 'scaffold base not found in MIUViG merged fasta'
                        n_unmapped += 1

                rows.append(out)

            status = 'ok'
            note = ''
            if idx_note:
                status = 'missing_miuvig_sample'
                note = idx_note

            summary_rows.append(
                {
                    'sample': sample,
                    'source_table': source,
                    'n_hits_total': int(n_total),
                    'n_mapped': int(n_mapped),
                    'n_unmapped': int(n_unmapped),
                    'n_ambiguous': int(n_amb),
                    'map_rate_pct': round((100.0 * n_mapped / n_total) if n_total else 0.0, 4),
                    'status': status,
                    'note': note,
                }
            )

    all_df = pd.DataFrame(rows)
    sum_df = pd.DataFrame(summary_rows).sort_values(['sample', 'source_table'])

    key_cols = ['sample', 'source_table', 'gene', 'scaffold_raw']
    for c in key_cols:
        if c not in all_df.columns:
            all_df[c] = 'NA'
    all_df = all_df.drop_duplicates(subset=key_cols, keep='first')

    all_df.to_csv(OUT_ALL, sep='\t', index=False)
    sum_df.to_csv(OUT_SUM, sep='\t', index=False)

    print(f'[OK] all hits: {OUT_ALL}')
    print(f'[OK] summary: {OUT_SUM}')


if __name__ == '__main__':
    main()
