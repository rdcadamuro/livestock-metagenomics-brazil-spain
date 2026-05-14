#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import pandas as pd

SAMPLES = [str(x) for x in range(1035, 1058)]

# USER: set your base path here

MAP_ALL = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/DRamV/DRAMv_MIUViG_scaffold_map_all_hits_1035_1057.tsv')
# USER: set your base path here
DRAM_BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/DRamV')
SYNC_MANIFEST = DRAM_BASE / 'DRAMv_sync_manifest_1035_1057.tsv'
# USER: set your base path here
TAXON_TOTAL = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/vRhyme_test_P19109_1035/Anotation/Taxon_Total_Phages_vmags.tsv')
# USER: set your base path here
OUT_BASE = Path('/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Article2_Brazil')
OUT_TSV = OUT_BASE / 'TSV'


def sample_to_context(sample: str) -> dict:
    s = int(sample) % 100
    if 35 <= s <= 41:
        return {
            'PlotGroup': 'Broiler litter W/A',
            'Context': 'Broiler litter',
            'Host': 'CHICKEN',
            'Stage': 'NA',
            'Digester': 'NA',
        }
    if s == 42:
        return {
            'PlotGroup': 'Broiler litter No/A',
            'Context': 'Broiler litter',
            'Host': 'CHICKEN',
            'Stage': 'NA',
            'Digester': 'NA',
        }
    if s == 43:
        return {
            'PlotGroup': 'Biodigestor Swine Entrance No/A',
            'Context': 'Biodigestor',
            'Host': 'SWINE',
            'Stage': 'IN',
            'Digester': 'BD_NOAB',
        }
    if s >= 44 and s % 2 == 0:
        pair = (s - 44) // 2 + 1
        return {
            'PlotGroup': 'Biodigestor Swine In',
            'Context': 'Biodigestor',
            'Host': 'SWINE',
            'Stage': 'IN',
            'Digester': f'BD{pair:02d}',
        }
    if s >= 45 and s % 2 == 1:
        pair = (s - 45) // 2 + 1
        return {
            'PlotGroup': 'Biodigestor Swine Out',
            'Context': 'Biodigestor',
            'Host': 'SWINE',
            'Stage': 'OUT',
            'Digester': f'BD{pair:02d}',
        }
    return {
        'PlotGroup': 'NA',
        'Context': 'Other',
        'Host': 'NA',
        'Stage': 'NA',
        'Digester': 'NA',
    }


def tsv_count(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        return max(sum(1 for _ in path.open()) - 1, 0)
    except Exception:
        return 0


def main() -> None:
    OUT_TSV.mkdir(parents=True, exist_ok=True)
    sync = pd.read_csv(SYNC_MANIFEST, sep='\t')
    sync['sample'] = sync['sample'].astype(str)
    sync_status = {r['sample']: r['status'] for _, r in sync.iterrows()}

    m = pd.read_csv(MAP_ALL, sep='\t', low_memory=False)
    m['sample'] = m['sample'].astype(str)

    tax = pd.read_csv(TAXON_TOTAL, sep='\t', low_memory=False)
    tax['sample'] = tax['SampleName'].astype(str)
    tax['miuvig_bin'] = tax['Bin'].astype(str)

    keep_tax_cols = [
        'sample', 'miuvig_bin', 'N_contigs',
        'Taxon_level', 'Taxo_name', 'Taxon_source', 'Consensus_4tools',
        'Tool1_Pharokka', 'Tool2_PhaGenus', 'Tool3_PhaBOX', 'Tool4_BLASTn', 'Tool5_geNomad',
        'Tool6_TaxMyPhage', 'Tool7_VIRIDIC_RefAware',
    ]
    keep_tax_cols = [c for c in keep_tax_cols if c in tax.columns]
    tax2 = tax[keep_tax_cols].drop_duplicates(subset=['sample', 'miuvig_bin'], keep='first')

    mapped = m[m['map_status'] == 'mapped'].copy()
    mapped = mapped.merge(tax2, on=['sample', 'miuvig_bin'], how='left')

    ctx_df = pd.DataFrame([{'sample': s, **sample_to_context(s)} for s in SAMPLES])
    mapped = mapped.merge(ctx_df, on='sample', how='left')

    potential = mapped[mapped['source_table'] == 'potential'].copy()
    supporting = mapped[mapped['source_table'] == 'supporting'].copy()
    strict = mapped[mapped['source_table'] == 'strict'].copy()

    potential.to_csv(OUT_TSV / 'AMG_by_MIUViG_main_potential_transfer.tsv', sep='\t', index=False)
    supporting.to_csv(OUT_TSV / 'AMG_by_MIUViG_supp_supporting_evidence.tsv', sep='\t', index=False)
    strict.to_csv(OUT_TSV / 'AMG_by_MIUViG_strict_high_confidence.tsv', sep='\t', index=False)

    # Funnel counts by sample
    by_sample = []
    for s in SAMPLES:
        d = DRAM_BASE / f'P19109_{s}'
        st = sync_status.get(s, 'missing_remote')
        if st not in {'found_synced', 'already_local'}:
            by_sample.append(
                {
                    'sample': s,
                    'dramv_status': 'missing_pending',
                    'amg_summary_n': 0,
                    'supporting_n': 0,
                    'strict_n': 0,
                    'potential_n': 0,
                }
            )
            continue

        by_sample.append(
            {
                'sample': s,
                'dramv_status': 'available',
                'amg_summary_n': tsv_count(d / 'amg_summary.tsv'),
                'supporting_n': tsv_count(d / 'AMG_supporting_evidence.tsv'),
                'strict_n': tsv_count(d / 'AMG_strict_high_confidence.tsv'),
                'potential_n': tsv_count(d / 'AMG_potential_transfer.tsv'),
            }
        )

    funil_sample = pd.DataFrame(by_sample).merge(ctx_df, on='sample', how='left')
    funil_sample.to_csv(OUT_TSV / 'AMG_funil_counts_by_sample.tsv', sep='\t', index=False)

    overall = pd.DataFrame(
        [
            {
                'scope': 'Brazil_current',
                'n_samples_total': len(SAMPLES),
                'n_samples_with_dramv': int((funil_sample['dramv_status'] == 'available').sum()),
                'n_samples_missing_pending': int((funil_sample['dramv_status'] == 'missing_pending').sum()),
                'amg_summary_n_total': int(funil_sample['amg_summary_n'].sum()),
                'supporting_n_total': int(funil_sample['supporting_n'].sum()),
                'strict_n_total': int(funil_sample['strict_n'].sum()),
                'potential_n_total': int(funil_sample['potential_n'].sum()),
            }
        ]
    )
    overall.to_csv(OUT_TSV / 'AMG_funil_counts_overall.tsv', sep='\t', index=False)

    unmapped = m[
        (m['source_table'].isin(['supporting', 'strict', 'potential']))
        & (m['map_status'] != 'mapped')
    ].copy()
    cols = [c for c in ['sample', 'source_table', 'gene', 'scaffold_raw', 'scaffold_base', 'map_status', 'map_note'] if c in unmapped.columns]
    unmapped = unmapped[cols]
    unmapped.to_csv(OUT_TSV / 'AMG_unmapped_scaffolds_audit.tsv', sep='\t', index=False)

    print('[OK] wrote: AMG_by_MIUViG_main_potential_transfer.tsv')
    print('[OK] wrote: AMG_by_MIUViG_supp_supporting_evidence.tsv')
    print('[OK] wrote: AMG_by_MIUViG_strict_high_confidence.tsv')
    print('[OK] wrote: AMG_funil_counts_by_sample.tsv')
    print('[OK] wrote: AMG_funil_counts_overall.tsv')
    print('[OK] wrote: AMG_unmapped_scaffolds_audit.tsv')


if __name__ == '__main__':
    main()
