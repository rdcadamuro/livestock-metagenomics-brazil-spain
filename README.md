# Genome-resolved metagenomics of livestock waste in Brazil and Spain: resistome structure, mobilome, and phage-host associations


## Data availability

| Dataset | Repository | Accession |
|---|---|---|
| Brazilian metagenomes (MAGs + raw reads) | ENA | [PRJEB63971](https://www.ebi.ac.uk/ena/browser/view/PRJEB63971) |
| Spanish uViGs (generated in this study) | ENA | [PRJEB110911](https://www.ebi.ac.uk/ena/browser/view/PRJEB110911) |
| Spanish raw metagenomes (source) | NCBI SRA | [PRJNA628671](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA628671) |

---

## Repository structure

```
livestock-metagenomics-brazil-spain/
│
├── 01_uvigs_reconstruction/     # geNomad + vRhyme + PhaBox2 pipeline
├── 02_uvigs_taxonomy/           # vConTACT3 inputs, taxonomy barplots
├── 03_mag_figures/              # MAG ARG, taxonomy, and function figures
├── 04_mobilome/                 # tRNA, IS, integron detection (MAGs and uViGs)
├── 05_antiphage_defense/        # CRISPR-Cas, PADLOC, spacer-host links
├── 06_dramv_amg/                # DRAM-v AMG annotation and figures
├── 07_uvigs_quality/            # uViG confidence scoring
├── 08_ena_submission/           # ENA MAG submission workflow
└── 09_pipeline_bash/            # Main bash pipeline (binning, CheckM2, GTDB-Tk)
```

---

## Script descriptions

### 01 — uViGs reconstruction

| Script | Description |
|---|---|
| `run_refresh_genomad_vrhyme_brazil.py` | Runs geNomad viral identification + vRhyme binning for Brazilian samples |
| `run_refresh_genomad_vrhyme_spain.py` | Same pipeline for Spanish samples |
| `run_phabox_task.py` | Runs PhaBox2 (end-to-end: lifestyle prediction + host prediction via CHERRY) |
| `move_not_phages_to_sample_dirs.py` | Organizes non-phage contigs after geNomad filtering |
| `run_provirus_flanks.py` | Extracts flanking regions of integrated proviruses |

### 02 — uViGs taxonomy

| Script | Description |
|---|---|
| `build_vcontact3_inputs_from_vrhyme_bins.py` | Prepares protein FASTA and gene-to-genome TSV inputs for vConTACT3 |
| `run_vcontact3_wrapper.sh` | Wrapper to run vConTACT3 in a conda environment |
| `build_spain_final_taxonomy_outputs.py` | Consolidates vConTACT3 + PhaBox2 outputs for Spanish uViGs |
| `build_uvig_taxonomy_barplots_vcontact3.py` | Generates ICTV-level taxonomy barplots (class/family/genus) |
| `build_uvig_taxonomy_barplots_article2.py` | Extended taxonomy barplots for article figures |

### 03 — MAG figures

| Script | Description |
|---|---|
| `build_fig2_mag_taxon_function_brazil.py` | Builds Figure 4 panels: ARG, VF, plasmid, and taxonomy profiles for Brazilian MAGs |
| `build_fig4_arg_mag_vs_miuvig_brazil.py` | Builds Figure comparing ARG profiles between MAGs and uViGs |
| `taxonomy_4tools_consensus.py` | Consensus taxonomic classification across four tools (GTDB-Tk, GUNC, etc.) |

### 04 — Mobilome

| Script | Description |
|---|---|
| `Loop_Element_Finder_tRNA_Mags.py` | Runs tRNAscan-SE and MobileElementFinder across all MAG FASTAs |
| `Loop_Element_Finder_tRNA_MIUViGs.py` | Same pipeline for uViG FASTAs |
| `loop_mobile_element_IS_TrnaScan.py` | Integron and IS detection loop using IntegronFinder and ISsaga |

### 05 — Anti-phage defense

| Script | Description |
|---|---|
| `anti_phage_intrasample_and_context_stats.py` | Aggregates PADLOC + CRISPRCasFinder outputs per MAG bin and sample context |
| `integrate_gtdb_antiphage_by_genus.py` | Links anti-phage defense profiles to GTDB-Tk genus-level taxonomy |
| `paper_ready_antiphage_figure_suite.py` | Generates Figure 6 panels (defense heatmaps + spacer-host matrix) |
| `assemble_fig3_antiphage_spacer_panel.py` | Assembles spacer-to-protospacer match panels |

### 06 — DRAM-v / AMG

| Script | Description |
|---|---|
| `build_dramv_miuvig_final_tables.py` | Builds final AMG annotation tables from DRAM-v outputs |
| `map_dramv_scaffold_to_miuvig_bins.py` | Maps DRAM-v scaffold annotations back to vRhyme bin IDs |
| `plot_figure5_dramv_miuvig.py` | Generates Supplementary Figure S2 (AMG evidence tiers and functions) |
| `sync_dramv_brazil_to_results.py` | Syncs DRAM-v outputs from compute node to results directory |

### 07 — uViGs quality

| Script | Description |
|---|---|
| `vmag_confidence_score.py` | Computes a composite confidence score for each uViG bin (completeness, contamination, taxonomy) |

### 08 — ENA submission

| Script | Description |
|---|---|
| `build_final_submission_729.py` | Builds final ENA submission tables for 729 bacterial MAGs |
| `generate_submission_files.py` | Generates manifest TSV, sample sheet, and assembly files for ENA upload |
| `download_fastas_batch.sh` | Downloads and organizes FASTA files in batch from ENA |
| `download_and_verify_fastas.py` | Verifies MD5 checksums and completeness of downloaded FASTAs |
| `submit_all_assemblies.sh` | Submits all MAG assemblies to ENA via Webin-CLI |

### 09 — Bash pipeline

| Script | Description |
|---|---|
| `run_35_57_vmags_full.sh` | Full bash pipeline: geNomad → vRhyme → CheckV → PhaBox2 → BLAST → VIRIDIC → taxmyphage → Pharokka per sample |

---

## Dependencies

| Tool | Version | Used in |
|---|---|---|
| geNomad | 1.12.0 | 01, 09 |
| vRhyme | 1.1.0 | 01, 09 |
| CheckV | (end-to-end) | 01, 09 |
| PhaBox2 | 2.1.11 | 01, 09 |
| vConTACT3 | 3.2.0 | 02 |
| GTDB-Tk | 2.6.1 | 03, 05 |
| GUNC | 1.0.6 | 09 |
| CheckM2 | 1.1.0 | 09 |
| tRNAscan-SE | 2.0.12 | 04 |
| IntegronFinder | 2.0.6 | 04 |
| MobileElementFinder | 1.0.3 | 04 |
| MOB-suite | 3.1.9 | 04 |
| CRISPRCasFinder | — | 05 |
| PADLOC | — | 05 |
| DRAM-v | 1.5.0 | 06 |
| ABRicate | 1.2.0 | 03 |
| Python | ≥ 3.9 | all .py |
| pandas, matplotlib, seaborn, scipy | — | all .py |

---

## Configuration — setting your base paths

All scripts use path variables defined at the top of each file (clearly marked with `# USER: set your base path here`). Before running any script, update these variables to match your local directory structure. No paths are embedded in the logic of the scripts.

Example (Python):
```python
# USER: set your base path here
BASE = Path("/your/local/path/to/results")
```

Example (bash):
```bash
# USER: set your base path here
ASSEMBLY_BASE="/your/local/path/to/assembly"
```

---

## Citation

> Cadamuro RD, Soratto T, Cañete Reyes Á, Wagner G, Andersson B, Viancelli A, Michelon W, Rogovski P, Rodríguez-Lázaro D, Fongaro G. Genome-resolved metagenomics of livestock waste in Brazil and Spain: resistome structure, mobilome, and phage-host associations. *Environmental Microbiome* (submitted).

---

## Contact

Rafael Dorighello Cadamuro — cadamuro.rafael@gmail.com
