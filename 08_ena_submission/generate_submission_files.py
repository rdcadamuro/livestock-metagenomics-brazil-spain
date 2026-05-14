#!/usr/bin/env python3
"""
Gera 3 ficheiros de submissão ENA para os 729 MAGs:

  1. ena_taxonomy_registration_request.tsv  — 4 nomes para registar no ENA
  2. ena_biosample_sheet_729.tsv            — sample sheet GSC MIMAGs (Webin Portal)
  3. ena_assembly_manifests/               — 1 manifest por MAG para webin-cli

Uso: python3 generate_submission_files.py
"""

import csv, os, re
from pathlib import Path
from collections import defaultdict

# USER: set your base path here

BASE     = Path("/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Mags_Vmags")
SUBDIR   = BASE / "ENA_final_submission_729"
TABLE    = SUBDIR / "ena_mag_candidates_final_729_taxfixed.tsv"
FASTA_BR = SUBDIR / "fastas/Brazil"
FASTA_SP = SUBDIR / "fastas/Spain"
OUT      = SUBDIR / "submission_files"
MANIFEST_DIR = OUT / "manifests"

OUT.mkdir(exist_ok=True)
MANIFEST_DIR.mkdir(exist_ok=True)

# ── Load table ──────────────────────────────────────────────────────────────
rows = []
with open(TABLE) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)
print(f"MAGs carregados: {len(rows)}")

# ── Helper: resolve local FASTA path ───────────────────────────────────────
def local_fasta(mag_id, dataset):
    if dataset == "Brazil":
        p = FASTA_BR / f"{mag_id}.fa"
    else:
        p = FASTA_SP / f"{mag_id}.fa"
    return str(p) if p.exists() else ""

# ── FILE 1: ENA Taxonomy Registration Request (4 names) ────────────────────
print("\n[1/3] Gerando taxonomy registration request...")

needs_reg = [r for r in rows if r["tax_id_status"] == "needs_ena_registration"]

tax_reg_file = OUT / "ena_taxonomy_registration_request.tsv"
tax_fields = [
    "proposed_name", "gtdb_taxonomy", "ncbi_phylum", "ncbi_phylum_taxid",
    "description", "affected_mag_ids", "n_mags"
]

# Group by proposed name
by_name = defaultdict(list)
for r in needs_reg:
    by_name[r["scientific_name_for_ena"]].append(r)

with open(tax_reg_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=tax_fields, delimiter="\t")
    writer.writeheader()
    for name, mags in sorted(by_name.items()):
        r0 = mags[0]
        # Parse phylum from GTDB
        gtdb = r0.get("gtdb_taxonomy", "")
        phylum_gtdb = ""
        for part in gtdb.split(";"):
            if part.strip().startswith("p__"):
                phylum_gtdb = part.strip()[3:]
        # GTDB→NCBI phylum mapping
        gtdb_to_ncbi_phylum = {
            "Bacteroidota": ("Bacteroidetes", "976"),
            "Bacillota": ("Firmicutes", "1239"),
            "Pseudomonadota": ("Proteobacteria", "1224"),
            "Actinomycetota": ("Actinobacteria", "201174"),
        }
        ncbi_phylum, ncbi_phylum_taxid = gtdb_to_ncbi_phylum.get(phylum_gtdb, ("", ""))
        writer.writerow({
            "proposed_name":      name,
            "gtdb_taxonomy":      gtdb,
            "ncbi_phylum":        ncbi_phylum,
            "ncbi_phylum_taxid":  ncbi_phylum_taxid,
            "description":        (
                f"Uncultured representative of GTDB genus {r0.get('gtdb_genus','?')}. "
                f"Recovered from metagenome-assembled genome (MAG) from livestock wastewater. "
                f"GTDB classification: {gtdb}."
            ),
            "affected_mag_ids":   "; ".join(r["mag_id"] for r in mags),
            "n_mags":             len(mags),
        })

print(f"  -> {tax_reg_file.name} ({len(by_name)} nomes)")
for name, mags in sorted(by_name.items()):
    print(f"     {name}: {len(mags)} MAG(s)")

# ── FILE 2: BioSample sheet (GSC MIMAGs checklist ERC000047) ───────────────
print("\n[2/3] Gerando BioSample sheet...")

# ENA GSC MIMAGs checklist fields (ERC000047)
bs_fields = [
    "sample_alias",
    "tax_id",
    "scientific_name",
    "sample_title",
    "sample_description",
    "collection date",
    "geographic location (country and/or sea)",
    "geographic location (latitude)",
    "geographic location (longitude)",
    "isolation_source",
    "broad-scale environmental context",
    "local environmental context",
    "environmental medium",
    "completeness score",
    "contamination score",
    "assembly quality",
    "assembly software",
    "sequencing method",
    "investigation type",
    "project name",
    "derived from",
    "binning software",
    "taxonomic identity marker",
    "taxonomy reference database name",
    "taxonomy reference database version",
]

# Context by dataset
CONTEXT = {
    "Brazil": {
        "collection_date": "2018/2019",
        "country": "Brazil",
        "lat": "-15.7801",
        "lon": "-47.9292",
        "isolation_source": "broiler litter (poultry production facility)",
        "biome": "anthropogenic terrestrial biome [ENVO:01000219]",
        "feature": "poultry house [ENVO:01001874]",
        "material": "poultry litter [ENVO:00002191]",
        "study": "PRJEB110925",
    },
    "Spain": {
        "collection_date": "2020",
        "country": "Spain",
        "lat": "40.4637",
        "lon": "-3.7492",
        "isolation_source": "pig farm slurry (swine production facility)",
        "biome": "anthropogenic terrestrial biome [ENVO:01000219]",
        "feature": "animal farm [ENVO:00003040]",
        "material": "slurry [ENVO:00002090]",
        "study": "PRJEB110911",
    },
}

bs_file = OUT / "ena_biosample_sheet_729.tsv"
skipped = []

with open(bs_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=bs_fields, delimiter="\t")
    writer.writeheader()
    for r in rows:
        mag_id  = r["mag_id"]
        dataset = r["dataset"]
        ctx     = CONTEXT[dataset]

        tax_id  = r.get("tax_id", "").strip()
        sci_name = r.get("scientific_name_for_ena", "").strip() or r.get("taxon_name_for_submission","").strip()

        if r["tax_id_status"] == "needs_ena_registration":
            # Use placeholder taxid — must be updated after ENA registration
            tax_id   = "PENDING_ENA_REGISTRATION"
            skipped.append(mag_id)

        # derived_from: use source secondary sample accession
        derived_from = r.get("source_secondary_sample_accession","").strip() or \
                       r.get("source_primary_sample_accession","").strip()

        completeness   = r.get("original_checkm2_completeness","").strip() or r.get("final_css","")
        contamination  = r.get("original_checkm2_contamination","").strip() or r.get("final_rrs","")
        mimags_quality = r.get("mimags_quality","").strip() or "Medium-quality draft"

        writer.writerow({
            "sample_alias":         f"MAG_{mag_id}",
            "tax_id":               tax_id,
            "scientific_name":      sci_name,
            "sample_title":         f"MAG {mag_id} from {dataset} livestock metagenome",
            "sample_description":   (
                f"Metagenome-assembled genome (MAG) recovered from {ctx['isolation_source']}. "
                f"GTDB taxonomy: {r.get('gtdb_taxonomy','unclassified')}. "
                f"CheckM2 completeness: {completeness}%, contamination: {contamination}%. "
                f"GUNC passed: {r.get('final_gunc_status','')}. "
                f"MIMAGs quality: {mimags_quality}."
            ),
            "collection date":      ctx["collection_date"],
            "geographic location (country and/or sea)": ctx["country"],
            "geographic location (latitude)":           ctx["lat"],
            "geographic location (longitude)":          ctx["lon"],
            "isolation_source":     ctx["isolation_source"],
            "broad-scale environmental context": ctx["biome"],
            "local environmental context":       ctx["feature"],
            "environmental medium":              ctx["material"],
            "completeness score":   completeness,
            "contamination score":  contamination,
            "assembly quality":     mimags_quality,
            "assembly software":    r.get("assembly_program","metaWRAP; MAGpurify2; DeepPurify"),
            "sequencing method":    "Illumina",
            "investigation type":   "metagenome-assembled genome",
            "project name":         ctx["study"],
            "derived from":         derived_from,
            "binning software":     "metaWRAP (MaxBin2, MetaBAT2)",
            "taxonomic identity marker": "GTDB-Tk marker genes (bac120/ar53)",
            "taxonomy reference database name":    "GTDB",
            "taxonomy reference database version": "r220",
        })

print(f"  -> {bs_file.name} ({len(rows)} amostras)")
if skipped:
    print(f"  ATENCAO: {len(skipped)} MAGs com PENDING_ENA_REGISTRATION — submeter DEPOIS do registo:")
    for m in skipped:
        print(f"    {m}")

# ── FILE 3: Assembly manifests (webin-cli) ──────────────────────────────────
print("\n[3/3] Gerando manifests webin-cli...")

# Also generate a batch submit script
manifest_count = 0
missing_fasta  = []

for r in rows:
    mag_id  = r["mag_id"]
    dataset = r["dataset"]
    fasta   = local_fasta(mag_id, dataset)

    if not fasta:
        missing_fasta.append(mag_id)
        continue

    study = r.get("study_accession","").strip()
    if not study or study.startswith("TODO"):
        study = CONTEXT[dataset]["study"]

    # SAMPLE: will be filled after BioSample registration — use alias for now
    sample_alias = f"MAG_{mag_id}"

    # RUN_REF: from source run accession
    run_ref = r.get("source_run_accession","").strip() or \
              r.get("run_accession","").strip()

    # Coverage
    cov = r.get("mean_coverage","").strip() or "1"
    try:
        cov = f"{float(cov):.1f}"
    except ValueError:
        cov = "1"

    # Assembly name (no spaces, ≤50 chars)
    asm_name = mag_id.replace("__", "_").replace(" ", "_")[:50]

    manifest_path = MANIFEST_DIR / f"{mag_id}_manifest.txt"
    with open(manifest_path, "w") as mf:
        mf.write(f"STUDY\t{study}\n")
        mf.write(f"SAMPLE\t{sample_alias}\n")   # UPDATE after BioSample registration
        mf.write(f"ASSEMBLYNAME\t{asm_name}\n")
        mf.write(f"ASSEMBLY_TYPE\tMetagenome-Assembled Genome (MAG)\n")
        mf.write(f"COVERAGE\t{cov}\n")
        mf.write(f"PROGRAM\t{r.get('assembly_program','metaWRAP; MAGpurify2; DeepPurify')}\n")
        mf.write(f"PLATFORM\tILLUMINA\n")
        mf.write(f"MOLECULETYPE\tgenomic DNA\n")
        if run_ref and not run_ref.startswith("TODO"):
            mf.write(f"RUN_REF\t{run_ref}\n")
        mf.write(f"FASTA\t{fasta}\n")
    manifest_count += 1

print(f"  -> manifests/  ({manifest_count} ficheiros)")
if missing_fasta:
    print(f"  SEM FASTA ({len(missing_fasta)}): {missing_fasta[:5]}")

# ── Batch submit script ──────────────────────────────────────────────────────
submit_script = OUT / "submit_all_assemblies.sh"
webin_jar = BASE.parent.parent / "ENA_submission/webin-cli.jar"

with open(submit_script, "w") as sh:
    sh.write("#!/bin/bash\n")
    sh.write("# Submete os 729 MAGs via webin-cli\n")
    sh.write("# ANTES de correr: actualizar campo SAMPLE nos manifests com SAMEA accessions\n")
    sh.write("# Uso: bash submit_all_assemblies.sh 'senha_webin'\n\n")
    sh.write("PASSWORD=${1:?'Fornecer senha: bash submit_all_assemblies.sh SENHA'}\n")
    sh.write("WEBIN_JAR=/path/to/webin-cli.jar\n")
    sh.write("USERNAME=Webin-XXXXX  # substituir pelo username ENA\n\n")
    sh.write(f"MANIFEST_DIR={MANIFEST_DIR}\n\n")
    sh.write("ok=0; fail=0\n")
    sh.write('for manifest in "$MANIFEST_DIR"/*_manifest.txt; do\n')
    sh.write('    mag=$(basename "$manifest" _manifest.txt)\n')
    sh.write('    echo "Submetendo $mag..."\n')
    sh.write('    java -jar "$WEBIN_JAR" \\\n')
    sh.write('        -context genome \\\n')
    sh.write('        -manifest "$manifest" \\\n')
    sh.write('        -username "$USERNAME" \\\n')
    sh.write('        -password "$PASSWORD" \\\n')
    sh.write('        -submit \\\n')
    sh.write('        2>&1 | tee -a submit_log_${mag}.txt\n')
    sh.write('    if [ $? -eq 0 ]; then ((ok++)); else ((fail++)); fi\n')
    sh.write('done\n')
    sh.write('echo "Concluido: $ok OK, $fail erros"\n')
os.chmod(submit_script, 0o755)

# ── Summary ──────────────────────────────────────────────────────────────────
print(f"""
{'='*60}
CONCLUÍDO — ficheiros em: {OUT}

  ena_taxonomy_registration_request.tsv   ({len(by_name)} nomes)
  ena_biosample_sheet_729.tsv             ({len(rows)} amostras)
  manifests/                              ({manifest_count} manifests)
  submit_all_assemblies.sh

PRÓXIMOS PASSOS:
  1. Registar os {len(by_name)} nomes no ENA Taxonomy Request
     https://www.ebi.ac.uk/ena/browser/support
     -> copiar conteúdo de ena_taxonomy_registration_request.tsv
     -> aguardar resposta com taxids (dias)

  2. Actualizar ena_biosample_sheet_729.tsv:
     -> substituir PENDING_ENA_REGISTRATION pelos taxids recebidos
     -> nas {len(skipped)} linhas: {[r for r in rows if r['tax_id_status']=='needs_ena_registration']}

  3. Registar BioSamples:
     -> Webin Portal -> Submit -> Samples -> Upload TSV
     -> Receber SAMEA accessions para cada MAG

  4. Actualizar campo SAMPLE nos manifests com os SAMEA:
     python3 update_manifests_samea.py samea_accessions.tsv

  5. Submeter assemblies:
     bash submit_all_assemblies.sh 'senha'
{'='*60}
""")
