#!/usr/bin/env python3
"""
Gera tabelas de submissão ENA finais para 729 MAGs bacterianos.

Operações:
  1. Remove 9 bins sem classificação GTDB (Unclassified/MSA insuficiente)
  2. Integra taxonomia Galaxy215 nos 41 bins que agora têm classificação
  3. Corrige caminhos de FASTAs (/media -> /mnt/nvme2/... onde disponíveis)
  4. Gera:
       ena_mag_candidates_final_729.tsv     — tabela master com 729 MAGs
       ena_mag_sample_sheet_729.tsv         — sample sheet para registo ENA
       ena_mag_assembly_manifests_729.tsv   — manifests para webin-cli
       ena_mag_fasta_availability_729.tsv   — auditoria de disponibilidade de FASTAs

Uso: python3 build_final_submission_729.py
"""

import csv
import re
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
# USER: set your base path here
BASE_SUBMISSION   = Path("/mnt/nvme2/Metagenome_Spain_Miuvigs/ENA_MAGS_Submission")
# USER: set your base path here
BASE_RESULTS      = Path("/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Mags_Vmags")
OUT_DIR           = BASE_RESULTS / "ENA_final_submission_729"

INPUT_TABLE       = BASE_SUBMISSION / "ena_mag_submit_candidates_finalized.tsv"
GALAXY215         = BASE_RESULTS   / "Galaxy215-[gtdbtk.bac120.summary].tsv"
RESCUED_FASTAS_BR = BASE_SUBMISSION / "rescued_fastas" / "Brazil"
RESCUED_FASTAS_SP = BASE_SUBMISSION / "rescued_fastas" / "Spain"
GALAXY_RERUN_SP   = BASE_SUBMISSION / "galaxy_gtdb_rerun_spain"
GALAXY_RERUN_BR   = BASE_SUBMISSION / "galaxy_gtdb_rerun_brazil"

# ── 9 bins a excluir (sem classificação GTDB no Galaxy215) ────────────────────
EXCLUIR = {
    "Brazil__P19109_1048__bin.3",
    "Brazil__P19109_1048__bin.5",
    "Brazil__P19109_1052__bin.3",
    "Brazil__P19109_1052__bin.4",
    "Brazil__P19109_1052__bin.5",
    "Brazil__P19109_1052__bin.6",
    "Spain__SRR11615166__bin_45",
    "Spain__SRR11615255__bin_27",
    "Spain__SRR11615256__bin_42_sub",
}

SUFFIX_PATTERN = re.compile(r"_[A-Z]$")


# ── GTDB helpers ───────────────────────────────────────────────────────────────
def parse_gtdb_classification(classification: str) -> dict:
    """Devolve dict {d, p, c, o, f, g, s} da string GTDB."""
    result = {}
    for token in classification.split(";"):
        token = token.strip()
        if "__" in token:
            rank, value = token.split("__", 1)
            result[rank.strip()] = value.strip()
    return result


def gtdb_species_name(parsed: dict) -> str:
    """Nome de espécie (ou fallback) do GTDB para submeter ao ENA."""
    s = parsed.get("s", "").strip()
    if s and s not in ("", " "):
        return SUFFIX_PATTERN.sub("", s)
    g = parsed.get("g", "").strip()
    if g:
        return SUFFIX_PATTERN.sub("", g) + " sp."
    f = parsed.get("f", "").strip()
    if f:
        return SUFFIX_PATTERN.sub("", f) + " bacterium"
    o = parsed.get("o", "").strip()
    if o:
        return SUFFIX_PATTERN.sub("", o) + " bacterium"
    p = parsed.get("p", "").strip()
    if p:
        return SUFFIX_PATTERN.sub("", p) + " bacterium"
    return "bacterium"


def gtdb_scientific_name_status(parsed: dict) -> str:
    if parsed.get("s", "").strip():
        return "non_binomial_gtdb" if re.search(r" sp\d", parsed["s"]) else "exact_species_gtdb"
    if parsed.get("g", ""):
        return "fallback_genus_sp"
    return "fallback_higher_taxon"


# ── FASTA resolver ─────────────────────────────────────────────────────────────
def resolve_fasta(row: dict, is_galaxy215_bin: bool) -> tuple[str, str]:
    """
    Devolve (fasta_path, fasta_source).
    fasta_source: 'rescued_fastas' | 'galaxy_rerun' | 'original_offline' | 'not_found'
    """
    mag_id   = row["mag_id"]
    dataset  = row["dataset"]
    rescued  = row.get("rescued_after_cleaning", "") == "yes"

    # 1. Bins do Galaxy215 não-rescued: copias em galaxy_rerun_dirs
    if is_galaxy215_bin and not rescued:
        d = GALAXY_RERUN_SP if dataset == "Spain" else GALAXY_RERUN_BR
        p = d / f"{mag_id}.fa"
        if p.exists():
            return str(p), "galaxy_rerun"

    # 2. Rescued bins: copias em rescued_fastas/
    if rescued:
        d = RESCUED_FASTAS_SP if dataset == "Spain" else RESCUED_FASTAS_BR
        p = d / f"{mag_id}.final.fa"
        if p.exists():
            return str(p), "rescued_fastas"

    # 3. Galaxy215 rescued: também devem estar em rescued_fastas
    if is_galaxy215_bin and rescued:
        d = RESCUED_FASTAS_SP if dataset == "Spain" else RESCUED_FASTAS_BR
        p = d / f"{mag_id}.final.fa"
        if p.exists():
            return str(p), "rescued_fastas"

    # 4. Fallback: path original (provavelmente offline)
    original = row.get("final_fasta", "")
    if original:
        return original, "original_offline"

    return "", "not_found"


# ── Main ───────────────────────────────────────────────────────────────────────
def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Carregar tabela principal
    with INPUT_TABLE.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    print(f"Tabela principal: {len(rows)} MAGs")

    # 2. Carregar Galaxy215 GTDB results
    galaxy215: dict[str, dict] = {}
    with GALAXY215.open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            mag_id = row["user_genome"].removesuffix(".fa")
            galaxy215[mag_id] = row
    print(f"Galaxy215: {len(galaxy215)} entradas (51 bins)")

    # 3. Filtrar: excluir 9 bins sem classificação
    rows = [r for r in rows if r["mag_id"] not in EXCLUIR]
    print(f"Após exclusão dos 9 sem classificação: {len(rows)} MAGs")
    br = sum(1 for r in rows if r["dataset"] == "Brazil")
    sp = sum(1 for r in rows if r["dataset"] == "Spain")
    print(f"  Brazil: {br}  Spain: {sp}")

    # 4. Integrar taxonomia Galaxy215 nos 41 bins classificados
    galaxy215_classified = {k: v for k, v in galaxy215.items() if k not in EXCLUIR}
    updated_tax = 0
    for row in rows:
        mid = row["mag_id"]
        if mid not in galaxy215_classified:
            continue
        g215 = galaxy215_classified[mid]
        classification = g215.get("classification", "")
        if not classification or "Unclassified" in classification:
            continue
        parsed = parse_gtdb_classification(classification)
        species_name = gtdb_species_name(parsed)
        status = gtdb_scientific_name_status(parsed)

        # Actualizar campos de taxonomia na linha
        row["gtdb_taxonomy"]         = classification
        row["gtdb_species"]          = parsed.get("s", "")
        row["gtdb_genus"]            = parsed.get("g", "")
        row["gtdb_domain"]           = "Bacteria"
        row["taxonomy_source"]       = "Galaxy215_GTDB"
        row["taxon_name_for_submission"] = species_name
        row["scientific_name_status"] = status
        row["scientific_name_for_ena"] = species_name
        # Taxid deixa-se "unresolved" para o utilizador registar via ENA portal
        # (ENA para MIMAGS non-binomial: registar como environmental name)
        row["tax_id_status"] = "needs_ena_registration" if "sp." in species_name or "bacterium" in species_name else "unresolved"
        # gtdb_classification_method e closest_reference do Galaxy215
        row["gtdb_classification_method"] = g215.get("classification_method", "")
        row["gtdb_closest_reference"]     = g215.get("closest_genome_reference", "")
        row["gtdb_closest_taxonomy"]      = g215.get("closest_genome_taxonomy", "")
        row["gtdb_closest_ani"]           = g215.get("closest_genome_ani", "")
        row["gtdb_closest_af"]            = g215.get("closest_genome_af", "")
        updated_tax += 1
    print(f"Taxonomia Galaxy215 integrada em {updated_tax} bins")

    # 5. Resolver FASTAs e corrigir caminhos
    fasta_audit = []
    for row in rows:
        mid = row["mag_id"]
        is_g215 = mid in galaxy215_classified
        fasta_path, fasta_source = resolve_fasta(row, is_g215)
        row["final_fasta"]        = fasta_path
        row["fasta_source"]       = fasta_source
        fasta_audit.append({
            "mag_id":       mid,
            "dataset":      row["dataset"],
            "fasta_source": fasta_source,
            "fasta_path":   fasta_path,
            "available":    "yes" if "offline" not in fasta_source and fasta_path else "no",
        })

    fasta_available   = sum(1 for a in fasta_audit if a["available"] == "yes")
    fasta_offline     = sum(1 for a in fasta_audit if a["fasta_source"] == "original_offline")
    fasta_not_found   = sum(1 for a in fasta_audit if a["fasta_source"] == "not_found")
    print(f"FASTAs disponíveis: {fasta_available} | offline (disco externo): {fasta_offline} | não encontrado: {fasta_not_found}")

    # 6. Actualizar submission_missing_reasons_final
    for row in rows:
        missing = [m for m in (row.get("submission_missing_reasons_final") or "").split(";") if m]
        # Remove missing_tax_id se agora tem taxon_name (mesmo que ainda precise de registo ENA)
        if row.get("taxon_name_for_submission", "").strip():
            missing = [m for m in missing if m not in ("missing_tax_id", "missing_binomial_scientific_name")]
        # Marca FASTA offline
        if row.get("fasta_source", "") == "original_offline":
            if "fasta_offline" not in missing:
                missing.append("fasta_offline")
        row["submission_missing_reasons_final"] = ";".join(dict.fromkeys(missing))

    # ── OUTPUT 1: tabela master ────────────────────────────────────────────────
    out_candidates = OUT_DIR / "ena_mag_candidates_final_729.tsv"
    fieldnames = list(rows[0].keys())
    # Garantir campos novos estão incluídos
    for f in ("fasta_source",):
        if f not in fieldnames:
            fieldnames.append(f)
    with out_candidates.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nOutput 1: {out_candidates}")

    # ── OUTPUT 2: sample sheet para registo ENA BioSamples ────────────────────
    sample_sheet_fields = [
        "alias", "sample_type", "dataset", "mag_id", "sample_id",
        "checklist", "scientific_name", "tax_id", "tax_id_status",
        "gtdb_taxonomy", "gtdb_classification_method",
        "gtdb_closest_reference", "gtdb_closest_taxonomy",
        "sample_derived_from", "assembly_quality",
        "completeness_score", "contamination_score",
        "gunc_final_status", "mean_coverage", "covered_fraction",
        "final_fasta", "fasta_source", "study_accession", "notes",
    ]
    out_sample_sheet = OUT_DIR / "ena_mag_sample_sheet_729.tsv"
    with out_sample_sheet.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=sample_sheet_fields)
        w.writeheader()
        for row in rows:
            comp = row.get("final_checkm2_completeness", "")
            cont = row.get("final_checkm2_contamination", "")
            try:
                c, x = float(comp), float(cont)
                if c >= 90 and x <= 5:
                    quality = "High-quality draft"
                elif c >= 50 and x <= 10:
                    quality = "Medium-quality draft"
                else:
                    quality = "Low-quality draft"
            except ValueError:
                quality = "Medium-quality draft"

            w.writerow({
                "alias":                    f"MAG_{row['mag_id']}",
                "sample_type":              "MAG",
                "dataset":                  row["dataset"],
                "mag_id":                   row["mag_id"],
                "sample_id":                row["sample_id"],
                "checklist":                "GSC MIMAGS",
                "scientific_name":          row.get("scientific_name_for_ena", ""),
                "tax_id":                   row.get("tax_id", ""),
                "tax_id_status":            row.get("tax_id_status", ""),
                "gtdb_taxonomy":            row.get("gtdb_taxonomy", ""),
                "gtdb_classification_method": row.get("gtdb_classification_method", ""),
                "gtdb_closest_reference":   row.get("gtdb_closest_reference", ""),
                "gtdb_closest_taxonomy":    row.get("gtdb_closest_taxonomy", ""),
                "sample_derived_from":      row.get("source_primary_sample_accession", ""),
                "assembly_quality":         quality,
                "completeness_score":       comp,
                "contamination_score":      cont,
                "gunc_final_status":        row.get("final_gunc_status", ""),
                "mean_coverage":            row.get("mean_coverage", ""),
                "covered_fraction":         row.get("covered_fraction", ""),
                "final_fasta":              row.get("final_fasta", ""),
                "fasta_source":             row.get("fasta_source", ""),
                "study_accession":          row.get("study_accession", ""),
                "notes":                    row.get("submission_missing_reasons_final", ""),
            })
    print(f"Output 2: {out_sample_sheet}")

    # ── OUTPUT 3: manifest table para webin-cli ────────────────────────────────
    out_manifests = OUT_DIR / "ena_mag_assembly_manifests_729.tsv"
    manifest_fields = [
        "mag_id", "dataset", "sample_id",
        "study_accession", "sample_alias",
        "run_ref", "assembly_name", "assembly_type",
        "coverage", "program", "platform",
        "molecule_type", "description", "fasta", "fasta_available",
        "manifest_file",
    ]
    with out_manifests.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=manifest_fields)
        w.writeheader()
        for row in rows:
            fasta_avail = "yes" if row.get("fasta_source","") not in ("original_offline","not_found") else "no"
            w.writerow({
                "mag_id":          row["mag_id"],
                "dataset":         row["dataset"],
                "sample_id":       row["sample_id"],
                "study_accession": row.get("study_accession", ""),
                "sample_alias":    f"MAG_{row['mag_id']}",
                "run_ref":         row.get("source_run_accession", ""),
                "assembly_name":   row["mag_id"],
                "assembly_type":   "Metagenome-Assembled Genome (MAG)",
                "coverage":        row.get("mean_coverage", ""),
                "program":         row.get("assembly_program", ""),
                "platform":        row.get("platform", ""),
                "molecule_type":   "genomic DNA",
                "description":     row.get("assembly_description", f"MAG from {row['dataset']} metagenome"),
                "fasta":           row.get("final_fasta", ""),
                "fasta_available": fasta_avail,
                "manifest_file":   f"{row['mag_id']}.manifest",
            })
    print(f"Output 3: {out_manifests}")

    # ── OUTPUT 4: auditoria de FASTAs ─────────────────────────────────────────
    out_fasta_audit = OUT_DIR / "ena_mag_fasta_availability_729.tsv"
    with out_fasta_audit.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(fasta_audit[0].keys()))
        w.writeheader()
        w.writerows(fasta_audit)
    print(f"Output 4: {out_fasta_audit}")

    # ── Sumário final ─────────────────────────────────────────────────────────
    total = len(rows)
    br_total = sum(1 for r in rows if r["dataset"] == "Brazil")
    sp_total = sum(1 for r in rows if r["dataset"] == "Spain")

    tax_resolved  = sum(1 for r in rows if r.get("tax_id_status","") == "resolved")
    tax_unresolved= sum(1 for r in rows if r.get("tax_id_status","") in ("unresolved","needs_ena_registration","taxonomy_missing"))
    tax_g215      = sum(1 for r in rows if r.get("taxonomy_source","") == "Galaxy215_GTDB")

    missing_final = sum(1 for r in rows if r.get("submission_missing_reasons_final","") and
                        r["submission_missing_reasons_final"] not in ("","missing_derived_sample_accession;missing_run_accession;missing_source_run_accession"))

    print(f"""
╔══════════════════════════════════════════════════════════╗
║         SUBMISSÃO ENA — {total} MAGs BACTERIANOS         ║
╠══════════════════════════════════════════════════════════╣
║  Brazil: {br_total:<5}   Spain: {sp_total:<5}   Total: {total:<5}             ║
╠══════════════════════════════════════════════════════════╣
║  TAXONOMIA                                               ║
║   Tax_id resolvido (NCBI): {tax_resolved:<5}                        ║
║   Taxonomia Galaxy215 nova: {tax_g215:<5}                       ║
║   Sem tax_id (precisam registo ENA): {tax_unresolved:<5}              ║
╠══════════════════════════════════════════════════════════╣
║  FASTAs disponíveis: {fasta_available:<5}                              ║
║  FASTAs offline (disco externo): {fasta_offline:<5}                  ║
╠══════════════════════════════════════════════════════════╣
║  Outputs:                                                ║
║   {str(OUT_DIR):<54}  ║
╚══════════════════════════════════════════════════════════╝
""")

    # ── Sumário de bins Galaxy215 integrados ──────────────────────────────────
    g215_rows = [r for r in rows if r.get("taxonomy_source","") == "Galaxy215_GTDB"]
    print("Bins com taxonomia Galaxy215 integrada:")
    for r in sorted(g215_rows, key=lambda x: x["mag_id"]):
        print(f"  {r['mag_id']}")
        print(f"    GTDB: {r['gtdb_taxonomy'][:80]}")
        print(f"    ENA name: {r['scientific_name_for_ena']}  | tax_id_status: {r['tax_id_status']}")

    # ── Bins que precisam de registo de taxon no ENA ──────────────────────────
    needs_reg = [r for r in rows if r.get("tax_id_status","") == "needs_ena_registration"]
    if needs_reg:
        print(f"\nBins que precisam de registo de taxon no ENA ({len(needs_reg)}):")
        for r in needs_reg:
            print(f"  {r['mag_id']}: {r.get('scientific_name_for_ena','')}")


if __name__ == "__main__":
    main()
