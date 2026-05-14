#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from pathlib import Path

import pandas as pd


# USER: set your base path here


ROOT = Path("/mnt/nvme2/Metagenome_Spain_Miuvigs/Refresh_Genomad")
MASTER_QC = ROOT / "Spain_vRhyme_phage_qc_master.tsv"

FINAL_ALL = ROOT / "Spain_vRhyme_final_taxonomy_unified.tsv"
FINAL_CONFIRMED = ROOT / "Spain_vRhyme_final_taxonomy_lifestyle_confirmed.tsv"
FINAL_TEMPERATE = ROOT / "Spain_vRhyme_final_taxonomy_temperate.tsv"
FINAL_VIRULENT = ROOT / "Spain_vRhyme_final_taxonomy_virulent.tsv"
FINAL_LIFESTYLE_NA = ROOT / "Spain_vRhyme_final_taxonomy_lifestyle_na.tsv"
FINAL_EXCLUDED = ROOT / "Spain_vRhyme_final_taxonomy_excluded_audit.tsv"
FINAL_XLSX = ROOT / "Spain_vRhyme_final_taxonomy_unified.xlsx"

AMBIGUOUS = {
    "",
    "-",
    "unknown",
    "unclassified",
    "viruses",
    "no hits to database",
    "hits not found in taxonomy files",
    "filtered",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def normalize(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[\s\-]+", "_", value)
    return value


def parse_lineage(lineage: str) -> list[str]:
    if normalize(lineage) in AMBIGUOUS:
        return []
    tokens: list[str] = []
    for raw in lineage.split(";"):
        part = raw.strip()
        if not part or part == "-":
            continue
        if ":" in part:
            part = part.split(":", 1)[1].strip()
        if not part:
            continue
        if normalize(part) in AMBIGUOUS:
            continue
        tokens.append(part)
    return tokens


def common_prefix(left: list[str], right: list[str]) -> list[str]:
    prefix: list[str] = []
    for lval, rval in zip(left, right):
        if normalize(lval) != normalize(rval):
            break
        prefix.append(lval)
    return prefix


def choose_unified_taxonomy(row: dict[str, str]) -> tuple[str, str, str]:
    taxonomy_class = row.get("taxonomy_class", "")
    if taxonomy_class in {"taxonomy_conflict", "non_phage_excluded"}:
        return "", "excluded", ""

    genomad = parse_lineage(row.get("consensus_taxonomy", ""))
    phabox = parse_lineage(row.get("phabox2_taxonomy", ""))
    prefix = common_prefix(genomad, phabox)

    if phabox and genomad:
        if len(prefix) == min(len(genomad), len(phabox)):
            chosen = phabox if len(phabox) >= len(genomad) else genomad
            source = "phabox2" if len(phabox) >= len(genomad) else "genomad"
            return ";".join(chosen), source, ";".join(prefix)
        if prefix:
            chosen = phabox if len(phabox) >= len(genomad) else genomad
            source = "phabox2_partial" if len(phabox) >= len(genomad) else "genomad_partial"
            return ";".join(chosen), source, ";".join(prefix)
        return ";".join(phabox), "phabox2", ""
    if phabox:
        return ";".join(phabox), "phabox2", ""
    if genomad:
        return ";".join(genomad), "genomad", ""
    return "", "unclassified", ""


def remap_lifestyle_raw(mapped: str) -> str:
    if mapped == "lytic":
        return "virulent"
    if mapped == "temperate":
        return "temperate"
    return "-"


def article_lifestyle(raw: str) -> str:
    if raw == "virulent":
        return "Virulent"
    if raw == "temperate":
        return "Temperate"
    return "-"


def main() -> int:
    rows = read_tsv(MASTER_QC)
    final_rows: list[dict[str, str]] = []
    for row in rows:
        unified_taxonomy, taxonomy_source, agreed_prefix = choose_unified_taxonomy(row)
        merged = dict(row)
        merged["phabox2_lifestyle_raw"] = remap_lifestyle_raw(row.get("phabox2_lifestyle", ""))
        merged["Lifestyle"] = article_lifestyle(merged["phabox2_lifestyle_raw"])
        merged["unified_taxonomy"] = unified_taxonomy
        merged["unified_taxonomy_source"] = taxonomy_source
        merged["agreed_taxonomy_prefix"] = agreed_prefix
        final_rows.append(merged)

    fieldnames = list(final_rows[0].keys()) if final_rows else []
    phage_rows = [row for row in final_rows if row.get("final_keep_phage") == "yes"]
    confirmed_rows = [row for row in phage_rows if row.get("Lifestyle") in {"Temperate", "Virulent"}]
    temperate_rows = [row for row in confirmed_rows if row.get("Lifestyle") == "Temperate"]
    virulent_rows = [row for row in confirmed_rows if row.get("Lifestyle") == "Virulent"]
    lifestyle_na_rows = [row for row in phage_rows if row.get("Lifestyle") == "-"]
    excluded_rows = [row for row in final_rows if row.get("final_bucket") == "excluded"]

    write_tsv(FINAL_ALL, phage_rows, fieldnames)
    write_tsv(FINAL_CONFIRMED, confirmed_rows, fieldnames)
    write_tsv(FINAL_TEMPERATE, temperate_rows, fieldnames)
    write_tsv(FINAL_VIRULENT, virulent_rows, fieldnames)
    write_tsv(FINAL_LIFESTYLE_NA, lifestyle_na_rows, fieldnames)
    write_tsv(FINAL_EXCLUDED, excluded_rows, fieldnames)

    with pd.ExcelWriter(FINAL_XLSX, engine="openpyxl") as writer:
        pd.DataFrame(phage_rows).to_excel(writer, sheet_name="all_phage", index=False)
        pd.DataFrame(confirmed_rows).to_excel(writer, sheet_name="confirmed_lifestyle", index=False)
        pd.DataFrame(temperate_rows).to_excel(writer, sheet_name="temperate", index=False)
        pd.DataFrame(virulent_rows).to_excel(writer, sheet_name="virulent", index=False)
        pd.DataFrame(lifestyle_na_rows).to_excel(writer, sheet_name="lifestyle_na", index=False)
        pd.DataFrame(excluded_rows).to_excel(writer, sheet_name="excluded_audit", index=False)

    print(f"all_phage\t{FINAL_ALL}")
    print(f"confirmed\t{FINAL_CONFIRMED}")
    print(f"temperate\t{FINAL_TEMPERATE}")
    print(f"virulent\t{FINAL_VIRULENT}")
    print(f"lifestyle_na\t{FINAL_LIFESTYLE_NA}")
    print(f"excluded\t{FINAL_EXCLUDED}")
    print(f"xlsx\t{FINAL_XLSX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
