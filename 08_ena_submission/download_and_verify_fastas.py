#!/usr/bin/env python3
"""
Download all 729 MAG FASTAs from Google Drive (rclone remote: ubu:)
and verify taxonomy coverage.

Sources:
  1. ubu:home/Tati_Metagenome_Mags_Filtered/Brazil/rescued   -> 113 Brazil rescued
  2. ubu:home/Tati_Metagenome_Mags_Filtered/Spain/rescued    ->  59 Spain rescued
  3. ubu:home/MAGS_BINS_2  (metawrap_50_10_bins/)            -> 444 Brazil initial_gunc_true
  4. ubu:home/Metagenome_Spain/MAGs  (dastool/*_DASTOOL_bins/) -> 113 Spain initial_gunc_true

Output: ENA_final_submission_729/fastas/Brazil/<mag_id>.fa
        ENA_final_submission_729/fastas/Spain/<mag_id>.fa
"""

import subprocess
import csv
import sys
import re
from pathlib import Path

# USER: set your base path here

BASE     = Path("/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Mags_Vmags")
OUT_DIR  = BASE / "ENA_final_submission_729/fastas"
BRAZIL   = OUT_DIR / "Brazil"
SPAIN    = OUT_DIR / "Spain"
CANDIDATES = BASE / "ENA_final_submission_729/ena_mag_candidates_final_729.tsv"

BRAZIL.mkdir(parents=True, exist_ok=True)
SPAIN.mkdir(parents=True, exist_ok=True)

# ── Load 729 MAG list with taxonomy ────────────────────────────────────────
print("Carregando tabela de candidatos...")
rows = {}
with open(CANDIDATES) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows[row["mag_id"]] = row

print(f"  {len(rows)} MAGs na tabela")

# ── Define download jobs ────────────────────────────────────────────────────
# Each job: (mag_id, rclone_src, dest_path)
jobs = []
missing_taxonomy = []

def rclone_ls(remote_path, include_pattern=None):
    cmd = ["rclone", "ls", remote_path]
    if include_pattern:
        cmd += ["--include", include_pattern]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout.strip().split("\n") if result.stdout.strip() else []

# ── Source 1: Brazil rescued ────────────────────────────────────────────────
print("\n[1/4] Listando Brazil/rescued no Drive...")
src1 = "ubu:home/Tati_Metagenome_Mags_Filtered/Brazil/rescued"
lines = rclone_ls(src1)
for line in lines:
    parts = line.strip().split(None, 1)
    if len(parts) < 2:
        continue
    fname = parts[1]  # e.g. Brazil__P19109_1035__bin.1.final.fa
    mag_id = fname.removesuffix(".final.fa").removesuffix(".fa")
    if mag_id in rows:
        dest = BRAZIL / f"{mag_id}.fa"
        jobs.append((mag_id, f"{src1}/{fname}", dest))

# ── Source 2: Spain rescued ─────────────────────────────────────────────────
print("[2/4] Listando Spain/rescued no Drive...")
src2 = "ubu:home/Tati_Metagenome_Mags_Filtered/Spain/rescued"
lines = rclone_ls(src2)
for line in lines:
    parts = line.strip().split(None, 1)
    if len(parts) < 2:
        continue
    fname = parts[1]
    mag_id = fname.removesuffix(".final.fa").removesuffix(".fa")
    if mag_id in rows:
        dest = SPAIN / f"{mag_id}.fa"
        jobs.append((mag_id, f"{src2}/{fname}", dest))

# ── Source 3: Brazil MAGS_BINS_2 (metawrap) ────────────────────────────────
print("[3/4] Listando MAGS_BINS_2 no Drive...")
src3 = "ubu:home/MAGS_BINS_2"
lines = rclone_ls(src3, "*/metawrap_50_10_bins/*.fa")
for line in lines:
    parts = line.strip().split(None, 1)
    if len(parts) < 2:
        continue
    fpath = parts[1]  # e.g. P19109_1042/metawrap_50_10_bins/bin.10.fa
    m = re.match(r'^(P19109_\d+)/metawrap_50_10_bins/bin\.(\d+)\.fa$', fpath)
    if not m:
        continue
    sample, bnum = m.group(1), m.group(2)
    mag_id = f"Brazil__{sample}__bin.{bnum}"
    if mag_id in rows:
        dest = BRAZIL / f"{mag_id}.fa"
        if not any(j[0] == mag_id for j in jobs):  # skip if already from rescued
            jobs.append((mag_id, f"{src3}/{fpath}", dest))

# ── Source 4: Spain Metagenome_Spain/MAGs (dastool) ────────────────────────
print("[4/4] Listando Metagenome_Spain/MAGs no Drive...")
src4 = "ubu:home/Metagenome_Spain/MAGs"
lines = rclone_ls(src4, "*/dastool/*_DASTOOL_bins/*.fa")
for line in lines:
    parts = line.strip().split(None, 1)
    if len(parts) < 2:
        continue
    fpath = parts[1]  # e.g. SRR11615157/dastool/SRR11615157_DASTOOL_bins/bin_29.fa
    m = re.match(r'^(SRR\d+)/dastool/\S+_DASTOOL_bins/bin_(.+)\.fa$', fpath)
    if not m:
        continue
    srr, bname = m.group(1), m.group(2)
    mag_id = f"Spain__{srr}__bin_{bname}"
    if mag_id in rows:
        dest = SPAIN / f"{mag_id}.fa"
        if not any(j[0] == mag_id for j in jobs):
            jobs.append((mag_id, f"{src4}/{fpath}", dest))

# ── Summary before download ─────────────────────────────────────────────────
covered = {j[0] for j in jobs}
not_covered = [m for m in rows if m not in covered]

print(f"\n{'='*60}")
print(f"Jobs de download: {len(jobs)}")
print(f"MAGs cobertos:    {len(covered)} / {len(rows)}")
if not_covered:
    print(f"NAO ENCONTRADOS ({len(not_covered)}):")
    for m in sorted(not_covered):
        print(f"  {m}")

# ── Taxonomy check ──────────────────────────────────────────────────────────
no_tax = []
for mag_id, row in rows.items():
    tax    = row.get("tax_id", "").strip()
    status = row.get("tax_id_status", "").strip()
    org    = row.get("scientific_name_for_ena", "").strip()
    if not tax or not org or status not in ("resolved", "needs_ena_registration"):
        no_tax.append(mag_id)

print(f"\nMAGs SEM taxonomia completa: {len(no_tax)}")
if no_tax:
    for m in sorted(no_tax)[:20]:
        r = rows[m]
        print(f"  {m}: tax_id='{r.get('tax_id','')}' status='{r.get('tax_id_status','')}' org='{r.get('scientific_name_for_ena','')}'")

if not_covered or no_tax:
    resp = input(f"\nContinuar com o download dos {len(jobs)} FASTAs disponíveis? [s/N] ").strip().lower()
    if resp != 's':
        print("Abortado.")
        sys.exit(0)
else:
    print("\nTudo OK — iniciando download...")

# ── Download ────────────────────────────────────────────────────────────────
print(f"\nIniciando download de {len(jobs)} FASTAs...")
downloaded = 0
errors = 0

for i, (mag_id, src, dest) in enumerate(jobs):
    if dest.exists():
        downloaded += 1
        if i % 100 == 0:
            print(f"  [{i+1}/{len(jobs)}] JÁ EXISTE: {mag_id}")
        continue
    result = subprocess.run(
        ["rclone", "copyto", src, str(dest)],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        downloaded += 1
        if i % 50 == 0:
            print(f"  [{i+1}/{len(jobs)}] OK: {mag_id}")
    else:
        errors += 1
        print(f"  ERRO [{i+1}]: {mag_id} — {result.stderr[:120]}")
        if errors > 20:
            print("Muitos erros, parando.")
            break

print(f"\n{'='*60}")
print(f"Download concluído: {downloaded} OK, {errors} erros")
print(f"Brazil: {len(list(BRAZIL.glob('*.fa')))} FASTAs")
print(f"Spain:  {len(list(SPAIN.glob('*.fa')))} FASTAs")

# ── Final verification ──────────────────────────────────────────────────────
print("\n=== VERIFICAÇÃO FINAL ===")
brazil_local = {f.stem for f in BRAZIL.glob("*.fa")}
spain_local  = {f.stem for f in SPAIN.glob("*.fa")}
all_local    = brazil_local | spain_local

missing_local = [m for m in rows if m not in all_local]
if missing_local:
    print(f"Ainda faltando ({len(missing_local)}):")
    for m in sorted(missing_local):
        print(f"  {m}")
else:
    print(f"TODOS OS 729 MAGs BAIXADOS E VERIFICADOS.")
    print(f"  Brazil: {len(brazil_local)}")
    print(f"  Spain:  {len(spain_local)}")

# ── Taxonomy completeness ───────────────────────────────────────────────────
print("\n=== COBERTURA DE TAXONOMIA ===")
by_status = {"resolved": [], "needs_ena_registration": [], "unresolved": [], "other": []}
for mag_id, row in rows.items():
    status = row.get("tax_id_status", "other").strip() or "other"
    by_status.setdefault(status, []).append(mag_id)

for status, mags in by_status.items():
    print(f"  {status}: {len(mags)}")

unresolved = by_status.get("unresolved", [])
if unresolved:
    print(f"\nMAGs com taxonomia UNRESOLVED ({len(unresolved)}) — precisam de registro ENA ou revisão:")
    for mag_id in sorted(unresolved)[:20]:
        r = rows[mag_id]
        print(f"  {mag_id}: gtdb='{r.get('gtdb_species','')}' taxon='{r.get('taxon_name_for_submission','')}'")
    if len(unresolved) > 20:
        print(f"  ... e mais {len(unresolved)-20}")

needs_reg = by_status.get("needs_ena_registration", [])
if needs_reg:
    print(f"\nMAGs que precisam de registro de taxon no ENA ({len(needs_reg)}):")
    for mag_id in sorted(needs_reg):
        r = rows[mag_id]
        print(f"  {mag_id}: taxon='{r.get('taxon_name_for_submission','')}'")

print(f"\nTotal com taxonomia pronta (resolved + needs_ena_registration): {len(by_status.get('resolved',[])) + len(needs_reg)}")
print(f"Total sem taxonomia (unresolved): {len(unresolved)}")
