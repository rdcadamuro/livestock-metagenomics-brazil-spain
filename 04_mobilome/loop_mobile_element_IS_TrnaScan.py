#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MOBILOME PIPELINE v2.2 (BIOLOGICALLY ACCURATE CATEGORIES)
---------------------------------------------------------
Changes:
- Interprets 'gene' with 'tnpA/transposase' as 'Mobilization_Machinery' -> 'Transposase'.
- Interprets 'tnpR' as 'Resolvase'.
- Separates Structural parts (IRs, attC) from Enzmyes and Elements.
"""

import os
import sys
import time
import csv
import json
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from shutil import which as shutil_which
import psutil

# =========================
# CONFIGURAÇÕES
# =========================
# USER: set your base path here
IN_BASE_DIR = Path("/mnt/nvme2/E_coli_colistin_Spain/Fasta")
# USER: set your base path here
OUT_BASE = Path("/mnt/nvme2/E_coli_colistin_Spain/Results/Mobilome")

# Pastas de Saída
OUT_TRNA_BASE = OUT_BASE / "tRNA"
OUT_MEF_BASE  = OUT_BASE / "MEF"
OUT_INT_BASE  = OUT_BASE / "IntegronFinder"
OUT_MOB_BASE  = OUT_BASE / "MOBsuite"
OUT_TMP_BASE  = OUT_BASE / "_TMP"

# ARQUIVO FINAL
# USER: set your base path here
OUT_FINAL_TSV = Path("/mnt/nvme2/E_coli_colistin_Spain/Results/Mobilome_Detailed_v2.2.tsv")

MAX_SAMPLES_PARALLEL = 3
RAM_LIMIT_PERCENT = 80

# Timeouts & Threads
TIMEOUT_TRNA, TIMEOUT_MEF, TIMEOUT_INT, TIMEOUT_MOB = 900, 3600, 3600, 7200
MEF_THREADS, INT_THREADS, MOB_THREADS = "16", "12", "12"

DOCKER_BIN = "docker"
_trna = Path(os.environ.get("CONDA_PREFIX", "")) / "bin" / "tRNAscan-SE"
TRNA_BIN = _trna if _trna.exists() else Path(shutil_which("tRNAscan-SE") or "")
# USER: set your base path here
INTEGRON_FINDER_BIN = Path("/mnt/nvme/conda_envs/bioinfo/bin/integron_finder")

MEF_IMAGE = "mkhj/mobile_element_finder:latest"
MOB_IMAGE = "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1"

ENV = os.environ.copy()
ENV.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
DOCKER_UID = str(os.getuid()) if hasattr(os, "getuid") else None
DOCKER_GID = str(os.getgid()) if hasattr(os, "getgid") else None

# =========================
# HELPERS
# =========================
def wait_ram(limit=RAM_LIMIT_PERCENT):
    while psutil.virtual_memory().percent > limit:
        time.sleep(10)

def run_cmd(cmd, log_file: Path, timeout=None, cwd=None, extra_env=None):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    env = ENV.copy()
    if extra_env: env.update(extra_env)
    try:
        with open(log_file, "w") as lf:
            subprocess.run(cmd, stdout=lf, stderr=lf, env=env, check=True, timeout=timeout, cwd=cwd)
        return True
    except: return False

def looks_like_fasta(path: Path) -> bool:
    try: return path.stat().st_size >= 100 and open(path).read(100).strip().startswith(">")
    except: return False

def get_sample_id(fa_path: Path) -> str:
    name = fa_path.name
    for suf in [".fasta", ".fa", ".fna"]:
        if name.endswith(suf): return name[:-len(suf)]
    return fa_path.stem

def discover_fastas() -> list[Path]:
    files = []
    for pat in ["*.fna", "*.fa", "*.fasta"]:
        files.extend([p for p in IN_BASE_DIR.glob(pat) if looks_like_fasta(p)])
    return sorted(list(set(files)), key=lambda p: p.stat().st_size)

def is_done(sid: str) -> dict:
    return {
        "trna": (OUT_TRNA_BASE / f"{sid}.trnascan.txt").exists(),
        "mef": (OUT_MEF_BASE / sid / f"{sid}.gff").exists(),
        "int": (OUT_INT_BASE / sid).exists() and any((OUT_INT_BASE / sid).rglob("*.summary")),
        "mob": (OUT_MOB_BASE / sid).exists() and any((OUT_MOB_BASE / sid).rglob("contig_report*"))
    }

# =========================
# RUNNERS (Omitidos detalhes repetitivos, mantendo lógica)
# =========================
def run_trnascan(fa, sid):
    run_cmd([str(TRNA_BIN), "-B", "-o", str(OUT_TRNA_BASE/f"{sid}.trnascan.txt"), "-f", str(OUT_TRNA_BASE/f"{sid}.ss"), str(fa)], OUT_TRNA_BASE/f"{sid}.log", TIMEOUT_TRNA)
    return "OK"

def run_mefinder(fa, sid):
    out = OUT_MEF_BASE/sid
    out.mkdir(parents=True, exist_ok=True)
    cmd = [DOCKER_BIN, "run", "--rm", "-v", f"{fa.parent}:/in:ro", "-v", f"{out.resolve()}:/out"]
    if DOCKER_UID: cmd += ["--user", f"{DOCKER_UID}:{DOCKER_GID}"]
    cmd += [MEF_IMAGE, "mefinder", "find", "-c", f"/in/{fa.name}", "--threads", MEF_THREADS, "-g", f"/out/{sid}", "--gff", "--json"]
    run_cmd(cmd, out/f"{sid}.log", TIMEOUT_MEF)
    return "OK"

def run_integron(fa, sid):
    out = OUT_INT_BASE/sid
    out.mkdir(parents=True, exist_ok=True)
    run_cmd([str(INTEGRON_FINDER_BIN), str(fa), "--outdir", str(out), "--cpu", INT_THREADS], out/f"{sid}.log", TIMEOUT_INT)
    return "OK"

def run_mobsuite(fa, sid):
    out = OUT_MOB_BASE/sid
    out.mkdir(parents=True, exist_ok=True)
    tmp = OUT_TMP_BASE/sid
    tmp.mkdir(parents=True, exist_ok=True)
    cmd = [DOCKER_BIN, "run", "--rm", "-v", f"{fa.parent}:/in:ro", "-v", f"{out.resolve()}:/out", "-v", f"{tmp.resolve()}:/tmp", "-e", "TMPDIR=/tmp", "-w", "/out", MOB_IMAGE, "mob_recon", "-i", f"/in/{fa.name}", "-o", f"/out/{sid}", "--force", "--num_threads", MOB_THREADS]
    run_cmd(cmd, out/f"{sid}.log", TIMEOUT_MOB)
    return "OK"

def process_sample(fa_path: Path):
    sid = get_sample_id(fa_path)
    s = is_done(sid)
    wait_ram()
    if not s["trna"]: run_trnascan(fa_path, sid)
    if not s["mef"]: wait_ram(); run_mefinder(fa_path, sid)
    if not s["int"]: wait_ram(); run_integron(fa_path, sid)
    if not s["mob"]: wait_ram(); run_mobsuite(fa_path, sid)
    return sid

# =========================================================
# PARSERS "INTELIGENTES" (Aqui está a lógica das categorias)
# =========================================================

def parse_gff_attributes(attr_string):
    attrs = {}
    for item in attr_string.split(';'):
        if '=' in item:
            key, val = item.split('=', 1)
            attrs[key.strip().lower()] = val.strip() # Chaves em lowercase para facilitar
    return attrs

def parse_mef_rich(sid):
    results = []
    gff_path = OUT_MEF_BASE / sid / f"{sid}.gff"
    if not gff_path.exists(): return results

    try:
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                
                contig = parts[0]
                feature_type = parts[2].lower() # insertion_sequence, mite, gene, irl...
                start, end = int(parts[3]), int(parts[4])
                strand = parts[6]
                
                # Parse Atributos
                attr_dict = parse_gff_attributes(parts[8])
                name_val = attr_dict.get("name", "").lower()
                id_val = attr_dict.get("id", "")
                
                # --- LÓGICA DE CATEGORIZAÇÃO AVANÇADA ---
                category = "Mobile_Element"
                spec_type = "Unknown"
                
                # 1. É o Elemento Inteiro?
                if feature_type in ["insertion_sequence", "mite", "mge_domain", "transposon", "composite_transposon"]:
                    category = "Mobile_Element"
                    # MEF costuma por o tipo em 'mobile_element_type' nos atributos, ou usa o feature_type
                    spec_type = attr_dict.get("mobile_element_type", feature_type).replace("mge_", "")
                    details = f"ID:{id_val};Family:{attr_dict.get('name', 'NA')}"

                # 2. É uma Parte Estrutural? (IRL, IRR)
                elif feature_type in ["irl", "irr", "repeat_region"]:
                    category = "ME_Structure"
                    spec_type = "Inverted_Repeat"
                    details = f"Parent:{attr_dict.get('parent', 'NA')}"
                    
                # 3. É um Gene? (Aqui resolvemos Transposase e Resolvase)
                elif feature_type == "gene":
                    # Checamos o nome do gene
                    if "transposase" in name_val or "tnpa" in name_val:
                        category = "Mobilization_Machinery"
                        spec_type = "Transposase"
                    elif "resolvase" in name_val or "tnpr" in name_val:
                        category = "Mobilization_Machinery"
                        spec_type = "Resolvase"
                    elif "integrase" in name_val:
                        category = "Mobilization_Machinery"
                        spec_type = "Integrase"
                    else:
                        category = "Passenger_Gene"
                        spec_type = "Other_Protein"
                    
                    details = f"GeneName:{attr_dict.get('name', 'NA')};Parent:{attr_dict.get('parent', 'NA')}"

                else:
                    # Qualquer outra coisa
                    category = "Other_Feature"
                    spec_type = feature_type
                    details = parts[8]

                results.append({
                    "Sample_ID": sid, "Tool": "MobileElementFinder", "Contig_ID": contig,
                    "Start": start, "End": end, "Strand": strand,
                    "Category": category, 
                    "Specific_Type": spec_type,
                    "Details": details
                })
    except Exception as e: pass
    return results

def parse_trna_rich(sid):
    results = []
    fpath = OUT_TRNA_BASE / f"{sid}.trnascan.txt"
    if not fpath.exists(): return results
    try:
        with open(fpath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 8 or parts[0] in ["Sequence", "Name"] or parts[0].startswith("---"): continue
                try:
                    start, end = int(parts[2]), int(parts[3])
                    strand = "+" if start < end else "-"
                    if start > end: start, end = end, start
                    results.append({
                        "Sample_ID": sid, "Tool": "tRNAscan-SE", "Contig_ID": parts[0],
                        "Start": start, "End": end, "Strand": strand,
                        "Category": "tRNA",
                        "Specific_Type": parts[4], # Amino acid
                        "Details": f"Anticodon:{parts[5]};Score:{parts[8]}"
                    })
                except: continue
    except: pass
    return results

def parse_int_rich(sid):
    results = []
    out_dir = OUT_INT_BASE / sid
    if not out_dir.exists(): return results
    int_files = list(out_dir.rglob("*.integrons"))
    if not int_files: return results
    
    try:
        with open(int_files[0], 'r') as f:
            lines = [l for l in f if not l.startswith("#")]
            reader = csv.DictReader(lines, delimiter='\t')
            for row in reader:
                try:
                    annot = row.get("annotation", "NA")
                    elt_type = row.get("type_elt", "unknown")
                    
                    cat = "Integron_Feature"
                    spec = elt_type
                    
                    # MELHORA A CATEGORIA DO INTEGRON
                    if "intI" in annot or "integrase" in annot:
                        cat = "Mobilization_Machinery"
                        spec = "Integrase"
                    elif elt_type == "attC":
                        cat = "ME_Structure"
                        spec = "attC_Site"
                    elif elt_type == "protein":
                        cat = "Passenger_Gene"
                        spec = "Gene_Cassette"
                        
                    s_val = row.get("strand", "0")
                    strand = "+" if s_val == "1" else "-" if s_val == "-1" else "."
                    
                    results.append({
                        "Sample_ID": sid, "Tool": "IntegronFinder", "Contig_ID": row.get("ID_replicon"),
                        "Start": int(row.get("pos_beg")), "End": int(row.get("pos_end")), "Strand": strand,
                        "Category": cat, "Specific_Type": spec,
                        "Details": f"Annot:{annot}"
                    })
                except: continue
    except: pass
    return results

def parse_mob_rich(sid):
    results = []
    report = OUT_MOB_BASE / sid / "contig_report.txt"
    if not report.exists(): return results
    try:
        with open(report, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                rep = row.get("rep_type(s)", "-")
                if not rep or rep == "-": rep = "Unclassified"
                results.append({
                    "Sample_ID": sid, "Tool": "MOB-suite", "Contig_ID": row.get("contig_id"),
                    "Start": 1, "End": int(row.get("size", 0)), "Strand": ".",
                    "Category": "Plasmid", "Specific_Type": rep,
                    "Details": f"Relaxase:{row.get('relaxase_type(s)')};Cluster:{row.get('primary_cluster_id')}"
                })
    except: pass
    return results

def compile_full_tsv(fasta_files):
    print(f"\n📑 Compilando TSV V2.2 (Categorias Corrigidas) em: {OUT_FINAL_TSV}")
    fieldnames = ["Sample_ID", "Tool", "Contig_ID", "Start", "End", "Strand", "Category", "Specific_Type", "Details"]
    
    try:
        OUT_FINAL_TSV.parent.mkdir(parents=True, exist_ok=True)
        with open(OUT_FINAL_TSV, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            count = 0
            for fa in fasta_files:
                sid = get_sample_id(fa)
                rows = []
                rows.extend(parse_trna_rich(sid))
                rows.extend(parse_mef_rich(sid))
                rows.extend(parse_int_rich(sid))
                rows.extend(parse_mob_rich(sid))
                if rows:
                    rows.sort(key=lambda x: (x["Contig_ID"], x["Start"]))
                    writer.writerows(rows)
                    count += len(rows)
            print(f"✅ Sucesso! Total de linhas: {count}")
    except Exception as e: print(f"❌ Erro: {e}")

# =========================
# MAIN
# =========================
def main():
    print("="*80)
    print("MOBILOME PIPELINE v2.2 (SMART CATEGORIES)")
    print("="*80)
    
    if not TRNA_BIN.exists() or not INTEGRON_FINDER_BIN.exists():
        sys.exit("Erro: Binários não encontrados.")

    fastas = discover_fastas()
    print(f"Processando {len(fastas)} amostras...")
    
    with ProcessPoolExecutor(max_workers=MAX_SAMPLES_PARALLEL) as ex:
        list(ex.map(process_sample, fastas)) # Força execução
            
    compile_full_tsv(fastas)

if __name__ == "__main__":
    main()