#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline MAGs: tRNAscan-SE (local) + MobileElementFinder (Docker)
-----------------------------------------------------------------
Entrada (MAGs):
  /mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/Fasta/<SAMPLE>/.../bin.N.(fa|fna|fasta)[.gz]

Saídas:
- tRNAscan:
  OUT_TRNA_BASE/<SAMPLE>__<MAG>.trnascan.txt
  OUT_TRNA_BASE/<SAMPLE>__<MAG>.trnascan.ss
  OUT_TRNA_BASE/<SAMPLE>__<MAG>.trnascan.log

- MobileElementFinder (Docker):
  OUT_MEF_BASE/<SAMPLE>/<MAG>/*
    (gff/json/csv + log)

- TSV unificado:
  OUT_SUMMARY_TSV
  colunas: Sample, MAG, contig, start, end, type, program

Patches:
- Descoberta para MAGs (bin.*.fa/fna/fasta)
- Sample/MAG derivados do path
- TRNA_BIN valida executável (nunca usa '.')
- Suporte a .gz via TMP
- Parser MEF procura *.gff no diretório (não assume nome fixo)
"""

import os
import re
import sys
import time
import gzip
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from shutil import which as shutil_which
from typing import Optional, Tuple, List

import psutil

# =========================
# CONFIG
# =========================

# USER: set your base path here

IN_BASE_DIR = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/Fasta")

# USER: set your base path here

OUT_TRNA_BASE = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/tRNA")
# USER: set your base path here
OUT_MEF_BASE  = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/MEF")
# USER: set your base path here
OUT_SUMMARY_TSV = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/Mobilome_All_Results.tsv")

# USER: set your base path here

TMP_BASE = Path("/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Mags/_TMP_MOBILOME")

# Docker
DOCKER_BIN = "docker"
MEF_IMAGE = "mkhj/mobile_element_finder:latest"
MEF_THREADS = "16"

# Paralelismo por MAG (job = 1 MAG)
MAX_JOBS_PARALLEL = 3

# RAM
RAM_LIMIT_PERCENT = 80

# Timeouts (s)
TIMEOUT_TRNA = 900
TIMEOUT_MEF  = 3600

# Ambiente (limitar threads implícitas)
ENV = os.environ.copy()
ENV.update({
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
})

DOCKER_UID = str(os.getuid()) if hasattr(os, "getuid") else None
DOCKER_GID = str(os.getgid()) if hasattr(os, "getgid") else None

# =========================
# HELPERS
# =========================

SAMPLE_RE = re.compile(r"^P\d+_\d+$", re.I)
MAG_RE = re.compile(r"(bin\.\d+)", re.I)

def find_trnascan_bin() -> Optional[Path]:
    """
    Retorna Path executável do tRNAscan-SE, ou None.
    Nunca retorna '.'.
    """
    candidates: List[Path] = []
    conda = os.environ.get("CONDA_PREFIX", "")
    if conda:
        candidates.append(Path(conda) / "bin" / "tRNAscan-SE")

    w = shutil_which("tRNAscan-SE")
    if w:
        candidates.append(Path(w))

    for cand in candidates:
        if cand.exists() and cand.is_file() and os.access(str(cand), os.X_OK):
            return cand

    return None

def which(cmd: str) -> Optional[str]:
    return shutil_which(cmd)

def get_ram_percent() -> float:
    return psutil.virtual_memory().percent

def wait_ram(limit=RAM_LIMIT_PERCENT):
    while get_ram_percent() > limit:
        print(f"    ⚠️  RAM: {get_ram_percent():.1f}% (aguardando...)")
        time.sleep(10)

def run_cmd(cmd, log_file: Path, timeout=None, cwd=None) -> bool:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(log_file, "w") as lf:
            subprocess.run(
                cmd, stdout=lf, stderr=lf,
                env=ENV, check=True, timeout=timeout, cwd=cwd
            )
        return True
    except subprocess.TimeoutExpired:
        print(f"      ⏱ TIMEOUT: {log_file.name}")
        return False
    except subprocess.CalledProcessError as e:
        print(f"      ✗ ERRO (exit {e.returncode}): {log_file.name}")
        return False
    except Exception as e:
        print(f"      ✗ ERRO: {e}")
        return False

def derive_sample_mag(fa_path: Path) -> Tuple[str, str]:
    """
    Deriva Sample e MAG do path do MAG fasta.
    Exemplo:
      .../Mags/Fasta/P19109_1035/metawrap_50_10_bins/bin.1.fa
      -> Sample=P19109_1035, MAG=bin.1
    """
    sample = ""
    for parent in fa_path.parents:
        if SAMPLE_RE.match(parent.name):
            sample = parent.name
            break
    if not sample:
        # fallback: tenta achar no path
        m = re.search(r"(P\d+_\d+)", str(fa_path), flags=re.I)
        sample = m.group(1) if m else "UNKNOWN_SAMPLE"

    m = MAG_RE.search(fa_path.name)
    if not m:
        # tenta no path inteiro
        m = MAG_RE.search(str(fa_path))
    mag = m.group(1).lower() if m else fa_path.stem.lower()

    return sample, mag

def ensure_plain_fasta(fa_path: Path, sample: str, mag: str) -> Path:
    """
    Se for .gz, descomprime para TMP_BASE/<sample>/<mag>/<original_sem_gz>.
    Retorna path do FASTA plain.
    """
    if not fa_path.name.endswith(".gz"):
        return fa_path

    out_dir = TMP_BASE / sample / mag
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / fa_path.name[:-3]  # remove .gz
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    with gzip.open(fa_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_path

def is_done(sample: str, mag: str) -> dict:
    key = f"{sample}__{mag}"
    trna_txt = OUT_TRNA_BASE / f"{key}.trnascan.txt"

    mef_dir = OUT_MEF_BASE / sample / mag
    # MEF outputs variam, então basta existir algo relevante
    gffs = list(mef_dir.glob("*.gff")) if mef_dir.exists() else []
    jsons = list(mef_dir.glob("*.json")) if mef_dir.exists() else []
    csvs = list(mef_dir.glob("*.csv")) if mef_dir.exists() else []

    return {
        "trna": trna_txt.exists() and trna_txt.stat().st_size > 0,
        "mef": any(p.stat().st_size > 0 for p in (gffs + jsons + csvs)) if mef_dir.exists() else False,
    }

# =========================
# PARSERS
# =========================

def parse_trnascan(trna_txt: Path, sample: str, mag: str):
    """
    Retorna:
      (Sample, MAG, contig, start, end, type, program)
    """
    rows = []
    if not trna_txt.exists() or trna_txt.stat().st_size == 0:
        return rows

    with open(trna_txt, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Sequence") or line.startswith("Name") or line.startswith("-") or line.startswith("#"):
                continue

            parts = re.split(r"\s+", line)
            if len(parts) < 5:
                continue

            contig = parts[0]
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue

            s = min(start, end)
            e = max(start, end)
            ttype = parts[4]
            rows.append((sample, mag, contig, s, e, ttype, "tRNAscan-SE"))

    return rows

def parse_any_gff_in_dir(mef_dir: Path, sample: str, mag: str):
    """
    Lê todos *.gff no diretório do MEF.
    Retorna:
      (Sample, MAG, contig, start, end, type, program)
    """
    rows = []
    if not mef_dir.exists():
        return rows

    gffs = sorted(mef_dir.glob("*.gff"))
    for gff_path in gffs:
        if gff_path.stat().st_size == 0:
            continue
        with open(gff_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                contig = parts[0]
                ftype = parts[2] if len(parts) >= 3 else "MGE"
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    continue
                s = min(start, end)
                e = max(start, end)
                rows.append((sample, mag, contig, s, e, ftype, "MobileElementFinder"))
    return rows

# =========================
# RUNNERS
# =========================

def run_trnascan(TRNA_BIN: Path, fasta_plain: Path, sample: str, mag: str) -> str:
    key = f"{sample}__{mag}"
    trna_txt = OUT_TRNA_BASE / f"{key}.trnascan.txt"
    trna_ss  = OUT_TRNA_BASE / f"{key}.trnascan.ss"
    trna_log = OUT_TRNA_BASE / f"{key}.trnascan.log"

    cmd = [str(TRNA_BIN), "-B", "-o", str(trna_txt), "-f", str(trna_ss), str(fasta_plain)]
    print(f"    → tRNAscan: {sample} | {mag}")
    ok = run_cmd(cmd, trna_log, timeout=TIMEOUT_TRNA)
    return "OK" if ok else "FAIL"

def run_mefinder_docker(fasta_plain: Path, sample: str, mag: str) -> str:
    out_dir = OUT_MEF_BASE / sample / mag
    out_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / f"{sample}__{mag}.mefinder.log"

    in_dir = fasta_plain.parent.resolve()
    in_file = fasta_plain.name

    # prefixo dentro do container (arquivos serão /out/<prefix>.*)
    out_prefix = f"/out/{sample}__{mag}"

    cmd = [
        DOCKER_BIN, "run", "--rm",
        "-v", f"{str(in_dir)}:/in:ro",
        "-v", f"{str(out_dir.resolve())}:/out",
    ]
    if DOCKER_UID is not None and DOCKER_GID is not None:
        cmd += ["--user", f"{DOCKER_UID}:{DOCKER_GID}"]

    cmd += [
        MEF_IMAGE,
        "mefinder", "find",
        "-c", f"/in/{in_file}",
        "--threads", str(MEF_THREADS),
        "-g", out_prefix,
        "--gff", "--json"
    ]

    print(f"    → MEFinder: {sample} | {mag}")
    ok = run_cmd(cmd, log, timeout=TIMEOUT_MEF)
    return "OK" if ok else "FAIL"

def process_one_mag(TRNA_BIN: Path, fa_path: Path):
    sample, mag = derive_sample_mag(fa_path)
    status = is_done(sample, mag)

    wait_ram()
    results = {"sample": sample, "mag": mag, "trna": "SKIP", "mef": "SKIP"}

    fasta_plain = ensure_plain_fasta(fa_path, sample, mag)

    if not status["trna"]:
        results["trna"] = run_trnascan(TRNA_BIN, fasta_plain, sample, mag)

    if not status["mef"]:
        wait_ram(max(10, RAM_LIMIT_PERCENT - 10))
        results["mef"] = run_mefinder_docker(fasta_plain, sample, mag)

    return results

# =========================
# DISCOVERY (MAGs)
# =========================

def discover_mag_fastas(base: Path):
    """
    Descobre FASTAs de MAG:
      bin.*.(fa|fna|fasta) e versões .gz
    Evita lixo como 'viruses', '1', etc.
    """
    patterns = [
        "bin.*.fa", "bin.*.fna", "bin.*.fasta",
        "bin.*.fa.gz", "bin.*.fna.gz", "bin.*.fasta.gz",
    ]
    files = []
    for pat in patterns:
        files.extend(base.rglob(pat))

    out = []
    seen = set()
    for p in files:
        if not p.is_file():
            continue
        try:
            if p.stat().st_size <= 0:
                continue
        except OSError:
            continue
        rp = str(p.resolve())
        if rp in seen:
            continue
        seen.add(rp)
        out.append(p)
    return sorted(out)

# =========================
# SUMMARY TSV
# =========================

def write_unified_tsv():
    OUT_SUMMARY_TSV.parent.mkdir(parents=True, exist_ok=True)
    rows_all = []

    # tRNAscan
    for trna_txt in sorted(OUT_TRNA_BASE.glob("*.trnascan.txt")):
        key = trna_txt.name.replace(".trnascan.txt", "")
        if "__" in key:
            sample, mag = key.split("__", 1)
        else:
            sample, mag = "UNKNOWN_SAMPLE", key
        rows_all.extend(parse_trnascan(trna_txt, sample, mag))

    # MEFinder
    if OUT_MEF_BASE.exists():
        for sample_dir in sorted(OUT_MEF_BASE.iterdir()):
            if not sample_dir.is_dir():
                continue
            sample = sample_dir.name
            for mag_dir in sorted(sample_dir.iterdir()):
                if not mag_dir.is_dir():
                    continue
                mag = mag_dir.name
                rows_all.extend(parse_any_gff_in_dir(mag_dir, sample, mag))

    with open(OUT_SUMMARY_TSV, "w", encoding="utf-8") as out:
        out.write("Sample\tMAG\tcontig\tstart\tend\ttype\tprogram\n")
        for r in rows_all:
            out.write("\t".join(map(str, r)) + "\n")

    return len(rows_all)

# =========================
# MAIN
# =========================

def main():
    print("=" * 70)
    print("Pipeline MAGs: tRNAscan-SE (local) + MobileElementFinder (Docker)")
    print(f"Jobs simultâneos: {MAX_JOBS_PARALLEL}")
    print(f"Threads: MEF({MEF_THREADS})")
    print(f"RAM total: {psutil.virtual_memory().total / (1024**3):.1f} GB")
    print("=" * 70)

    OUT_TRNA_BASE.mkdir(parents=True, exist_ok=True)
    OUT_MEF_BASE.mkdir(parents=True, exist_ok=True)
    TMP_BASE.mkdir(parents=True, exist_ok=True)

    TRNA_BIN = find_trnascan_bin()
    if TRNA_BIN is None:
        print("❌ ERRO: tRNAscan-SE não encontrado como executável.")
        print("   Faça:")
        print("     conda activate bioinfo")
        print("     which tRNAscan-SE")
        sys.exit(1)
    print(f"[INFO] tRNAscan-SE = {TRNA_BIN}")

    if which(DOCKER_BIN) is None:
        print(f"❌ ERRO: '{DOCKER_BIN}' não encontrado no PATH.")
        sys.exit(1)

    mag_fastas = discover_mag_fastas(IN_BASE_DIR)
    if not mag_fastas:
        print(f"❌ Nenhum MAG FASTA (bin.*.fa|fna|fasta[.gz]) encontrado em {IN_BASE_DIR}")
        sys.exit(1)

    # Ordena por tamanho
    mag_fastas.sort(key=lambda p: p.stat().st_size if p.exists() else 0)
    total = len(mag_fastas)

    print(f"📂 Encontrados {total} MAG FASTAs")
    print(f"▶ Iniciando (RAM: {get_ram_percent():.1f}%)\n")

    stats = {
        "trna_ok": 0, "trna_skip": 0, "trna_fail": 0,
        "mef_ok": 0,  "mef_skip": 0,  "mef_fail": 0,
    }

    completed = 0
    with ProcessPoolExecutor(max_workers=MAX_JOBS_PARALLEL) as executor:
        futures = {executor.submit(process_one_mag, TRNA_BIN, fp): fp for fp in mag_fastas}
        for future in as_completed(futures):
            completed += 1
            fp = futures[future]
            try:
                res = future.result()
                sample = res["sample"]
                mag = res["mag"]

                for tool in ["trna", "mef"]:
                    st = res[tool]
                    if st == "OK":
                        stats[f"{tool}_ok"] += 1
                    elif st == "SKIP":
                        stats[f"{tool}_skip"] += 1
                    else:
                        stats[f"{tool}_fail"] += 1

                print(f"  [{completed}/{total}] {sample} | {mag} → T:{res['trna']} MEF:{res['mef']} | RAM:{get_ram_percent():.1f}%")

            except Exception as e:
                sample, mag = derive_sample_mag(fp)
                print(f"  [{completed}/{total}] ✗ ERRO: {sample} | {mag} → {e}")

    print("\n" + "=" * 70)
    print("PIPELINE CONCLUÍDO")
    print("=" * 70)
    print(f"RAM final: {get_ram_percent():.1f}%\n")

    print("📊 ESTATÍSTICAS:")
    print(f"  tRNAscan → OK: {stats['trna_ok']} | SKIP: {stats['trna_skip']} | FAIL: {stats['trna_fail']}")
    print(f"  MEFinder → OK: {stats['mef_ok']} | SKIP: {stats['mef_skip']} | FAIL: {stats['mef_fail']}")

    total_fails = stats["trna_fail"] + stats["mef_fail"]
    if total_fails > 0:
        print(f"\n⚠️  {total_fails} tarefas falharam. Logs em:")
        print(f"    {OUT_TRNA_BASE}/*.trnascan.log")
        print(f"    {OUT_MEF_BASE}/*/*/*.mefinder.log")
    else:
        print("\n✅ Todas as tarefas concluídas com sucesso!")

    nrows = write_unified_tsv()
    print("\n📦 TSV UNIFICADO GERADO:")
    print(f"    {OUT_SUMMARY_TSV}")
    print(f"    Linhas: {nrows}")

    print("\n📁 RESULTADOS:")
    print(f"    tRNAscan → {OUT_TRNA_BASE}")
    print(f"    MEFinder → {OUT_MEF_BASE}/<Sample>/<MAG>/")
    print(f"    TMP     → {TMP_BASE} (somente .gz descomprimidos)")
    print("=" * 70)

if __name__ == "__main__":
    main()
