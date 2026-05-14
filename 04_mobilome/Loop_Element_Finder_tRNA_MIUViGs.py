#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline MIUViGs: tRNAscan-SE (local) + MobileElementFinder (Docker)
--------------------------------------------------------------------
Entrada (MIUViGs por amostra):
  /mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/vRhyme_test_P19109_1035/Anotation/<SAMPLE>/_sample_level/P19109_<SAMPLE>.vMAGs.clean.merged.fasta

Saídas (Results):
  /mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/MIUViG_ElementFinder_tRNA/
    - tRNA/<SAMPLE>__<UNIT>.trnascan.*
    - MEF/<SAMPLE>/<UNIT>/*
    - MIUViG_Mobilome_All_Results.tsv

Notas:
- Usa 20 threads no MobileElementFinder.
- Paralelismo por unidade definido para 1 (20 threads totais efetivos).
- Inclui coluna Bin derivada do contig: vRhyme_<n>__* -> vRhyme_bin_<n>.
"""

from __future__ import annotations

import os
import re
import sys
import time
import gzip
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from shutil import which as shutil_which
from typing import Optional, Tuple, List

import psutil

# =========================
# CONFIG
# =========================

IN_BASE_DIR = Path(
    os.environ.get(
        "MIUVIG_MOBILOME_IN",
        "/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/vRhyme_test_P19109_1035/Anotation",
    )
)

OUT_BASE = Path(
    os.environ.get(
        "MIUVIG_MOBILOME_OUT",
        "/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly/Results/MIUViG_ElementFinder_tRNA",
    )
)
OUT_TRNA_BASE = OUT_BASE / "tRNA"
OUT_MEF_BASE = OUT_BASE / "MEF"
OUT_SUMMARY_TSV = OUT_BASE / "MIUViG_Mobilome_All_Results.tsv"
TMP_BASE = OUT_BASE / "_TMP_MOBILOME"

# Docker
DOCKER_BIN = "docker"
MEF_IMAGE = "mkhj/mobile_element_finder:latest"
MEF_THREADS = "20"

# Paralelismo (1 job x 20 threads)
MAX_JOBS_PARALLEL = 1

# RAM
RAM_LIMIT_PERCENT = 80

# Timeouts (s)
TIMEOUT_TRNA = 1200
TIMEOUT_MEF = 7200

# Ambiente (limitar threads implícitas fora do MEF)
ENV = os.environ.copy()
ENV.update(
    {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
)

DOCKER_UID = str(os.getuid()) if hasattr(os, "getuid") else None
DOCKER_GID = str(os.getgid()) if hasattr(os, "getgid") else None

SAMPLE_RE = re.compile(r"^P\d+_(\d+)$", re.I)
SAMPLE_NUM_RE = re.compile(r"^\d{4,}$")
SAMPLE_FROM_NAME_RE = re.compile(r"P\d+_(\d+)", re.I)
BIN_RE = re.compile(r"^vRhyme_(\d+)__")


# =========================
# HELPERS
# =========================
def find_trnascan_bin() -> Optional[Path]:
    candidates: List[Path] = []
    conda = os.environ.get("CONDA_PREFIX", "")
    if conda:
        candidates.append(Path(conda) / "bin" / "tRNAscan-SE")
    # Common local envs on this workstation
    candidates.append(Path("/mnt/nvme/conda_envs/bioinfo/bin/tRNAscan-SE"))
    candidates.append(Path("/mnt/nvme/conda_envs/pharokka_env/bin/tRNAscan-SE"))
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


def wait_ram(limit: int = RAM_LIMIT_PERCENT):
    while get_ram_percent() > limit:
        print(f"    ⚠️ RAM: {get_ram_percent():.1f}% (aguardando...)")
        time.sleep(10)


def run_cmd(cmd, log_file: Path, timeout=None, cwd=None) -> bool:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(log_file, "w") as lf:
            subprocess.run(
                cmd,
                stdout=lf,
                stderr=lf,
                env=ENV,
                check=True,
                timeout=timeout,
                cwd=cwd,
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


def infer_bin(contig: str) -> str:
    m = BIN_RE.match(contig or "")
    if not m:
        return "NA"
    return f"vRhyme_bin_{m.group(1)}"


def derive_sample_unit(fa_path: Path) -> Tuple[str, str]:
    sample = "UNKNOWN_SAMPLE"
    for parent in fa_path.parents:
        m = SAMPLE_RE.match(parent.name)
        if m:
            sample = m.group(1)
            break
        if SAMPLE_NUM_RE.match(parent.name):
            sample = parent.name
            break
    if sample == "UNKNOWN_SAMPLE":
        m2 = SAMPLE_FROM_NAME_RE.search(fa_path.name)
        if m2:
            sample = m2.group(1)
    unit = fa_path.stem.lower().replace(".", "_")
    return sample, unit


def ensure_plain_fasta(fa_path: Path, sample: str, unit: str) -> Path:
    if not fa_path.name.endswith(".gz"):
        return fa_path
    out_dir = TMP_BASE / sample / unit
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / fa_path.name[:-3]
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    with gzip.open(fa_path, "rb") as fin, open(out_path, "wb") as fout:
        fout.write(fin.read())
    return out_path


def is_done(sample: str, unit: str) -> dict:
    key = f"{sample}__{unit}"
    trna_txt = OUT_TRNA_BASE / f"{key}.trnascan.txt"
    mef_dir = OUT_MEF_BASE / sample / unit
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
def parse_trnascan(trna_txt: Path, sample: str, unit: str):
    rows = []
    if not trna_txt.exists() or trna_txt.stat().st_size == 0:
        return rows
    with open(trna_txt, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(("Sequence", "Name", "-", "#")):
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
            rows.append((sample, unit, infer_bin(contig), contig, s, e, ttype, "tRNAscan-SE"))
    return rows


def parse_any_gff_in_dir(mef_dir: Path, sample: str, unit: str):
    rows = []
    if not mef_dir.exists():
        return rows
    for gff_path in sorted(mef_dir.glob("*.gff")):
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
                rows.append((sample, unit, infer_bin(contig), contig, s, e, ftype, "MobileElementFinder"))
    return rows


# =========================
# RUNNERS
# =========================
def run_trnascan(trna_bin: Path, fasta_plain: Path, sample: str, unit: str) -> str:
    key = f"{sample}__{unit}"
    trna_txt = OUT_TRNA_BASE / f"{key}.trnascan.txt"
    trna_ss = OUT_TRNA_BASE / f"{key}.trnascan.ss"
    trna_log = OUT_TRNA_BASE / f"{key}.trnascan.log"
    cmd = [str(trna_bin), "-B", "-o", str(trna_txt), "-f", str(trna_ss), str(fasta_plain)]
    print(f"    → tRNAscan: {sample} | {unit}")
    ok = run_cmd(cmd, trna_log, timeout=TIMEOUT_TRNA)
    return "OK" if ok else "FAIL"


def run_mefinder_docker(fasta_plain: Path, sample: str, unit: str) -> str:
    out_dir = OUT_MEF_BASE / sample / unit
    out_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / f"{sample}__{unit}.mefinder.log"
    in_dir = fasta_plain.parent.resolve()
    in_file = fasta_plain.name
    out_prefix = f"/out/{sample}__{unit}"
    cmd = [
        DOCKER_BIN,
        "run",
        "--rm",
        "-v",
        f"{str(in_dir)}:/in:ro",
        "-v",
        f"{str(out_dir.resolve())}:/out",
    ]
    if DOCKER_UID is not None and DOCKER_GID is not None:
        cmd += ["--user", f"{DOCKER_UID}:{DOCKER_GID}"]
    cmd += [
        MEF_IMAGE,
        "mefinder",
        "find",
        "-c",
        f"/in/{in_file}",
        "--threads",
        str(MEF_THREADS),
        "-g",
        out_prefix,
        "--gff",
        "--json",
    ]
    print(f"    → MEFinder: {sample} | {unit}")
    ok = run_cmd(cmd, log, timeout=TIMEOUT_MEF)
    return "OK" if ok else "FAIL"


def process_one_unit(trna_bin: Path, fa_path: Path):
    sample, unit = derive_sample_unit(fa_path)
    status = is_done(sample, unit)
    wait_ram()
    results = {"sample": sample, "unit": unit, "trna": "SKIP", "mef": "SKIP"}
    fasta_plain = ensure_plain_fasta(fa_path, sample, unit)
    if not status["trna"]:
        results["trna"] = run_trnascan(trna_bin, fasta_plain, sample, unit)
    if not status["mef"]:
        wait_ram(max(10, RAM_LIMIT_PERCENT - 10))
        results["mef"] = run_mefinder_docker(fasta_plain, sample, unit)
    return results


# =========================
# DISCOVERY
# =========================
def discover_miuvig_fastas(base: Path):
    pats = ["*/_sample_level/P19109_*.vMAGs.clean.merged.fasta", "*/_sample_level/P19109_*.vMAGs.clean.merged.fasta.gz"]
    files = []
    for pat in pats:
        files.extend(base.glob(pat))
    out = []
    seen = set()
    for p in files:
        if not p.is_file():
            continue
        if p.stat().st_size <= 0:
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
    for trna_txt in sorted(OUT_TRNA_BASE.glob("*.trnascan.txt")):
        key = trna_txt.name.replace(".trnascan.txt", "")
        if "__" in key:
            sample, unit = key.split("__", 1)
        else:
            sample, unit = "UNKNOWN_SAMPLE", key
        rows_all.extend(parse_trnascan(trna_txt, sample, unit))
    if OUT_MEF_BASE.exists():
        for sample_dir in sorted(OUT_MEF_BASE.iterdir()):
            if not sample_dir.is_dir():
                continue
            sample = sample_dir.name
            for unit_dir in sorted(sample_dir.iterdir()):
                if not unit_dir.is_dir():
                    continue
                unit = unit_dir.name
                rows_all.extend(parse_any_gff_in_dir(unit_dir, sample, unit))
    with open(OUT_SUMMARY_TSV, "w", encoding="utf-8") as out:
        out.write("Sample\tUnit\tBin\tcontig\tstart\tend\ttype\tprogram\n")
        for r in rows_all:
            out.write("\t".join(map(str, r)) + "\n")
    return len(rows_all)


# =========================
# MAIN
# =========================
def main():
    print("=" * 70)
    print("Pipeline MIUViGs: tRNAscan-SE (local) + MobileElementFinder (Docker)")
    print(f"Jobs simultâneos: {MAX_JOBS_PARALLEL}")
    print(f"Threads MEF: {MEF_THREADS}")
    print(f"RAM total: {psutil.virtual_memory().total / (1024**3):.1f} GB")
    print("=" * 70)

    OUT_TRNA_BASE.mkdir(parents=True, exist_ok=True)
    OUT_MEF_BASE.mkdir(parents=True, exist_ok=True)
    TMP_BASE.mkdir(parents=True, exist_ok=True)

    trna_bin = find_trnascan_bin()
    if trna_bin is None:
        print("❌ ERRO: tRNAscan-SE não encontrado.")
        sys.exit(1)
    print(f"[INFO] tRNAscan-SE = {trna_bin}")

    if which(DOCKER_BIN) is None:
        print(f"❌ ERRO: '{DOCKER_BIN}' não encontrado no PATH.")
        sys.exit(1)

    miuvig_fastas = discover_miuvig_fastas(IN_BASE_DIR)
    if not miuvig_fastas:
        print(f"❌ Nenhum MIUViG merged FASTA encontrado em {IN_BASE_DIR}")
        sys.exit(1)

    miuvig_fastas.sort(key=lambda p: p.stat().st_size if p.exists() else 0)
    total = len(miuvig_fastas)
    print(f"📂 Encontrados {total} FASTAs MIUViG")
    print(f"▶ Iniciando (RAM: {get_ram_percent():.1f}%)\n")

    stats = {"trna_ok": 0, "trna_skip": 0, "trna_fail": 0, "mef_ok": 0, "mef_skip": 0, "mef_fail": 0}
    completed = 0
    with ProcessPoolExecutor(max_workers=MAX_JOBS_PARALLEL) as executor:
        futures = {executor.submit(process_one_unit, trna_bin, fp): fp for fp in miuvig_fastas}
        for future in as_completed(futures):
            completed += 1
            fp = futures[future]
            try:
                res = future.result()
                sample = res["sample"]
                unit = res["unit"]
                for tool in ["trna", "mef"]:
                    st = res[tool]
                    if st == "OK":
                        stats[f"{tool}_ok"] += 1
                    elif st == "SKIP":
                        stats[f"{tool}_skip"] += 1
                    else:
                        stats[f"{tool}_fail"] += 1
                print(f"  [{completed}/{total}] {sample} | {unit} → T:{res['trna']} MEF:{res['mef']} | RAM:{get_ram_percent():.1f}%")
            except Exception as e:
                sample, unit = derive_sample_unit(fp)
                print(f"  [{completed}/{total}] ✗ ERRO: {sample} | {unit} → {e}")

    print("\n" + "=" * 70)
    print("PIPELINE CONCLUÍDO")
    print("=" * 70)
    print(f"RAM final: {get_ram_percent():.1f}%\n")
    print("📊 ESTATÍSTICAS:")
    print(f"  tRNAscan → OK: {stats['trna_ok']} | SKIP: {stats['trna_skip']} | FAIL: {stats['trna_fail']}")
    print(f"  MEFinder → OK: {stats['mef_ok']} | SKIP: {stats['mef_skip']} | FAIL: {stats['mef_fail']}")

    nrows = write_unified_tsv()
    print("\n📦 TSV UNIFICADO GERADO:")
    print(f"    {OUT_SUMMARY_TSV}")
    print(f"    Linhas: {nrows}")
    print("\n📁 RESULTADOS:")
    print(f"    tRNAscan → {OUT_TRNA_BASE}")
    print(f"    MEFinder → {OUT_MEF_BASE}/<Sample>/<Unit>/")
    print(f"    TMP      → {TMP_BASE}")
    print("=" * 70)


if __name__ == "__main__":
    main()
