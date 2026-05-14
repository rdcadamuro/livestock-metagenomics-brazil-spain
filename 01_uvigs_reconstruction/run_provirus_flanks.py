#!/usr/bin/env python3
"""
Extract 5 kb flanks from geNomad proviruses (coordinates != NA) and run
MobileElementFinder (docker) + IntegronFinder on them.
"""
import os, subprocess, textwrap
from pathlib import Path

# USER: set your base path here

GENOMAD_BASE = Path("/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Mags_Vmags/Genomad_Brazil")
# USER: set your base path here
OUT_BASE     = Path("/mnt/nvme2/RECOVERY_m2/Metagenome/Tati/Assembly/Results/Mags_Vmags/provirus_flanks")
FLANK_KB     = 5000
# USER: set your base path here
INTEGRON_BIN = Path("/mnt/nvme/conda_envs/bioinfo/bin/integron_finder")
MEF_IMAGE    = "mkhj/mobile_element_finder:latest"

ENV = os.environ.copy()
ENV["PATH"] = "/mnt/nvme/conda_envs/bioinfo/bin:" + ENV.get("PATH", "")

# ── collect proviruses ────────────────────────────────────────────────────────
def collect_proviruses():
    records = []
    for summary in sorted(GENOMAD_BASE.glob("P19109_*/downloaded_genomad/*/P19109_*_virus_summary.tsv")):
        sample = summary.parts[-3].split("_")[0] + "_" + summary.parts[-3].split("_")[1]  # P19109_XXXX
        with open(summary) as fh:
            hdr = fh.readline().strip().split("\t")
            ci  = hdr.index("coordinates")
            for line in fh:
                cols = line.strip().split("\t")
                coords = cols[ci]
                if coords == "NA" or not coords:
                    continue
                seq_name = cols[0]  # e.g. NODE_78_...|provirus_151_43221
                # parse start/end from coordinates field
                start, end = map(int, coords.split("-"))
                # original contig name (strip |provirus_...)
                contig = seq_name.split("|")[0]
                contig_len = int(cols[1])
                records.append(dict(sample=sample, contig=contig,
                                    prov_start=start, prov_end=end,
                                    contig_len=contig_len))
    return records

# ── get contig fasta from assembly ───────────────────────────────────────────
def find_assembly(sample):
    # look for megahit/spades final assembly
    for pattern in [
        f"{GENOMAD_BASE}/{sample}/assembly/final.contigs.fa",
        f"{GENOMAD_BASE}/{sample}/assembly/scaffolds.fasta",
        f"{GENOMAD_BASE}/{sample}/*.fasta",
        f"{GENOMAD_BASE}/{sample}/*.fa",
    ]:
        hits = list(Path("/").glob(pattern.lstrip("/")))
        if hits:
            return hits[0]
    # fallback: genomad input fasta stored alongside
    hits = list(GENOMAD_BASE.glob(f"{sample}/downloaded_genomad/*_virus_genes.tsv"))
    if hits:
        # find the fasta used as input
        d = hits[0].parent.parent
        for fa in d.glob("*.fa") or d.glob("*.fasta"):
            return fa
    return None

def extract_contig_seq(fa_path, contig_name):
    """Return full sequence of a named contig from a FASTA file."""
    seq_lines, capture = [], False
    with open(fa_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                capture = line[1:].split()[0] == contig_name
                continue
            if capture:
                seq_lines.append(line)
    return "".join(seq_lines)

def write_flanks(seq, prov_start, prov_end, contig_len, out_fa, label):
    left_start  = max(0, prov_start - 1 - FLANK_KB)
    left_end    = prov_start - 1
    right_start = prov_end
    right_end   = min(contig_len, prov_end + FLANK_KB)

    records = []
    if left_end > left_start:
        records.append((f">{label}__left_flank_{left_start+1}_{left_end}",
                        seq[left_start:left_end]))
    if right_end > right_start:
        records.append((f">{label}__right_flank_{right_start+1}_{right_end}",
                        seq[right_start:right_end]))
    if not records:
        return False
    with open(out_fa, "w") as fh:
        for hdr, s in records:
            fh.write(hdr + "\n")
            fh.write("\n".join(textwrap.wrap(s, 60)) + "\n")
    return True

# ── runners ──────────────────────────────────────────────────────────────────
def run_mef(fa, out_dir):
    fa = Path(fa).resolve()
    out_dir = Path(out_dir).resolve()
    stem = fa.stem
    cmd = ["docker", "run", "--rm",
           "-v", f"{fa.parent}:/input:ro",
           "-v", f"{out_dir}:/output",
           MEF_IMAGE,
           "mefinder", "find",
           "--complete", "/input/" + fa.name,
           "/output/" + stem]
    log = out_dir / f"{stem}_mef.log"
    with open(log, "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=lf)
    return r.returncode

def run_integron(fa, out_dir):
    fa = Path(fa).resolve()
    out_dir = Path(out_dir).resolve()
    cmd = [str(INTEGRON_BIN), str(fa),
           "--outdir", str(out_dir),
           "--cpu", "4", "--linear"]
    log = out_dir / f"{fa.stem}_integron.log"
    with open(log, "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=lf, env=ENV)
    return r.returncode

# ── main ─────────────────────────────────────────────────────────────────────
def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    records = collect_proviruses()
    print(f"Proviruses found: {len(records)}")

    # find assemblies — look in geNomad input dir
    asm_cache = {}
    for sample in sorted({r["sample"] for r in records}):
        # the geNomad input fasta is the file fed to genomad; find it
        candidates = list(GENOMAD_BASE.glob(f"{sample}/downloaded_genomad/../*.fa")) + \
                     list(GENOMAD_BASE.glob(f"{sample}/*.fa")) + \
                     list(GENOMAD_BASE.glob(f"{sample}/*.fasta"))
        # prefer megahit output naming
        for c in candidates:
            if "final.contigs" in c.name or "scaffolds" in c.name or "assembly" in c.name:
                asm_cache[sample] = c; break
        if sample not in asm_cache and candidates:
            asm_cache[sample] = candidates[0]
        if sample not in asm_cache:
            print(f"  WARNING: no assembly fasta for {sample}")

    results = []
    for i, rec in enumerate(records, 1):
        sample  = rec["sample"]
        contig  = rec["contig"]
        ps, pe  = rec["prov_start"], rec["prov_end"]
        clen    = rec["contig_len"]
        label   = f"{sample}__{contig.replace('|','__')}__prov_{ps}_{pe}"

        out_dir = OUT_BASE / label
        out_dir.mkdir(parents=True, exist_ok=True)
        fa_out  = out_dir / f"{label}_flanks.fa"

        print(f"[{i}/{len(records)}] {label}", flush=True)

        if sample not in asm_cache:
            print("   SKIP: no assembly"); continue

        seq = extract_contig_seq(asm_cache[sample], contig)
        if not seq:
            print("   SKIP: contig not found in fasta"); continue

        ok = write_flanks(seq, ps, pe, clen, fa_out, label)
        if not ok:
            print("   SKIP: flanks too short"); continue

        rc_mef = run_mef(fa_out, out_dir)
        rc_int = run_integron(fa_out, out_dir)
        print(f"   MEF rc={rc_mef}  Integron rc={rc_int}")
        results.append(dict(label=label, mef_rc=rc_mef, integron_rc=rc_int))

    # summary
    print("\n=== DONE ===")
    print(f"Processed: {len(results)}/{len(records)}")

    # collect MEF hits
    print("\n--- MEF hits ---")
    for d in OUT_BASE.iterdir():
        for csv in d.glob("*.csv"):
            with open(csv) as fh:
                lines = [l for l in fh if not l.startswith("#") and l.strip()]
            if len(lines) > 1:  # header + at least one hit
                print(f"  {csv.parent.name}: {len(lines)-1} hit(s)")

    # collect integron hits
    print("\n--- Integron hits ---")
    for d in OUT_BASE.iterdir():
        for ig in d.rglob("*.integrons"):
            with open(ig) as fh:
                content = fh.read()
            if "No Integron found" not in content and "ID_integron" in content:
                hit_lines = [l for l in content.splitlines()
                             if l and not l.startswith("#") and "ID_integron" not in l]
                if hit_lines:
                    print(f"  {d.name}: {len(hit_lines)} row(s)")

if __name__ == "__main__":
    main()
