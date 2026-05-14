#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def fasta_iter(path: Path):
    header = None
    seq_parts = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def collect_sample_dirs(root: Path):
    return sorted(root.glob("*/vRhyme/vRhyme_best_bins_fasta"))


def write_inputs(roots: list[Path], outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    proteins_out = outdir / "all_vrhyme_bins.faa"
    mapping_out = outdir / "all_vrhyme_bins.gene2genome.tsv"
    lengths_out = outdir / "all_vrhyme_bins.len_nucleotide.tsv"
    manifest_out = outdir / "all_vrhyme_bins.manifest.tsv"

    total_bins = 0
    total_proteins = 0
    total_samples = 0

    with (
        proteins_out.open("w") as faa_handle,
        mapping_out.open("w") as map_handle,
        lengths_out.open("w") as len_handle,
        manifest_out.open("w") as manifest_handle,
    ):
        map_handle.write("protein_id\tgenome_id\tkeywords\n")
        len_handle.write("genome_id\tlength\n")
        manifest_handle.write("dataset\tsample_id\tgenome_id\tfaa_path\tfasta_path\tprotein_count\tlength_bp\n")

        for root in roots:
            dataset = root.name.replace("Genomad_", "")
            sample_dirs = collect_sample_dirs(root)
            total_samples += len(sample_dirs)

            for sample_dir in sample_dirs:
                sample_id = sample_dir.parts[-3]
                for faa_path in sorted(sample_dir.glob("*.faa")):
                    bin_name = faa_path.stem
                    fasta_path = sample_dir / f"{bin_name}.fasta"
                    if not fasta_path.exists():
                        raise FileNotFoundError(f"Missing nucleotide FASTA for {faa_path}")

                    genome_id = f"{sample_id}|{bin_name}"
                    nuc_length = 0
                    for _, seq in fasta_iter(fasta_path):
                        nuc_length += len(seq)
                    len_handle.write(f"{genome_id}\t{nuc_length}\n")

                    protein_count = 0
                    for protein_header, seq in fasta_iter(faa_path):
                        protein_count += 1
                        protein_id = f"{genome_id}|{protein_header}"
                        faa_handle.write(f">{protein_id}\n{seq}\n")
                        map_handle.write(f"{protein_id}\t{genome_id}\tNone\n")

                    manifest_handle.write(
                        f"{dataset}\t{sample_id}\t{genome_id}\t{faa_path}\t{fasta_path}\t{protein_count}\t{nuc_length}\n"
                    )
                    total_bins += 1
                    total_proteins += protein_count

    print(f"sample_dirs={total_samples}")
    print(f"total_bins={total_bins}")
    print(f"total_proteins={total_proteins}")
    print(f"proteins={proteins_out}")
    print(f"gene2genome={mapping_out}")
    print(f"len_nucleotide={lengths_out}")
    print(f"manifest={manifest_out}")


def main():
    parser = argparse.ArgumentParser(description="Build vConTACT3 protein-mode inputs from vRhyme viral bin FASTAs.")
    parser.add_argument(
        "--root",
        action="append",
        required=True,
        help="Root directory containing per-sample vRhyme/vRhyme_best_bins_fasta folders. Repeat for multiple datasets.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory for combined vConTACT3 inputs.")
    args = parser.parse_args()

    roots = [Path(p).resolve() for p in args.root]
    outdir = Path(args.outdir).resolve()
    write_inputs(roots, outdir)


if __name__ == "__main__":
    main()
