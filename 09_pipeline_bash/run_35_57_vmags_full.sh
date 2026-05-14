#!/usr/bin/env bash
set -euo pipefail

THREADS="${1:-22}"
START_SAMPLE="${2:-1035}"
END_SAMPLE="${3:-1057}"
DELETE_FASTQ_ON_SUCCESS="${DELETE_FASTQ_ON_SUCCESS:-0}"
ALLOW_FASTQ_DELETE="${ALLOW_FASTQ_DELETE:-0}"
FORCE_REPROCESS_COMPLETED="${FORCE_REPROCESS_COMPLETED:-0}"

# USER: set your base path here

ASSEMBLY_BASE="/mnt/nvme/RECOVERY_m2/Metagenome/Tati/Assembly"
VIBRANT_BASE="$ASSEMBLY_BASE/Vibrant"
# USER: set your base path here
FASTQ_NEW_BASE="/mnt/nvme2/Tati/FastQ"
FASTQ_FALLBACK_BASE="$ASSEMBLY_BASE/FastQ"

RUN_BASE="$ASSEMBLY_BASE/vRhyme_batch_35_57"
ANOT_ROOT="$ASSEMBLY_BASE/vRhyme_test_P19109_1035/Anotation"
# USER: set your base path here
CHECKV_DB="/mnt/nvme/Taxonomia/Proteins/Checkv/DB_CLEAN/checkv-db-v1.5"
# USER: set your base path here
PHAROKKA_DB="/mnt/nvme/Taxonomia/Proteins/Pharokka_DB"
# USER: set your base path here
PHABOX_DB="/mnt/nvme/Taxonomia/phabox_db_v2_1"
# USER: set your base path here
BLAST_DB="/mnt/nvme/blastdb_viruses/ref_viruses_rep_genomes"
PHAROKKA_TIMEOUT_SECS="${PHAROKKA_TIMEOUT_SECS:-5400}"
PHAROKKA_RETRY_FAILED="${PHAROKKA_RETRY_FAILED:-0}"
PHAROKKA_OPTIONAL_ON_STALL="${PHAROKKA_OPTIONAL_ON_STALL:-1}"
SAMPLE_MAX_SECONDS="${SAMPLE_MAX_SECONDS:-18000}"
SAMPLE_MAX_PHAROKKA_TIMEOUTS="${SAMPLE_MAX_PHAROKKA_TIMEOUTS:-5}"
PHAGENUS_TIMEOUT_SECS="${PHAGENUS_TIMEOUT_SECS:-3600}"
PHAGENUS_RETRY_FAILED="${PHAGENUS_RETRY_FAILED:-0}"
PHAGENUS_ROOT="/home/rafael/Documentos/Phage_genus/phagenus-main"
# USER: set your base path here
PHAGENUS_PY="/mnt/nvme/conda_envs/phagenus/bin/python"
PHAGENUS_SCRIPT="$PHAGENUS_ROOT/phagenus.py"
CONSENSUS_PY="$ASSEMBLY_BASE/taxonomy_4tools_consensus.py"
VMAG_CONFIDENCE_PY="$ASSEMBLY_BASE/vmag_confidence_score.py"
VMAG_KEEP_MIN_CONF="${VMAG_KEEP_MIN_CONF:-Medium}"
RUN_GENOMAD="${RUN_GENOMAD:-1}"
# USER: set your base path here
GENOMAD_BIN="/mnt/nvme/conda_envs/mvip/bin/genomad"
# USER: set your base path here
GENOMAD_DB="/mnt/nvme/Taxonomia/Proteins/genomad_db/genomad_db"
GENOMAD_EXTRACT_PY="$ASSEMBLY_BASE/genomad_taxonomy_extract.py"
GENOMAD_INPUT_SET="${GENOMAD_INPUT_SET:-selected}" # selected|merged
RUN_TAXMYPHAGE="${RUN_TAXMYPHAGE:-1}"
# USER: set your base path here
TAXMYPHAGE_BIN="/mnt/nvme/conda_envs/taxmyphage/bin/taxmyphage"
# USER: set your base path here
TAXMYPHAGE_DB="/mnt/nvme/Taxonomia/Proteins/taxmyphage_db"
TAXMYPHAGE_EXTRACT_PY="$ASSEMBLY_BASE/taxmyphage_extract.py"
TAXMYPHAGE_TIMEOUT_SECS="${TAXMYPHAGE_TIMEOUT_SECS:-10800}"
TAXMYPHAGE_THREADS="${TAXMYPHAGE_THREADS:-$THREADS}"
VMAG_CONCAT_EXPECTED_PY="$ASSEMBLY_BASE/vmag_concat_checkv_expected.py"
RUN_CONCAT_CHECKV_EXPECTED="${RUN_CONCAT_CHECKV_EXPECTED:-1}"
CONCAT_GAP_NS="${CONCAT_GAP_NS:-0}"
CHECKV_CONCAT_THREADS="${CHECKV_CONCAT_THREADS:-8}"
NCBI_CACHE_JSON="$RUN_BASE/ncbi_genome_size_cache.json"
NCBI_EMAIL="${NCBI_EMAIL:-}"
NCBI_API_KEY="${NCBI_API_KEY:-}"
RUN_VIRIDIC="${RUN_VIRIDIC:-1}"
VIRIDIC_PY="$ASSEMBLY_BASE/vmag_viridic.py"
# USER: set your base path here
VIRIDIC_RSCRIPT="/mnt/nvme/conda_envs/viridic/bin/Rscript"
# USER: set your base path here
VIRIDIC_SCRIPTS_DIR="/mnt/nvme/tools/VIRIDIC/VIRIDIC/stand_alone/viridic_scripts"
VIRIDIC_TIMEOUT_SECS="${VIRIDIC_TIMEOUT_SECS:-5400}"
VIRIDIC_THREADS="${VIRIDIC_THREADS:-$THREADS}"

# BLAST thresholds used in 4-tool taxonomy consensus
BLAST_TOPN=10
BLAST_MAX_TARGET_SEQS=50
BLAST_EVALUE_MAX=1e-20
BLAST_PIDENT_MIN=85.0
BLAST_MIN_ALN_BP=500
BLAST_MIN_ALN_FRAC_OF_CONTIG=0.05
BLAST_MIN_VALID_HITS=3
BLAST_MAJORITY_FRAC=0.60
BLAST_ANI_PROXY_PIDENT_MIN="${BLAST_ANI_PROXY_PIDENT_MIN:-95.0}"
BLAST_ANI_PROXY_QCOV_MIN="${BLAST_ANI_PROXY_QCOV_MIN:-0.85}"
REP_CONTIG_MIN_BP="${REP_CONTIG_MIN_BP:-10000}"

# USER: set your base path here

VRHYME_BIN="/mnt/nvme/conda_envs/mvip/bin/vRhyme"
# USER: set your base path here
CHECKV_BIN="/mnt/nvme/conda_envs/checkv_env/bin/checkv"
# USER: set your base path here
PHABOX_BIN="/mnt/nvme/conda_envs/phabox-env/bin/phabox2"
# USER: set your base path here
BLASTN_BIN="/mnt/nvme/conda_envs/bioinfo/bin/blastn"
# USER: set your base path here
PHAROKKA_BIN="/mnt/nvme/conda_envs/pharokka_env/bin/pharokka.py"
SAMTOOLS_BIN="/usr/bin/samtools"

mkdir -p "$RUN_BASE" "$ANOT_ROOT"
MASTER_SUMMARY="$RUN_BASE/master_summary.tsv"
if [ ! -s "$MASTER_SUMMARY" ]; then
  echo -e "sample\tstatus\tn_bins\tn_clean_bins\tmerged_clean_contigs\tnote" > "$MASTER_SUMMARY"
fi
OVERALL_BIN="$RUN_BASE/taxonomy_4tools_bin_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_CONTIG="$RUN_BASE/taxonomy_4tools_contig_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_SUMMARY="$RUN_BASE/taxonomy_4tools_overall_${START_SAMPLE}_${END_SAMPLE}.summary.tsv"
OVERALL_CONF="$RUN_BASE/vmag_confidence_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_CONF_KEEP="$RUN_BASE/vmag_confidence_overall_${START_SAMPLE}_${END_SAMPLE}_${VMAG_KEEP_MIN_CONF,,}plus.tsv"
OVERALL_CONF_SUMMARY="$RUN_BASE/vmag_confidence_overall_${START_SAMPLE}_${END_SAMPLE}.summary.tsv"
OVERALL_GENOMAD_BIN="$RUN_BASE/genomad_bin_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_GENOMAD_CONTIG="$RUN_BASE/genomad_contig_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_GENOMAD_SUMMARY="$RUN_BASE/genomad_overall_${START_SAMPLE}_${END_SAMPLE}.summary.tsv"
OVERALL_CONCAT_EXPECTED="$RUN_BASE/vmag_concat_checkv_expected_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_CONCAT_EXPECTED_SUMMARY="$RUN_BASE/vmag_concat_checkv_expected_overall_${START_SAMPLE}_${END_SAMPLE}.summary.tsv"
OVERALL_VIRIDIC_BIN="$RUN_BASE/viridic_bin_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_VIRIDIC_SUMMARY="$RUN_BASE/viridic_overall_${START_SAMPLE}_${END_SAMPLE}.summary.tsv"
OVERALL_TAXMYPHAGE_BIN="$RUN_BASE/taxmyphage_bin_overall_${START_SAMPLE}_${END_SAMPLE}.tsv"
OVERALL_TAXMYPHAGE_SUMMARY="$RUN_BASE/taxmyphage_overall_${START_SAMPLE}_${END_SAMPLE}.summary.tsv"

find_fastq_pair() {
  local sample="$1"
  local r1 r2
  r1=$(ls "$FASTQ_NEW_BASE"/P19109_"$sample"_*_R1_001.fastq.gz 2>/dev/null | head -n1 || true)
  r2=$(ls "$FASTQ_NEW_BASE"/P19109_"$sample"_*_R2_001.fastq.gz 2>/dev/null | head -n1 || true)
  if [ -n "$r1" ] && [ -n "$r2" ]; then
    echo "$r1|$r2"
    return
  fi
  r1=$(ls "$FASTQ_FALLBACK_BASE"/P19109_"$sample"_*_R1_001.fastq.gz 2>/dev/null | head -n1 || true)
  r2=$(ls "$FASTQ_FALLBACK_BASE"/P19109_"$sample"_*_R2_001.fastq.gz 2>/dev/null | head -n1 || true)
  if [ -n "$r1" ] && [ -n "$r2" ]; then
    echo "$r1|$r2"
    return
  fi
  echo "|"
}

for sample in $(seq "$START_SAMPLE" "$END_SAMPLE"); do
  s="P19109_${sample}"
  echo "================ $s ================"

  sample_dir="$RUN_BASE/$s"
  vr_out="$sample_dir/vRhyme"
  checkv_root="$sample_dir/CheckV_bins"
  bins_dir="$vr_out/vRhyme_best_bins_fasta"
  sample_anot="$ANOT_ROOT/$sample"
  sample_tax="$sample_anot/_sample_level"
  mkdir -p "$sample_dir" "$checkv_root" "$sample_anot" "$sample_tax"
  sample_done_marker="$sample_tax/.pipeline_complete.ok"
  merged="$sample_tax/${s}.vMAGs.clean.merged.fasta"
  phabox_out="$sample_tax/phabox"
  phagenus_out="$sample_tax/phagenus/prediction_output.csv"
  blast_out="$sample_tax/${s}.vmags.blastn.tsv"
  tax_bin_tsv="$sample_tax/taxonomy_4tools_bin_${sample}.tsv"
  tax_contig_tsv="$sample_tax/taxonomy_4tools_contig_${sample}.tsv"
  conf_tsv="$sample_tax/vmag_confidence_${sample}.tsv"
  conf_keep_tsv="$sample_tax/vmag_confidence_filtered_${sample}_${VMAG_KEEP_MIN_CONF,,}plus.tsv"
  conf_selected_merged="$sample_tax/${s}.vMAGs.selected_${VMAG_KEEP_MIN_CONF,,}plus.merged.fasta"
  genomad_bin_tsv="$sample_tax/genomad_bin_${sample}.tsv"
  genomad_contig_tsv="$sample_tax/genomad_contig_${sample}.tsv"
  concat_expected_tsv="$sample_tax/vmag_concat_checkv_expected_${sample}.tsv"
  concat_expected_summary_tsv="$sample_tax/vmag_concat_checkv_expected_${sample}.summary.tsv"
  viridic_bin_tsv="$sample_tax/viridic_bin_${sample}.tsv"
  viridic_summary_tsv="$sample_tax/viridic_summary_${sample}.tsv"
  taxmyphage_dir="$sample_tax/taxmyphage"
  taxmyphage_raw_summary="$taxmyphage_dir/Summary_taxonomy.tsv"
  taxmyphage_bin_tsv="$sample_tax/taxmyphage_bin_${sample}.tsv"
  taxmyphage_summary_tsv="$sample_tax/taxmyphage_summary_${sample}.tsv"

  if [ "$FORCE_REPROCESS_COMPLETED" != "1" ]; then
    if [ -s "$sample_done_marker" ]; then
      echo "[$s] checkpoint: sample completo, pulando."
      continue
    fi
    if grep -q "^${sample}[[:space:]]\\+ok[[:space:]]" "$MASTER_SUMMARY" 2>/dev/null \
      && [ -s "$merged" ] \
      && [ -s "$phabox_out/final_prediction/phagcn_prediction.tsv" ] \
      && { [ -s "$phagenus_out" ] || [ -s "$sample_tax/phagenus/.phagenus_timeout" ] || [ -s "$sample_tax/phagenus/.phagenus_failed" ]; } \
      && [ -s "$blast_out" ] \
      && [ -s "$tax_bin_tsv" ] \
      && [ -s "$tax_contig_tsv" ] \
      && [ -s "$conf_tsv" ] \
      && [ -s "$conf_keep_tsv" ] \
      && [ -s "$conf_selected_merged" ] \
      && { [ "$RUN_CONCAT_CHECKV_EXPECTED" != "1" ] || { [ -s "$concat_expected_tsv" ] && [ -s "$concat_expected_summary_tsv" ]; }; } \
      && { [ "$RUN_VIRIDIC" != "1" ] || { [ -s "$viridic_bin_tsv" ] && [ -s "$viridic_summary_tsv" ]; }; } \
      && { [ "$RUN_TAXMYPHAGE" != "1" ] || { [ -s "$taxmyphage_bin_tsv" ] && [ -s "$taxmyphage_summary_tsv" ]; }; } \
      && { [ "$RUN_GENOMAD" != "1" ] || { [ -s "$genomad_bin_tsv" ] && [ -s "$genomad_contig_tsv" ]; }; }; then
      printf "timestamp\t%s\nsource\tmaster_summary+outputs\n" "$(date -Is)" > "$sample_done_marker"
      echo "[$s] checkpoint rápido: outputs já completos, pulando."
      continue
    fi
  fi

  sample_start_ts=$(date +%s)
  sample_stalled=0
  sample_stall_reason=""
  pharokka_timeout_count=0
  pharokka_incomplete=0
  pharokka_incomplete_reason=""

  phages="$VIBRANT_BASE/VIBRANT_${s}/VIBRANT_phages_${s}/${s}.phages_combined.fna"
  if [ ! -s "$phages" ]; then
    echo -e "$sample\tskipped\t0\t0\t0\tmissing_phages_combined" >> "$MASTER_SUMMARY"
    continue
  fi

  pair=$(find_fastq_pair "$sample")
  r1="${pair%%|*}"
  r2="${pair##*|}"
  if [ -z "$r1" ] || [ -z "$r2" ] || [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
    echo -e "$sample\twaiting_fastq\t0\t0\t0\tfastq_missing_or_incomplete" >> "$MASTER_SUMMARY"
    continue
  fi

  # Normalize read names for vRhyme pair detection.
  reads_dir="$sample_dir/reads"
  mkdir -p "$reads_dir"
  r1_norm="$reads_dir/${s}_R1.fastq.gz"
  r2_norm="$reads_dir/${s}_R2.fastq.gz"
  ln -sfn "$r1" "$r1_norm"
  ln -sfn "$r2" "$r2_norm"

  vr_membership=$(ls "$vr_out"/vRhyme_best_bins.*.membership.tsv 2>/dev/null | head -n1 || true)
  if [ -z "$vr_membership" ] || [ ! -s "$vr_membership" ]; then
    # vRhyme aborts if output dir already exists; remove stale/incomplete output.
    [ -d "$vr_out" ] && rm -rf "$vr_out"
    if ! PATH="/mnt/nvme/conda_envs/mvip/bin:$PATH" "$VRHYME_BIN" \
      -i "$phages" -r "$r1_norm" "$r2_norm" -t "$THREADS" -o "$vr_out" > "$sample_dir/vrhyme.log" 2>&1; then
      echo -e "$sample\tfailed\t0\t0\t0\tvrhyme_failed" >> "$MASTER_SUMMARY"
      continue
    fi
  fi

  if [ ! -d "$bins_dir" ]; then
    echo -e "$sample\tfailed\t0\t0\t0\tvrhyme_no_bins_dir" >> "$MASTER_SUMMARY"
    continue
  fi

  # CheckV + cleaning
  clean_count=0
  total_bins=0
  for fasta in "$bins_dir"/vRhyme_bin_*.fasta; do
    [ -e "$fasta" ] || continue
    [[ "$fasta" == *.clean.fasta ]] && continue
    total_bins=$((total_bins+1))
    bn=$(basename "$fasta" .fasta)
    qdir="$checkv_root/$bn"
    tsv="$qdir/quality_summary.tsv"
    mkdir -p "$qdir"

    clean="$bins_dir/${bn}.clean.fasta"
    if [ -s "$tsv" ] && [ -s "$clean" ] && grep -q '^>' "$clean"; then
      clean_count=$((clean_count+1))
      continue
    fi

    if [ ! -s "$tsv" ]; then
      # USER: set your base path here
      PATH="/mnt/nvme/conda_envs/checkv_env/bin:$PATH" "$CHECKV_BIN" end_to_end "$fasta" "$qdir" -d "$CHECKV_DB" -t 8 > "$qdir/run.log" 2>&1 || true
    fi

    if [ ! -s "$tsv" ]; then
      continue
    fi

    awk -F'\t' 'NR>1 && ($6+0)>0 {print $1}' "$tsv" > "$qdir/keep.ids"
    awk 'NR==FNR{keep[$1]=1; next} /^>/{id=substr($1,2); flag=(id in keep)} flag{print}' "$qdir/keep.ids" "$fasta" > "$clean"

    if grep -q '^>' "$clean"; then
      clean_count=$((clean_count+1))
    fi
  done

  : > "$merged"
  for c in "$bins_dir"/vRhyme_bin_*.clean.fasta; do
    [ -s "$c" ] || continue
    grep -q '^>' "$c" || continue
    cat "$c" >> "$merged"
  done
  merged_contigs=$(grep -c '^>' "$merged" 2>/dev/null || true)

  if [ "$merged_contigs" -eq 0 ]; then
    echo -e "$sample\tfailed\t$total_bins\t$clean_count\t0\tno_clean_merged_contigs" >> "$MASTER_SUMMARY"
    continue
  fi

  # PhaBOX
  if [ ! -s "$phabox_out/final_prediction/phagcn_prediction.tsv" ]; then
    mkdir -p "$phabox_out"
    # USER: set your base path here
    PATH="/mnt/nvme/conda_envs/phabox-env/bin:$PATH" "$PHABOX_BIN" --task end_to_end --dbdir "$PHABOX_DB" --contigs "$merged" --outpth "$phabox_out" --threads "$THREADS" > "$sample_tax/phabox.log" 2>&1 || true
  fi

  # PhaGenus
  if [ ! -s "$sample_tax/phagenus/prediction_output.csv" ]; then
    mkdir -p "$sample_tax/phagenus"
    phagenus_timeout_marker="$sample_tax/phagenus/.phagenus_timeout"
    phagenus_failed_marker="$sample_tax/phagenus/.phagenus_failed"
    if [ "$PHAGENUS_RETRY_FAILED" != "1" ] && { [ -s "$phagenus_timeout_marker" ] || [ -s "$phagenus_failed_marker" ]; }; then
      :
    else
      phagenus_rc=0
      (
        cd "$PHAGENUS_ROOT"
        # USER: set your base path here
        PATH="/mnt/nvme/conda_envs/mvip/bin:/mnt/nvme/conda_envs/phagenus/bin:$PATH" timeout --signal=TERM --kill-after=120 "${PHAGENUS_TIMEOUT_SECS}" \
          "$PHAGENUS_PY" phagenus.py \
            --contigs "$merged" --midfolder "$sample_tax/phagenus" --out prediction_output.csv --sim high --threads "$THREADS"
      ) > "$sample_tax/phagenus.log" 2>&1 || phagenus_rc=$?
      if [ -s "$sample_tax/phagenus/prediction_output.csv" ]; then
        rm -f "$phagenus_timeout_marker" "$phagenus_failed_marker"
      elif [ "$phagenus_rc" -eq 124 ] || [ "$phagenus_rc" -eq 137 ]; then
        printf "timestamp\t%s\nexit_code\t%s\nreason\ttimeout\n" "$(date -Is)" "$phagenus_rc" > "$phagenus_timeout_marker"
        rm -f "$phagenus_failed_marker"
      elif [ "$phagenus_rc" -ne 0 ]; then
        printf "timestamp\t%s\nexit_code\t%s\nreason\tfailed\n" "$(date -Is)" "$phagenus_rc" > "$phagenus_failed_marker"
        rm -f "$phagenus_timeout_marker"
      fi
    fi
  fi

  # BLASTn
  if [ ! -s "$blast_out" ]; then
    # USER: set your base path here
    BLASTDB="/mnt/nvme/blastdb_viruses" "$BLASTN_BIN" \
      -query "$merged" -db "$BLAST_DB" -out "$blast_out" \
      -outfmt "6 qseqid sseqid pident length evalue bitscore staxids sscinames qlen" \
      -max_target_seqs 100 -num_threads "$THREADS" > "$sample_tax/blastn.log" 2>&1 || true
  fi

  # Pharokka per clean bin
  for clean in "$bins_dir"/vRhyme_bin_*.clean.fasta; do
    now_ts=$(date +%s)
    sample_elapsed=$((now_ts - sample_start_ts))
    if [ "$sample_elapsed" -ge "$SAMPLE_MAX_SECONDS" ]; then
      sample_stalled=1
      sample_stall_reason="sample_elapsed_limit_${SAMPLE_MAX_SECONDS}s"
      pharokka_incomplete=1
      pharokka_incomplete_reason="$sample_stall_reason"
      break
    fi

    [ -s "$clean" ] || continue
    grep -q '^>' "$clean" || continue
    bn=$(basename "$clean" .clean.fasta)
    outdir="$sample_anot/$bn"
    mkdir -p "$outdir"
    gbk="$outdir/${bn}.gbk"
    mash="$outdir/${bn}_top_hits_mash_inphared.tsv"
    timeout_marker="$outdir/.pharokka_timeout"
    failed_marker="$outdir/.pharokka_failed"

    if [ -s "$gbk" ] && [ -s "$mash" ]; then
      continue
    fi
    if [ "$PHAROKKA_RETRY_FAILED" != "1" ] && [ ! -s "$gbk" ] && [ ! -s "$mash" ] && [ -s "$outdir/${bn}.pharokka.log" ] && [ ! -s "$timeout_marker" ] && [ ! -s "$failed_marker" ]; then
      printf "timestamp\t%s\nexit_code\tNA\nreason\tprevious_incomplete_attempt\n" "$(date -Is)" > "$failed_marker"
      continue
    fi
    if [ "$PHAROKKA_RETRY_FAILED" != "1" ] && { [ -s "$timeout_marker" ] || [ -s "$failed_marker" ]; }; then
      if [ -s "$timeout_marker" ]; then
        pharokka_timeout_count=$((pharokka_timeout_count + 1))
      fi
      continue
    fi

    ncontigs=$(grep -c '^>' "$clean" || true)
    pharokka_rc=0
    if [ "$ncontigs" -gt 1 ]; then
      # USER: set your base path here
      PATH="/mnt/nvme/conda_envs/pharokka_env/bin:$PATH" timeout --signal=TERM --kill-after=120 "${PHAROKKA_TIMEOUT_SECS}" \
        "$PHAROKKA_BIN" -i "$clean" -o "$outdir" -d "$PHAROKKA_DB" -t "$THREADS" --meta --skip_extra_annotations --fast --prefix "$bn" --force \
        > "$outdir/${bn}.pharokka.log" 2>&1 || pharokka_rc=$?
    else
      # USER: set your base path here
      PATH="/mnt/nvme/conda_envs/pharokka_env/bin:$PATH" timeout --signal=TERM --kill-after=120 "${PHAROKKA_TIMEOUT_SECS}" \
        "$PHAROKKA_BIN" -i "$clean" -o "$outdir" -d "$PHAROKKA_DB" -t "$THREADS" --skip_extra_annotations --fast --prefix "$bn" --force \
        > "$outdir/${bn}.pharokka.log" 2>&1 || pharokka_rc=$?
    fi
    [ -s "$outdir/pharokka.gbk" ] && cp -f "$outdir/pharokka.gbk" "$outdir/${bn}.gbk"
    [ -s "$outdir/pharokka.gff" ] && cp -f "$outdir/pharokka.gff" "$outdir/${bn}.gff"
    [ -s "$outdir/pharokka_top_hits_mash_inphared.tsv" ] && cp -f "$outdir/pharokka_top_hits_mash_inphared.tsv" "$outdir/${bn}_top_hits_mash_inphared.tsv"
    if [ -s "$gbk" ] && [ -s "$mash" ]; then
      rm -f "$timeout_marker" "$failed_marker"
    elif [ "$pharokka_rc" -eq 124 ] || [ "$pharokka_rc" -eq 137 ]; then
      printf "timestamp\t%s\nexit_code\t%s\nreason\ttimeout\n" "$(date -Is)" "$pharokka_rc" > "$timeout_marker"
      rm -f "$failed_marker"
      pharokka_timeout_count=$((pharokka_timeout_count + 1))
    elif [ "$pharokka_rc" -ne 0 ]; then
      printf "timestamp\t%s\nexit_code\t%s\nreason\tfailed\n" "$(date -Is)" "$pharokka_rc" > "$failed_marker"
      rm -f "$timeout_marker"
    elif [ ! -s "$gbk" ] || [ ! -s "$mash" ]; then
      printf "timestamp\t%s\nexit_code\t0\nreason\tmissing_expected_outputs\n" "$(date -Is)" > "$failed_marker"
      rm -f "$timeout_marker"
    fi

    if [ "$pharokka_timeout_count" -ge "$SAMPLE_MAX_PHAROKKA_TIMEOUTS" ]; then
      sample_stalled=1
      sample_stall_reason="pharokka_timeouts_${pharokka_timeout_count}"
      pharokka_incomplete=1
      pharokka_incomplete_reason="$sample_stall_reason"
      break
    fi
  done

  if [ "$sample_stalled" -eq 1 ]; then
    if [ "$PHAROKKA_OPTIONAL_ON_STALL" = "1" ] && [[ "$sample_stall_reason" == pharokka_timeouts_* || "$sample_stall_reason" == sample_elapsed_limit_* ]]; then
      sample_stalled=0
      pharokka_incomplete=1
      pharokka_incomplete_reason="$sample_stall_reason"
      echo "[$s] Pharokka incompleto (${sample_stall_reason}); seguindo com consenso sem bloquear amostra."
    else
      echo -e "$sample\tstalled\t$total_bins\t$clean_count\t$merged_contigs\t$sample_stall_reason" >> "$MASTER_SUMMARY"
      rm -f "$sample_done_marker"
      continue
    fi
  fi

  # 4-tool taxonomy consensus per sample (Pharokka + PhaGenus + PhaBOX + BLASTn)
  if [ -x "$CONSENSUS_PY" ]; then
    python3 "$CONSENSUS_PY" \
      --sample "$sample" \
      --annotation-root "$sample_anot" \
      --sample-tag "$s" \
      --blast-topn "$BLAST_TOPN" \
      --blast-max-target-seqs "$BLAST_MAX_TARGET_SEQS" \
      --blast-evalue-max "$BLAST_EVALUE_MAX" \
      --blast-pident-min "$BLAST_PIDENT_MIN" \
      --blast-min-aln-bp "$BLAST_MIN_ALN_BP" \
      --blast-min-aln-frac-of-contig "$BLAST_MIN_ALN_FRAC_OF_CONTIG" \
      --blast-min-valid-hits "$BLAST_MIN_VALID_HITS" \
      --blast-majority-frac "$BLAST_MAJORITY_FRAC" \
      --blast-ani-proxy-pident-min "$BLAST_ANI_PROXY_PIDENT_MIN" \
      --blast-ani-proxy-qcov-min "$BLAST_ANI_PROXY_QCOV_MIN" \
      --rep-contig-min-bp "$REP_CONTIG_MIN_BP" \
      > "$sample_tax/taxonomy_4tools.log" 2>&1 || true
  fi

  # vMAG confidence scoring and filtering (default keep: Medium+High)
  if [ -x "$VMAG_CONFIDENCE_PY" ]; then
    python3 "$VMAG_CONFIDENCE_PY" \
      --sample "$sample" \
      --annotation-root "$sample_anot" \
      --checkv-root "$checkv_root" \
      --bins-dir "$bins_dir" \
      --sample-tag "$s" \
      --keep-min-confidence "$VMAG_KEEP_MIN_CONF" \
      > "$sample_tax/vmag_confidence.log" 2>&1 || true
  fi

  # Concatenated vMAG-by-bin CheckV + NCBI expected genome size (species/genus/family)
  if [ "$RUN_CONCAT_CHECKV_EXPECTED" = "1" ] && [ -x "$VMAG_CONCAT_EXPECTED_PY" ]; then
    python3 "$VMAG_CONCAT_EXPECTED_PY" \
      --sample "$sample" \
      --sample-tag "$s" \
      --bins-dir "$bins_dir" \
      --sample-level "$sample_tax" \
      --checkv-bin "$CHECKV_BIN" \
      --checkv-db "$CHECKV_DB" \
      --threads "$CHECKV_CONCAT_THREADS" \
      --gap-n "$CONCAT_GAP_NS" \
      --cache-json "$NCBI_CACHE_JSON" \
      --ncbi-email "$NCBI_EMAIL" \
      --ncbi-api-key "$NCBI_API_KEY" \
      > "$sample_tax/vmag_concat_checkv_expected.log" 2>&1 || true
  fi

  # TaxMyPhage on concatenated bins
  concat_mode="no_gap"
  if [ "$CONCAT_GAP_NS" -gt 0 ]; then
    concat_mode="n_gap_${CONCAT_GAP_NS}"
  fi
  concat_input="$sample_tax/concat_checkv_expected/${s}.vmags.concat_${concat_mode}.fasta"

  if [ "$RUN_TAXMYPHAGE" = "1" ] && [ -x "$TAXMYPHAGE_BIN" ] && [ -d "$TAXMYPHAGE_DB" ]; then
    if [ -s "$concat_input" ] && grep -q '^>' "$concat_input"; then
      if [ ! -s "$taxmyphage_raw_summary" ]; then
        mkdir -p "$taxmyphage_dir"
        # USER: set your base path here
        PATH="/mnt/nvme/conda_envs/taxmyphage/bin:$PATH" timeout --signal=TERM --kill-after=120 "${TAXMYPHAGE_TIMEOUT_SECS}" \
          "$TAXMYPHAGE_BIN" run \
            -i "$concat_input" \
            -o "$taxmyphage_dir" \
            -t "$TAXMYPHAGE_THREADS" \
            -db "$TAXMYPHAGE_DB" \
            --no-figures \
            --mash /mnt/nvme/conda_envs/taxmyphage/bin/mash \
            --blastdbcmd /mnt/nvme/conda_envs/taxmyphage/bin/blastdbcmd \
            --blastn /mnt/nvme/conda_envs/taxmyphage/bin/blastn \
            --makeblastdb /mnt/nvme/conda_envs/taxmyphage/bin/makeblastdb \
          > "$sample_tax/taxmyphage.log" 2>&1 || true
      fi

      if [ -x "$TAXMYPHAGE_EXTRACT_PY" ] && [ -s "$taxmyphage_raw_summary" ]; then
        python3 "$TAXMYPHAGE_EXTRACT_PY" \
          --sample "$sample" \
          --sample-level "$sample_tax" \
          > "$sample_tax/taxmyphage_extract.log" 2>&1 || true
      fi
    fi
  fi

  # VIRIDIC clustering support on concatenated vMAGs (no taxonomy override)
  if [ "$RUN_VIRIDIC" = "1" ] && [ -x "$VIRIDIC_PY" ] && [ -x "$VIRIDIC_RSCRIPT" ] && [ -d "$VIRIDIC_SCRIPTS_DIR" ]; then
    # USER: set your base path here
    PATH="/mnt/nvme/conda_envs/viridic/bin:$PATH" timeout --signal=TERM --kill-after=120 "${VIRIDIC_TIMEOUT_SECS}" \
      python3 "$VIRIDIC_PY" \
        --sample "$sample" \
        --sample-tag "$s" \
        --sample-level "$sample_tax" \
        --rscript-bin "$VIRIDIC_RSCRIPT" \
        --viridic-scripts-dir "$VIRIDIC_SCRIPTS_DIR" \
        --threads "$VIRIDIC_THREADS" \
        --gap-n "$CONCAT_GAP_NS" \
      > "$sample_tax/viridic.log" 2>&1 || true
  fi

  # geNomad annotation/classification/taxonomy
  if [ "$RUN_GENOMAD" = "1" ] && [ -x "$GENOMAD_BIN" ] && [ -d "$GENOMAD_DB" ]; then
    if [ "$GENOMAD_INPUT_SET" = "selected" ]; then
      genomad_input="$conf_selected_merged"
      genomad_out="$sample_tax/genomad_selected_${VMAG_KEEP_MIN_CONF,,}plus"
    else
      genomad_input="$merged"
      genomad_out="$sample_tax/genomad_merged"
    fi
    genomad_done="$genomad_out/.genomad_complete.ok"

    if [ -s "$genomad_input" ] && grep -q '^>' "$genomad_input"; then
      if [ ! -s "$genomad_done" ]; then
        # USER: set your base path here
        PATH="/mnt/nvme/conda_envs/mvip/bin:$PATH" "$GENOMAD_BIN" end-to-end \
          "$genomad_input" "$genomad_out" "$GENOMAD_DB" -t "$THREADS" --conservative-taxonomy \
          > "$sample_tax/genomad.log" 2>&1 || true

        if find "$genomad_out" -type f -name '*_virus_summary.tsv' | grep -q .; then
          printf "timestamp\t%s\ninput\t%s\nout\t%s\n" "$(date -Is)" "$genomad_input" "$genomad_out" > "$genomad_done"
        fi
      fi

      if [ -x "$GENOMAD_EXTRACT_PY" ] && find "$genomad_out" -type f -name '*_virus_summary.tsv' | grep -q .; then
        python3 "$GENOMAD_EXTRACT_PY" \
          --sample "$sample" \
          --genomad-out "$genomad_out" \
          --sample-level "$sample_tax" \
          > "$sample_tax/genomad_extract.log" 2>&1 || true
      fi
    fi
  fi

  note="completed_sample_loop"
  if [ "$pharokka_incomplete" -eq 1 ]; then
    note="${note};pharokka_incomplete:${pharokka_incomplete_reason:-unknown}"
  fi
  if [ "$DELETE_FASTQ_ON_SUCCESS" = "1" ] && [ "$ALLOW_FASTQ_DELETE" = "1" ]; then
    if [[ "$r1" == "$FASTQ_NEW_BASE/"* ]] && [[ "$r2" == "$FASTQ_NEW_BASE/"* ]]; then
      rm -f "$r1" "$r2" || true
      note="${note};fastq_deleted"
    else
      note="${note};fastq_kept_non_new_base"
    fi
  elif [ "$DELETE_FASTQ_ON_SUCCESS" = "1" ] && [ "$ALLOW_FASTQ_DELETE" != "1" ]; then
    note="${note};fastq_delete_blocked"
  fi

  phagenus_ok=0
  if [ -s "$phagenus_out" ] || [ -s "$sample_tax/phagenus/.phagenus_timeout" ] || [ -s "$sample_tax/phagenus/.phagenus_failed" ]; then
    phagenus_ok=1
  fi
  genomad_ok=1
  if [ "$RUN_GENOMAD" = "1" ] && [ ! -s "$genomad_bin_tsv" ] && [ ! -s "$genomad_contig_tsv" ]; then
    genomad_ok=0
  fi

  if [ -s "$merged" ] \
    && [ -s "$phabox_out/final_prediction/phagcn_prediction.tsv" ] \
    && [ "$phagenus_ok" -eq 1 ] \
    && [ -s "$blast_out" ] \
    && [ -s "$tax_bin_tsv" ] \
    && [ -s "$tax_contig_tsv" ] \
    && [ -s "$conf_tsv" ] \
    && [ -s "$conf_keep_tsv" ] \
    && [ -s "$conf_selected_merged" ] \
    && { [ "$RUN_CONCAT_CHECKV_EXPECTED" != "1" ] || { [ -s "$concat_expected_tsv" ] && [ -s "$concat_expected_summary_tsv" ]; }; } \
    && { [ "$RUN_VIRIDIC" != "1" ] || { [ -s "$viridic_bin_tsv" ] && [ -s "$viridic_summary_tsv" ]; }; } \
    && { [ "$RUN_TAXMYPHAGE" != "1" ] || { [ -s "$taxmyphage_bin_tsv" ] && [ -s "$taxmyphage_summary_tsv" ]; }; } \
    && [ "$genomad_ok" -eq 1 ]; then
    echo -e "$sample\tok\t$total_bins\t$clean_count\t$merged_contigs\t$note" >> "$MASTER_SUMMARY"
    printf "timestamp\t%s\nthreads\t%s\nnote\t%s\n" "$(date -Is)" "$THREADS" "$note" > "$sample_done_marker"
  else
    echo -e "$sample\tpartial\t$total_bins\t$clean_count\t$merged_contigs\tmissing_required_outputs" >> "$MASTER_SUMMARY"
    rm -f "$sample_done_marker"
  fi
done

taxmyphage_first=1
: > "$OVERALL_TAXMYPHAGE_BIN"
for sample in $(seq "$START_SAMPLE" "$END_SAMPLE"); do
  tbin="$ANOT_ROOT/$sample/_sample_level/taxmyphage_bin_${sample}.tsv"
  if [ -s "$tbin" ]; then
    if [ "$taxmyphage_first" -eq 1 ]; then
      cat "$tbin" > "$OVERALL_TAXMYPHAGE_BIN"
      taxmyphage_first=0
    else
      tail -n +2 "$tbin" >> "$OVERALL_TAXMYPHAGE_BIN"
    fi
  fi
done

if [ -s "$OVERALL_TAXMYPHAGE_BIN" ]; then
  python3 - "$OVERALL_TAXMYPHAGE_BIN" "$OVERALL_TAXMYPHAGE_SUMMARY" <<'PY'
import csv
import sys
from collections import Counter

in_tsv = sys.argv[1]
out_tsv = sys.argv[2]
rows = list(csv.DictReader(open(in_tsv), delimiter="\t"))
levels = Counter((r.get("taxmyphage_taxon_level") or "").strip() for r in rows if (r.get("taxmyphage_taxon_level") or "").strip())

with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["n_bins", len(rows)])
    w.writerow(["n_bins_with_taxon", sum(1 for r in rows if (r.get("taxmyphage_taxon_name") or "").strip())])
    for k, v in sorted(levels.items()):
        w.writerow([f"taxmyphage_level::{k}", v])
PY
fi

viridic_first=1
: > "$OVERALL_VIRIDIC_BIN"
for sample in $(seq "$START_SAMPLE" "$END_SAMPLE"); do
  vbin="$ANOT_ROOT/$sample/_sample_level/viridic_bin_${sample}.tsv"
  if [ -s "$vbin" ]; then
    if [ "$viridic_first" -eq 1 ]; then
      cat "$vbin" > "$OVERALL_VIRIDIC_BIN"
      viridic_first=0
    else
      tail -n +2 "$vbin" >> "$OVERALL_VIRIDIC_BIN"
    fi
  fi
done

if [ -s "$OVERALL_VIRIDIC_BIN" ]; then
  python3 - "$OVERALL_VIRIDIC_BIN" "$OVERALL_VIRIDIC_SUMMARY" <<'PY'
import csv
import statistics
import sys

bin_tsv = sys.argv[1]
out_tsv = sys.argv[2]
rows = list(csv.DictReader(open(bin_tsv), delimiter="\t"))

best = []
for r in rows:
    try:
        best.append(float((r.get("viridic_best_similarity_pct") or "").strip()))
    except Exception:
        pass

species = set((r.get("viridic_species_cluster") or "").strip() for r in rows if (r.get("viridic_species_cluster") or "").strip())
genus = set((r.get("viridic_genus_cluster") or "").strip() for r in rows if (r.get("viridic_genus_cluster") or "").strip())

with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["n_bins", len(rows)])
    w.writerow(["n_bins_with_species_cluster", sum(1 for r in rows if (r.get("viridic_species_cluster") or "").strip())])
    w.writerow(["n_bins_with_genus_cluster", sum(1 for r in rows if (r.get("viridic_genus_cluster") or "").strip())])
    w.writerow(["n_species_clusters_total", len(species)])
    w.writerow(["n_genus_clusters_total", len(genus)])
    w.writerow(["median_best_similarity_pct", f"{statistics.median(best):.3f}" if best else ""])
PY
fi

# Build overall consolidated 4-tool tables for the processed range.
bin_first=1
contig_first=1
: > "$OVERALL_BIN"
: > "$OVERALL_CONTIG"
for sample in $(seq "$START_SAMPLE" "$END_SAMPLE"); do
  sbin="$ANOT_ROOT/$sample/_sample_level/taxonomy_4tools_bin_${sample}.tsv"
  scontig="$ANOT_ROOT/$sample/_sample_level/taxonomy_4tools_contig_${sample}.tsv"
  if [ -s "$sbin" ]; then
    if [ "$bin_first" -eq 1 ]; then
      cat "$sbin" > "$OVERALL_BIN"
      bin_first=0
    else
      tail -n +2 "$sbin" >> "$OVERALL_BIN"
    fi
  fi
  if [ -s "$scontig" ]; then
    if [ "$contig_first" -eq 1 ]; then
      cat "$scontig" > "$OVERALL_CONTIG"
      contig_first=0
    else
      tail -n +2 "$scontig" >> "$OVERALL_CONTIG"
    fi
  fi
done

# Build overall vMAG confidence tables for the processed range.
conf_first=1
conf_keep_first=1
: > "$OVERALL_CONF"
: > "$OVERALL_CONF_KEEP"
for sample in $(seq "$START_SAMPLE" "$END_SAMPLE"); do
  sconf="$ANOT_ROOT/$sample/_sample_level/vmag_confidence_${sample}.tsv"
  sconf_keep="$ANOT_ROOT/$sample/_sample_level/vmag_confidence_filtered_${sample}_${VMAG_KEEP_MIN_CONF,,}plus.tsv"
  if [ -s "$sconf" ]; then
    if [ "$conf_first" -eq 1 ]; then
      cat "$sconf" > "$OVERALL_CONF"
      conf_first=0
    else
      tail -n +2 "$sconf" >> "$OVERALL_CONF"
    fi
  fi
  if [ -s "$sconf_keep" ]; then
    if [ "$conf_keep_first" -eq 1 ]; then
      cat "$sconf_keep" > "$OVERALL_CONF_KEEP"
      conf_keep_first=0
    else
      tail -n +2 "$sconf_keep" >> "$OVERALL_CONF_KEEP"
    fi
  fi
done

if [ -s "$OVERALL_BIN" ]; then
  python3 - "$OVERALL_BIN" "$OVERALL_SUMMARY" <<'PY'
import csv
import sys
from collections import Counter

bin_tsv = sys.argv[1]
out_tsv = sys.argv[2]

rows = []
with open(bin_tsv) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

n_bins = len(rows)
agree = sum(1 for r in rows if r.get("all_tools_agree") == "True")
tools_ge2 = sum(1 for r in rows if int(r.get("n_tools_called", "0") or 0) >= 2)
any_call = sum(1 for r in rows if int(r.get("n_tools_called", "0") or 0) >= 1)
rules = Counter(r.get("consensus_rule", "") for r in rows)

with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["n_bins", n_bins])
    w.writerow(["n_bins_any_call", any_call])
    w.writerow(["n_bins_tools_ge2", tools_ge2])
    w.writerow(["n_bins_all_tools_agree", agree])
    for k, v in sorted(rules.items()):
        w.writerow([f"consensus_rule::{k}", v])
PY
fi

if [ -s "$OVERALL_CONF" ]; then
  python3 - "$OVERALL_CONF" "$OVERALL_CONF_SUMMARY" <<'PY'
import csv
import sys
from collections import Counter

conf_tsv = sys.argv[1]
summary_tsv = sys.argv[2]

rows = []
with open(conf_tsv) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

conf = Counter(r.get("final_confidence", "") for r in rows)
keep = sum(1 for r in rows if r.get("keep") == "True")

with open(summary_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["n_bins", len(rows)])
    w.writerow(["n_keep", keep])
    w.writerow(["n_high", conf.get("High", 0)])
    w.writerow(["n_medium", conf.get("Medium", 0)])
    w.writerow(["n_low", conf.get("Low", 0)])
PY
fi

# Build overall geNomad tables for the processed range.
genomad_bin_first=1
genomad_contig_first=1
: > "$OVERALL_GENOMAD_BIN"
: > "$OVERALL_GENOMAD_CONTIG"
for sample in $(seq "$START_SAMPLE" "$END_SAMPLE"); do
  gbin="$ANOT_ROOT/$sample/_sample_level/genomad_bin_${sample}.tsv"
  gcontig="$ANOT_ROOT/$sample/_sample_level/genomad_contig_${sample}.tsv"
  if [ -s "$gbin" ]; then
    if [ "$genomad_bin_first" -eq 1 ]; then
      cat "$gbin" > "$OVERALL_GENOMAD_BIN"
      genomad_bin_first=0
    else
      tail -n +2 "$gbin" >> "$OVERALL_GENOMAD_BIN"
    fi
  fi
  if [ -s "$gcontig" ]; then
    if [ "$genomad_contig_first" -eq 1 ]; then
      cat "$gcontig" > "$OVERALL_GENOMAD_CONTIG"
      genomad_contig_first=0
    else
      tail -n +2 "$gcontig" >> "$OVERALL_GENOMAD_CONTIG"
    fi
  fi
done

if [ -s "$OVERALL_GENOMAD_BIN" ]; then
  python3 - "$OVERALL_GENOMAD_BIN" "$OVERALL_GENOMAD_SUMMARY" <<'PY'
import csv
import sys

bin_tsv = sys.argv[1]
out_tsv = sys.argv[2]

rows = list(csv.DictReader(open(bin_tsv), delimiter="\t"))
n_bins = len(rows)
n_tax = sum(1 for r in rows if r.get("genomad_taxonomy_bin", ""))

with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["n_bins", n_bins])
    w.writerow(["n_bins_with_taxonomy", n_tax])
PY
fi

concat_first=1
: > "$OVERALL_CONCAT_EXPECTED"
for sample in $(seq "$START_SAMPLE" "$END_SAMPLE"); do
  cexp="$ANOT_ROOT/$sample/_sample_level/vmag_concat_checkv_expected_${sample}.tsv"
  if [ -s "$cexp" ]; then
    if [ "$concat_first" -eq 1 ]; then
      cat "$cexp" > "$OVERALL_CONCAT_EXPECTED"
      concat_first=0
    else
      tail -n +2 "$cexp" >> "$OVERALL_CONCAT_EXPECTED"
    fi
  fi
done

if [ -s "$OVERALL_CONCAT_EXPECTED" ]; then
  python3 - "$OVERALL_CONCAT_EXPECTED" "$OVERALL_CONCAT_EXPECTED_SUMMARY" <<'PY'
import csv
import statistics
import sys

in_tsv = sys.argv[1]
out_tsv = sys.argv[2]

rows = list(csv.DictReader(open(in_tsv), delimiter="\t"))
n = len(rows)
n_tax = sum(1 for r in rows if (r.get("taxonomy_consensus_taxon") or "").strip())
n_eligible = sum(1 for r in rows if (r.get("ncbi_rank") or "").strip().lower() in {"species", "genus", "family"})
n_expected = sum(1 for r in rows if (r.get("expected_genome_bp_median") or "").strip())

rec = []
comp = []
for r in rows:
    try:
        rec.append(float(r.get("recovered_bp_pct") or ""))
    except Exception:
        pass
    try:
        comp.append(float(r.get("checkv_concat_completeness") or ""))
    except Exception:
        pass

with open(out_tsv, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric", "value"])
    w.writerow(["n_bins", n])
    w.writerow(["n_bins_with_taxon", n_tax])
    w.writerow(["n_bins_rank_species_genus_family", n_eligible])
    w.writerow(["n_bins_with_expected_size", n_expected])
    w.writerow(["median_recovered_bp_pct", f"{statistics.median(rec):.2f}" if rec else ""])
    w.writerow(["median_checkv_concat_completeness", f"{statistics.median(comp):.2f}" if comp else ""])
PY
fi

echo "DONE: $MASTER_SUMMARY"
echo "DONE: $OVERALL_BIN"
echo "DONE: $OVERALL_CONTIG"
echo "DONE: $OVERALL_SUMMARY"
echo "DONE: $OVERALL_CONF"
echo "DONE: $OVERALL_CONF_KEEP"
echo "DONE: $OVERALL_CONF_SUMMARY"
echo "DONE: $OVERALL_GENOMAD_BIN"
echo "DONE: $OVERALL_GENOMAD_CONTIG"
echo "DONE: $OVERALL_GENOMAD_SUMMARY"
echo "DONE: $OVERALL_CONCAT_EXPECTED"
echo "DONE: $OVERALL_CONCAT_EXPECTED_SUMMARY"
echo "DONE: $OVERALL_VIRIDIC_BIN"
echo "DONE: $OVERALL_VIRIDIC_SUMMARY"
echo "DONE: $OVERALL_TAXMYPHAGE_BIN"
echo "DONE: $OVERALL_TAXMYPHAGE_SUMMARY"
