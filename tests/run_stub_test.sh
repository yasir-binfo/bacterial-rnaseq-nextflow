#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_dir="$(mktemp -d)"
trap 'rm -rf "${test_dir}"' EXIT

cd "${repo_dir}"

nextflow run workflow.nf \
  -profile test \
  -stub-run \
  -ansi-log false \
  -work-dir "${test_dir}/work" \
  --datasets tests/data/samples.csv \
  --ref_file tests/data/reference.fna \
  --anno_file tests/data/annotation.gff \
  --outdir "${test_dir}/results" \
  --skip_rrna true

test -e "${test_dir}/results/quantification/Expression.csv"
test -e "${test_dir}/results/bigwig/control_group.bw"
test -e "${test_dir}/results/bigwig/treated_group.bw"

nextflow run workflow.nf \
  -profile test \
  -stub-run \
  -ansi-log false \
  -work-dir "${test_dir}/rrna-work" \
  --datasets tests/data/samples.csv \
  --ref_file tests/data/reference.fna \
  --anno_file tests/data/annotation.gff \
  --outdir "${test_dir}/rrna-results" \
  --skip_rrna false \
  --rrna_16s tests/data/rrna_1.fasta \
  --rrna_23s tests/data/rrna_2.fasta \
  --rrna_5s tests/data/rrna_3.fasta

test -e "${test_dir}/rrna-results/clean_reads/pe_sample_1_non_rRNA_reads_fwd.fastq"
test -e "${test_dir}/rrna-results/clean_reads/se_sample_1_non_rRNA_reads.fastq"
