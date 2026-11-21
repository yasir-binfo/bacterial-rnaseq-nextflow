# Bacterial RNA-seq pipeline

A Nextflow workflow for bacterial RNA-seq preprocessing, alignment, gene-level
quantification, TPM normalization, and condition-level coverage tracks. The
workflow accepts paired-end (PE), single-end (SE), or mixed PE/SE sample sheets.

> **Status:** research workflow under review. Read [Known limitations](#known-limitations)
> before using the results in production analyses.

## Workflow overview

```text
sample sheet
    |
    +-- PE/SE trimming (fastp)
    |
    +-- optional rRNA filtering (SortMeRNA)
    |
    +-- alignment to reference (Bowtie 2)
    |
    +-- coordinate-sorted BAM (SAMtools)
    |
    +-- strandedness inference (RSeQC)
    |
    +-- gene/pseudogene counts (featureCounts)
    |
    +-- per-sample TPM + condition mean TPM
    |
    `-- condition-level CPM-normalized BigWig (deepTools)
```

## Requirements

- Nextflow (the workflow was reviewed with 25.10.2)
- Java compatible with the installed Nextflow release
- Python 3 with `pandas` and `numpy`
- fastp
- SortMeRNA (only when rRNA filtering is enabled)
- Bowtie 2
- SAMtools
- RSeQC (`infer_experiment.py`)
- Subread (`featureCounts`)
- deepTools (`bamCoverage`)

The current repository does not include a Conda environment, container profile,
or `nextflow.config`; all programs must therefore be available on `PATH`.

## Input files

### Sample sheet

Supply a CSV file with these exact headers:

```csv
sample_id,condition,fastq_1,fastq_2
control_1,control,/absolute/path/control_1_R1.fastq.gz,/absolute/path/control_1_R2.fastq.gz
treated_1,treated,/absolute/path/treated_1.fastq.gz,
```

- `sample_id` should be unique and safe for use in filenames (letters, numbers,
  underscores, and hyphens are recommended).
- `condition` labels are used to group BAM files and calculate mean TPM values.
- `fastq_1` is required.
- `fastq_2` is required for paired-end samples and must be empty for single-end
  samples.
- Absolute paths are recommended because Nextflow processes run in isolated work
  directories.

### Reference and annotation

- `--ref_file`: reference genome in FASTA format.
- `--anno_file`: annotation in GFF3 format. Quantification currently selects
  `gene` and `pseudogene` features and groups them by the `locus_tag` attribute.
  Confirm that those feature types and attributes exist in your annotation.

### SortMeRNA databases

When `--skip_rrna false` is used, the workflow expects six rRNA reference FASTA
files. Their locations can be supplied individually with `--database1` through
`--database6`; see [Parameters](#parameters).

## Running the pipeline

From the repository root:

```bash
nextflow run workflow.nf \
  --datasets /absolute/path/samples.csv \
  --ref_file /absolute/path/reference.fna \
  --anno_file /absolute/path/annotation.gff \
  --gff_to_bed_script "$(pwd)/Scripts/gff_to_bed.py" \
  --gff_to_saf_script "$(pwd)/Scripts/gff_to_saf.py" \
  --strandedness_script "$(pwd)/Scripts/infer_strandedness.py" \
  --tpm_normalization_script "$(pwd)/Scripts/tpm_calc.py" \
  --mean_tpm_calc "$(pwd)/Scripts/mean_tpm_calc.py" \
  --skip_rrna true \
  -resume
```

`--skip_rrna` defaults to `true`. To enable SortMeRNA filtering:

```bash
nextflow run workflow.nf \
  --datasets /absolute/path/samples.csv \
  --ref_file /absolute/path/reference.fna \
  --anno_file /absolute/path/annotation.gff \
  --gff_to_bed_script "$(pwd)/Scripts/gff_to_bed.py" \
  --strandedness_script "$(pwd)/Scripts/infer_strandedness.py" \
  --tpm_normalization_script "$(pwd)/Scripts/tpm_calc.py" \
  --mean_tpm_calc "$(pwd)/Scripts/mean_tpm_calc.py" \
  --skip_rrna false \
  --database1 /path/to/rfam-5.8s-database-id98.fasta \
  --database2 /path/to/rfam-5s-database-id98.fasta \
  --database3 /path/to/silva-bac-16s-id90.fasta \
  --database4 /path/to/silva-bac-23s-id98.fasta \
  --database5 /path/to/silva-arc-16s-id95.fasta \
  --database6 /path/to/silva-arc-23s-id98.fasta \
  -resume
```

## Parameters

| Parameter | Default | Purpose |
|---|---:|---|
| `--datasets` | `datasets.csv` | Input sample sheet |
| `--ref_file` | `ref` | Reference genome FASTA |
| `--anno_file` | `gff` | GFF3 annotation |
| `--q_score` | `5` | Mean-quality cutoff passed to fastp |
| `--skip_rrna` | `true` | Skip SortMeRNA when true |
| `--database1` … `--database6` | paths under `Databases/` | SortMeRNA reference FASTA files |
| `--gff_to_bed_script` | `gff_to_bed_script.py` | GFF3-to-BED helper for RSeQC |
| `--gff_to_saf_script` | `gff_to_saf_script.py` | GFF3-to-SAF helper (currently unused) |
| `--strandedness_script` | `scripts/infer_strandedness.py` | Converts RSeQC output to featureCounts strand code |
| `--tpm_normalization_script` | `tpm_calc.py` | Calculates per-sample TPM |
| `--mean_tpm_calc` | `mean_tpm_calc.py` | Combines samples and calculates condition means |

The workflow also defines internal `--anno_file_relative_path` and
`--index_path` parameters. Their defaults are coupled to Nextflow work-directory
depth; overriding them may be necessary if the directory layout changes.

## Outputs

The current workflow publishes files under `Omics/Expression/`:

| Directory | Contents |
|---|---|
| `Trimming/` | Trimmed FASTQ files and PE fastp HTML reports |
| `Clean_Reads/` | Non-rRNA reads when filtering is enabled |
| `Index/` | Bowtie 2 reference index |
| `Mapping/` | Coordinate-sorted BAM files |
| `Strandedness/` | RSeQC strand reports and BAM copies |
| `Quantification/` | featureCounts tables, length metrics, normalized tables, `Expression.csv`, and condition BigWigs |

Nextflow also creates `work/`, `.nextflow/`, `.nextflow.log`, and execution
metadata. These are excluded from Git.

## Helper scripts

- `Scripts/gff_to_bed.py`: converts selected GFF3 features (CDS by default) to
  BED12 for RSeQC.
- `Scripts/gff_to_saf.py`: converts GFF3 to SAF; it is not called by the current
  workflow.
- `Scripts/infer_strandedness.py`: maps RSeQC fractions to featureCounts codes
  (`0` unstranded, `1` forward, `2` reverse).
- `Scripts/tpm_calc.py`: calculates effective-length-adjusted TPM.
- `Scripts/mean_tpm_calc.py`: merges normalized samples and calculates mean TPM
  for each condition.

Each converter supports `--help`, for example:

```bash
python3 Scripts/gff_to_bed.py --help
```

## Known limitations

The review identified the following points that should be resolved before the
workflow is treated as reproducible or production-ready:

1. The Bowtie 2 index process is not connected to the mapping processes through
   a Nextflow channel. Mapping instead reads a fixed published path, so a fresh
   run may race the index build or reuse a stale index.
2. Annotation and index defaults use hard-coded `../../../` paths relative to
   process work directories. This is fragile and bypasses normal Nextflow file
   staging.
3. Output locations are hard-coded; the `--outdir` shown in the older
   `run_commands.md` is not used by `workflow.nf`.
4. SortMeRNA processes run `rm -rf ~/sortmerna/` before every task. Concurrent
   samples can interfere with one another and the command removes state outside
   the Nextflow task directory.
5. `errorStrategy 'ignore'` on quantification and BigWig creation can allow a run
   to finish despite missing scientific outputs.
6. Tool versions and resources (`cpus`, memory, time) are not pinned or defined.
7. Condition extraction in `mean_tpm_calc.py` splits column names on `_`; sample
   IDs or condition names containing underscores can therefore be misclassified.
8. The workflow has no automated test data or end-to-end validation test.

## Reproducibility recommendations

Before publishing releases, add a small test dataset, pin tool versions through
Conda or containers, connect all reference inputs through channels, expose a
single output directory, and fail explicitly when required outputs cannot be
created.

## License and citation

No license or citation metadata is currently included. Add a `LICENSE` file and
the appropriate workflow/data citations before distributing the repository.
