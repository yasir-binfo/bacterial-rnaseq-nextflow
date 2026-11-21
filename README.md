# Bacterial RNA-seq pipeline

A reproducible Nextflow DSL2 workflow for bacterial RNA-seq trimming, optional
rRNA depletion, alignment, strand inference, gene-level quantification, TPM
normalization, and condition-level coverage tracks. It supports paired-end,
single-end, and mixed-layout sample sheets.

## Workflow

```text
sample sheet
  -> fastp trimming
  -> optional SortMeRNA filtering
  -> Bowtie 2 alignment -> sorted BAM
  -> RSeQC strand inference
  -> featureCounts
  -> effective-length TPM -> condition mean TPM
  -> condition-level CPM BigWig
```

The Bowtie 2 index, annotation, reads, rRNA databases, and repository helper
scripts are all staged as Nextflow inputs. No process relies on `../../../`, a
published output from another process, or the directory from which a task runs.

## Requirements

- Nextflow 25.10.2 or a compatible release
- Java compatible with that Nextflow release
- Conda

All analysis tools and Python libraries are pinned in `environment.yml`.
Nextflow creates and uses that environment automatically. The first run can take
several minutes while Conda resolves and downloads packages.

For an existing environment that already provides every dependency, use
`-profile no_conda`.

## Inputs

### Sample sheet

The CSV must have these exact headers:

```csv
sample_id,condition,fastq_1,fastq_2
control_1,control_group,/data/control_1_R1.fastq.gz,/data/control_1_R2.fastq.gz
treated_1,treated_group,/data/treated_1.fastq.gz,
```

- `sample_id`, `condition`, and `fastq_1` are required.
- Leave `fastq_2` empty for single-end samples.
- IDs and conditions may contain letters, numbers, dots, underscores, and
  hyphens. Other characters are rejected because condition names become output
  filenames.
- Sample IDs must be unique.
- Absolute FASTQ paths are recommended. Relative paths are resolved from the
  directory in which Nextflow is launched.

### Reference and annotation

- `--ref_file`: reference genome FASTA.
- `--anno_file`: GFF3 annotation. Quantification selects `gene` and
  `pseudogene` features and groups them by `locus_tag`. Strand inference uses
  `CDS` features and their `Parent` attribute.

Confirm those feature types and attributes exist in the supplied annotation.

## Run

From the repository root:

```bash
nextflow run workflow.nf \
  --datasets /absolute/path/samples.csv \
  --ref_file /absolute/path/reference.fna \
  --anno_file /absolute/path/annotation.gff \
  --outdir results \
  --skip_rrna true \
  -resume
```

The helper scripts are loaded from `Scripts/` automatically; script-path
parameters are neither required nor supported.

### Enable rRNA filtering

`--skip_rrna` defaults to `true`. For bacteria-only samples, the required
references are bacterial 16S, bacterial 23S, and 5S rRNA. Download the pinned
legacy separated references distributed with SortMeRNA v4.2.0:

```bash
mkdir -p databases/sortmerna

curl -L \
  https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-bac-16s-id90.fasta \
  -o databases/sortmerna/silva-bac-16s-id90.fasta

curl -L \
  https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-bac-23s-id98.fasta \
  -o databases/sortmerna/silva-bac-23s-id98.fasta

curl -L \
  https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/rfam-5s-database-id98.fasta \
  -o databases/sortmerna/rfam-5s-database-id98.fasta
```

Direct downloads:

- [SILVA bacterial 16S, 90% identity](https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-bac-16s-id90.fasta)
- [SILVA bacterial 23S, 98% identity](https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/silva-bac-23s-id98.fasta)
- [Rfam 5S, 98% identity](https://raw.githubusercontent.com/biocore/sortmerna/v4.2.0/data/rRNA_databases/rfam-5s-database-id98.fasta)

These three references are suitable for pure bacterial isolates and
bacteria-only cultures. They are not comprehensive for archaeal, eukaryotic,
host-associated, or mixed environmental samples. Bacteria do not have 5.8S
rRNA, so the Rfam 5.8S database is not needed for a bacteria-only workflow.

Current SortMeRNA releases recommend their newer combined default database,
which includes bacterial, archaeal, and eukaryotic references. See the
[official SortMeRNA database section](https://github.com/sortmerna/sortmerna#databases)
if broad contamination screening is desired.

To run filtering with the three bacterial references:

```bash
nextflow run workflow.nf \
  --datasets /absolute/path/samples.csv \
  --ref_file /absolute/path/reference.fna \
  --anno_file /absolute/path/annotation.gff \
  --outdir results \
  --skip_rrna false \
  --rrna_16s databases/sortmerna/silva-bac-16s-id90.fasta \
  --rrna_23s databases/sortmerna/silva-bac-23s-id98.fasta \
  --rrna_5s databases/sortmerna/rfam-5s-database-id98.fasta \
  -resume
```

Each SortMeRNA task uses a private work directory inside its Nextflow task. It
does not delete or share `~/sortmerna` state.

## Parameters

| Parameter | Default | Purpose |
|---|---:|---|
| `--datasets` | `datasets.csv` | Input sample sheet |
| `--ref_file` | required | Reference genome FASTA |
| `--anno_file` | required | GFF3 annotation |
| `--outdir` | `results` | Published results directory |
| `--q_score` | `5` | fastp mean-quality cutoff |
| `--skip_rrna` | `true` | Skip SortMeRNA when true |
| `--rrna_16s` | none | Bacterial 16S rRNA FASTA required when filtering |
| `--rrna_23s` | none | Bacterial 23S rRNA FASTA required when filtering |
| `--rrna_5s` | none | 5S rRNA FASTA required when filtering |

## Outputs

All published files are written beneath `--outdir`:

| Directory | Contents |
|---|---|
| `trimming/` | Trimmed FASTQs and fastp reports |
| `clean_reads/` | Non-rRNA FASTQs when filtering is enabled |
| `index/` | Bowtie 2 index bundle |
| `mapping/` | Coordinate-sorted BAM files |
| `strandedness/` | RSeQC strand reports |
| `quantification/` | featureCounts results, length metrics, per-sample TPM, and `Expression.csv` |
| `bigwig/` | Condition-level CPM-normalized BigWigs |
| `pipeline_info/` | Nextflow report, timeline, trace, and DAG |

Processes now use Nextflow's default fail-fast behavior. Quantification or
BigWig failures therefore fail the workflow instead of silently dropping output.

## Bowtie 2 index channel

`BUILD_INDEX` receives the reference as a value channel and emits the complete
`index*.bt2*` set as a reusable value-channel bundle:

```groovy
index_bundle = BUILD_INDEX(reference).index_files
mapped_pe = MAP_PE(clean_pe, index_bundle).bam
mapped_se = MAP_SE(clean_se, index_bundle).bam
```

Because `reference` is a value channel, the process output is also a value
channel and can be reused by every mapping task. Each task receives all six
small-index files (or all six large-index files) in its work directory and calls
`bowtie2 -x index`. This is different from passing one queue emission per index
file, which causes incomplete indexes, or pairing a one-item queue with many
reads, which can schedule only one sample. Mapping cannot begin until
`BUILD_INDEX` completes.

## Reproducibility and resources

`nextflow.config` defines fail-fast execution and three resource classes:

| Label | CPUs | Memory | Time |
|---|---:|---:|---:|
| `small` | 2 | 4 GB | 2 hours |
| `medium` | 4 | 8 GB | 8 hours |
| `large` | 8 | 16 GB | 12 hours |

Adjust these values for the executor and dataset size. Tool versions are pinned
in `environment.yml`; update and retest that file deliberately when upgrading.

## Tests

The repository includes mixed PE/SE test metadata, tiny reference/annotation
fixtures, a Nextflow stub-run smoke test, and a Python unit test for condition
names containing underscores.

Run them from the repository root:

```bash
python -m unittest discover -s tests -v
bash tests/run_stub_test.sh
```

The stub test validates the complete dataflow and declared outputs without
running bioinformatics tools or creating the Conda environment. It is not a
scientific validation dataset; a real small-dataset integration test remains a
recommended future addition.

## Helper scripts

- `Scripts/gff_to_bed.py`: GFF3 to BED12 conversion for RSeQC.
- `Scripts/gff_to_saf.py`: optional GFF3 to SAF converter; not used by the
  current workflow.
- `Scripts/infer_strandedness.py`: maps RSeQC fractions to featureCounts strand
  codes.
- `Scripts/tpm_calc.py`: calculates effective-length-adjusted TPM.
- `Scripts/mean_tpm_calc.py`: uses the sample sheet to calculate condition means,
  avoiding ambiguous underscore-based parsing.

## Remaining limitations

- The committed test is a stub-run structural test, not an end-to-end scientific
  validation against expected biological results.
- Resource values are reasonable defaults, not guarantees for every genome or
  sample depth.
- Conda package availability depends on the configured channels and platform.

## License and citation

No license or citation metadata is currently included. Add a `LICENSE` file and
the appropriate workflow/data citations before distributing the repository.
