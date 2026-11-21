# Run commands

See [README.md](README.md#run) for the supported commands and parameters.

```bash
nextflow run workflow.nf \
  --datasets /absolute/path/samples.csv \
  --ref_file /absolute/path/reference.fna \
  --anno_file /absolute/path/annotation.gff \
  --outdir results \
  --skip_rrna true \
  -resume
```

The pinned Conda environment is enabled by default. Use `-profile no_conda`
only when every required tool is already installed and available on `PATH`.
