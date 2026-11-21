# Run commands

The supported command, sample-sheet schema, parameters, outputs, and current
review findings are documented in [README.md](README.md#running-the-pipeline).

From this directory, print the installed Nextflow version with:

```bash
nextflow -version
```

Use `-resume` when restarting an interrupted or previously cached run. The
workflow currently publishes to `Omics/Expression/`; `--outdir` is not
implemented.
