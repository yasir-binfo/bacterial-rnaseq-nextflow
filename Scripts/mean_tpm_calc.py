#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


BASE_COLUMNS = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
NORMALIZED_SUFFIX = "_normalized.results"


def mean_tpm_calc(normalized_files, samplesheet, output_file):
    metadata = pd.read_csv(samplesheet, dtype=str).fillna("")
    required = {"sample_id", "condition"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Sample sheet is missing columns: {', '.join(sorted(missing))}")
    if metadata["sample_id"].duplicated().any():
        duplicates = metadata.loc[metadata["sample_id"].duplicated(), "sample_id"]
        raise ValueError(f"Duplicate sample_id values: {', '.join(duplicates)}")

    sample_to_condition = dict(zip(metadata["sample_id"], metadata["condition"]))
    merged = None
    condition_columns = {}

    for normalized_file in normalized_files:
        filename = Path(normalized_file).name
        if not filename.endswith(NORMALIZED_SUFFIX):
            raise ValueError(f"Unexpected normalized filename: {filename}")

        sample_id = filename[: -len(NORMALIZED_SUFFIX)]
        if sample_id not in sample_to_condition:
            raise ValueError(f"Sample '{sample_id}' is absent from the sample sheet")

        condition = sample_to_condition[sample_id]
        tpm_column = f"{sample_id}_{condition}_TPM"
        frame = pd.read_csv(normalized_file)

        missing_base = set(BASE_COLUMNS).difference(frame.columns)
        if missing_base:
            raise ValueError(
                f"{filename} is missing columns: {', '.join(sorted(missing_base))}"
            )
        if tpm_column not in frame.columns:
            raise ValueError(f"{filename} is missing TPM column '{tpm_column}'")

        sample_frame = frame[BASE_COLUMNS + [tpm_column]].copy()
        if merged is None:
            merged = sample_frame
        else:
            merged = merged.merge(
                sample_frame[["Geneid", tpm_column]],
                on="Geneid",
                how="inner",
                validate="one_to_one",
            )
        condition_columns.setdefault(condition, []).append(tpm_column)

    if merged is None:
        raise ValueError("No normalized result files were supplied")

    result = merged[BASE_COLUMNS].copy()
    for condition, columns in condition_columns.items():
        result[condition] = merged[columns].mean(axis=1)

    result.to_csv(output_file, index=False)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate condition-level mean TPM from normalized samples"
    )
    parser.add_argument("normalized_files", nargs="+")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    mean_tpm_calc(args.normalized_files, args.samplesheet, args.output)
