import pandas as pd
from glob import glob
from functools import reduce

def mean_tpm_calc(datasets, outfile):

    # Load all datasets
    dfs = [pd.read_csv(file) for file in datasets]

    # Base columns to include only once
    base_cols = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']

    # Extract base dataframe from the first dataset
    df_base = dfs[0][base_cols].copy()

    # Now concatenate TPM and effective length (or any remaining columns) from each dataset
    # You can rename columns to reflect their dataset origin
    additional_dfs = []
    for i, df in enumerate(dfs):
        # Drop base columns
        extra = df.drop(columns=base_cols + list(df.columns[df.columns.str.contains('.bam')]) + list(df.columns[df.columns.str.contains('Effective')]))

        # Keep geneid for merging
        extra.insert(0, 'Geneid', df['Geneid'])

        additional_dfs.append(extra)

    # Merge all additional dataframes on geneid
    df_extra = reduce(lambda left, right: pd.merge(left, right, on='Geneid'), additional_dfs)

    # Final combined dataframe
    df_final = pd.merge(df_base[['Geneid', 'Chr','Start', 'End', 'Strand', 'Length']], df_extra, on='Geneid')

    conditions = [x.split('_')[1] for x in df_final.columns.tolist()[6:]]
    samples_cols =  [x for x in df_final.columns.tolist()[6:]]

    condition_and_sample = {}
    for x in range(len(conditions)):
        if conditions[x] not in condition_and_sample.keys():
            condition_and_sample[conditions[x]] = [samples_cols[x]]
        else:
            condition_and_sample[conditions[x]].append(samples_cols[x])

    df_mean_tpm_cond = pd.DataFrame({'Geneid':df_final['Geneid'], 'Chr':df_final['Chr'], 'Start':df_final['Start'], 'End':df_final['End'], 'Strand':df_final['Strand'], 'Length':df_final['Length']})
    for cond in condition_and_sample:
        df_mean_tpm_cond[cond] = df_final[condition_and_sample[cond]].mean(axis=1)

    df_mean_tpm_cond.to_csv(outfile, index=False)


# Usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python mean_tpm_calc.py <file1> <file2> ... <outfile>")
        sys.exit(1)

    *datasets, outfile = sys.argv[1:]
    mean_tpm_calc(datasets, outfile)
