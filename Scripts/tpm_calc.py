import pandas as pd
import numpy as np


def gene_specific_fragment_length(gene_length, frag_len_probs):
    total_prob = sum(prob for l, prob in frag_len_probs.items() if l <= gene_length)
    if total_prob == 0:
        return 0
    expected_len = sum(l * prob for l, prob in frag_len_probs.items() if l <= gene_length)
    return expected_len / total_prob


def tpm_calc(sample_id, condition, metric_file, featurecount_file, output_file):
    df_fragment_count = pd.read_table(metric_file, comment='#')

    df_fragment_count['Probability'] = df_fragment_count['Count'] / df_fragment_count['Count'].sum()

    frag_len_probs = dict(zip(df_fragment_count['Fragment_Length'], df_fragment_count['Probability']))

    df_feature = pd.read_table(featurecount_file, comment='#')

    df_feature['Fragment_Length'] = df_feature['Length'].apply(lambda l: gene_specific_fragment_length(l, frag_len_probs))

    df_feature[sample_id+'_Effective_Length'] = df_feature['Length'] - df_feature['Fragment_Length'] + 1

    df_feature.set_index('Geneid', inplace=True)

    df_feature[sample_id+'_'+condition+'_TPM'] = (df_feature[sample_id+'.sorted.bam'] / df_feature[sample_id+'_Effective_Length']) / (df_feature[sample_id+'.sorted.bam'] / df_feature[sample_id+'_Effective_Length']).sum() * 1e6

    df_feature.drop(['Fragment_Length'], axis=1, inplace=True)

    df_feature.to_csv(output_file)


# Usage
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 6:
        print('Usage: python tpm_calc.py <sample_id> <condition> <metric_file> <featurecount_file> <output_file>')
        sys.exit(1)

    sample_id = sys.argv[1]
    condition = sys.argv[2]
    metric_file = sys.argv[3]
    featurecount_file = sys.argv[4]
    output_file = sys.argv[5]
    tpm_calc(sample_id, condition, metric_file, featurecount_file, output_file)
