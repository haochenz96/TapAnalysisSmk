import argparse, sys
import pandas as pd
import re

def main(args):
    sample_name = args.sample_name
    AD_TXT = args.AD_TXT
    AD_TXT_FORMATTED_FOR_PY = args.AD_TXT_FORMATTED_FOR_PY

    AD_df = pd.read_csv(AD_TXT, sep='\t', header=0)

    # get index from the first 4 columns
    if not str(AD_df['# [1]CHROM'][0]).startswith('chr'):
        AD_df['# [1]CHROM'] = 'chr' + AD_df['# [1]CHROM'].astype(str)

    AD_df.index = AD_df['# [1]CHROM'].astype(str) + ':' + AD_df['[2]POS'].astype(str) + ':' + AD_df['[3]REF'].astype(str) + '/' + AD_df['[4]ALT'].astype(str)
    AD_df.drop(['# [1]CHROM', '[2]POS', '[3]REF', '[4]ALT', ], axis=1,inplace=True)
    
    # drop empty columns
    AD_df= AD_df.loc[:, ~AD_df.isna().all(axis=0)] 
    
    # format the column names
    AD_df.columns = AD_df.columns.map(lambda x: re.sub(r'\[\d+\]', '', x)).map(lambda x: re.sub(':AD','', x))
    
    if sample_name is not None:
        AD_df.columns = AD_df.columns.map(lambda x: sample_name + ':' + x)

    # fill in DP=0
    AD_df.replace('.', '0,0', inplace=True)

    # write to output
    AD_df.T.to_csv(AD_TXT_FORMATTED_FOR_PY, sep='\t', header=True, index=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type = str, default=None, help='sample name')
    parser.add_argument('--AD_TXT', type=str, help='AD matrix in TXT format, directly output by bcftools')
    parser.add_argument('--AD_TXT_FORMATTED_FOR_PY', type=str, help='formatted AD matrix in TXT format, for Python')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)