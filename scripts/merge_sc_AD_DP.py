import pandas as pd
import numpy as np
import sys, argparse, os
from pathlib import Path
import glob
import logging
import re


def main(args):

    ##### ----- process input ----- #####
    sample_name = args.sample_name
    output_dir = Path(args.output_dir)
    # read in candidate alleles (biallelic records)
    # index column needs to be set to CHROM:POS:REF/ALT
    alleles_f = args.candidate_alleles_df
    alleles_df = pd.read_csv(alleles_f, sep="\t", header=None, index_col = 0) 
    
    input_dir = Path(args.input_dir) # this should contain all the AD and DP dfs
    # note for AD dfs: the index needs to be in the format: CHROM:POS:REF/ALT
    sc_AD_fs = list(input_dir.glob("*AD.csv"))
    logging.info(f"Found {len(sc_AD_fs)} AD files")
    # note for DP dfs: the index needs to be in the format: CHROM:POS
    sc_DP_fs = list(input_dir.glob("*DP.csv"))
    logging.info(f"Found {len(sc_DP_fs)} DP files")
    # sanity checks
    if not len(sc_AD_fs) == len(sc_DP_fs):
        raise ValueError("Number of found AD and DP files do not match!")

    if args.sc_bams_dir is not None:
        sc_bams_dir = Path(args.sc_bams_dir) # this should contain all the sc_bams. for sanity check
        num_sc_bams = len([f for f in sc_bams_dir.iterdir() if f.is_file() and f.suffix == ".bam"])
        logging.info(f"Found {num_sc_bams} sc_bams")
        if not num_sc_bams == len(sc_AD_fs):
            raise ValueError("Number of found sc-bams does not match the number of AD and DP files!")
    
    ########################################
    # ----- 0. format the alleles df ----- #
    alleles_df['pos'] = alleles_df.index.str.split(":").str[:2].str.join(":") # CHROM:POS
    alleles_df['locus'] = alleles_df.index.str.split('/').str[0] # CHROM:POS:REF
    alleles_df['alleles'] = alleles_df.index # CHROM:POS:REF/ALT
    alleles_df['ref'] = alleles_df.index.str.split(':').str[2].str.split('/').str[0] # REF
    alleles_df['alt'] = alleles_df.index.str.split(':').str[2].str.split('/').str[1] # ALT

    # ----- 1. format AD and DP dfs ----- #
    # get the AD layer
    sc_AD_cell_ids = [re.findall(r'cell_\d+', f.name)[0] for f in sc_AD_fs] # infer cell_id from file names
    sc_AD_dfs_list = [pd.read_csv(f, index_col=0, names = ['AD']) for f in sc_AD_fs]
    sc_AD_dfs_list = [
        df_i[~df_i.index.duplicated(keep='first')] 
        for df_i in sc_AD_dfs_list
        ] # remove duplicated SNVs (CHR:POS:REF/ALT)
    sc_AD_dfs_map = dict(zip(
        sc_AD_cell_ids,
        sc_AD_dfs_list
    ))
    merged_sc_AD_df = pd.concat(sc_AD_dfs_map, axis=1)
    merged_sc_AD_df.columns = sc_AD_dfs_map.keys()
    ordered_sc_idx = merged_sc_AD_df.columns.sort_values() # order the single cell indices !!!!
    ordered_var = merged_sc_AD_df.index.sort_values() # order the variants !!!!
    merged_sc_AD_df = merged_sc_AD_df.loc[ordered_var, ordered_sc_idx]

    # get the DP layer
    sc_DP_cell_ids = [re.findall(r'cell_\d+', f.name)[0] for f in sc_DP_fs]
    sc_DP_dfs_list = [pd.read_csv(f, index_col=0, names = ['DP']) for f in sc_DP_fs] # index_col is now CHROM:POS:REF

    # sc_DP_dfs_list = [
    #     df_i[~df_i.index.duplicated(keep='first')] for df_i in sc_DP_dfs_list
    #     ] # drop duplicated loci (CHR:POS:REF)
    # # ^^^^^^^^ this is problematic for indels with len(REF) != 1
    
    # set index from CHROM:POS:REF to CHROM:POS
    sc_DP_dfs_list = [
        pd.DataFrame(data = df_i.values, index = df_i.index.str.split(':').str[:2].str.join(':'), columns = ['DP'])
        for df_i in sc_DP_dfs_list
    ]
    # aggregate duplicated positions with the `max` function
    sc_DP_dfs_list = [
        df_i.fillna(0).groupby(level=0).agg('max')
        for df_i in sc_DP_dfs_list
        ]

    sc_DP_dfs_map = dict(zip(
        sc_DP_cell_ids,
        sc_DP_dfs_list
    ))
    merged_sc_DP_df = pd.concat(sc_DP_dfs_map, axis=1)
    merged_sc_DP_df.columns = sc_DP_dfs_map.keys()
    ordered_sc_idx = merged_sc_DP_df.columns.sort_values() # order the single cell indices !!!!
    ordered_var = merged_sc_DP_df.index.sort_values() # order the variants !!!!
    merged_sc_DP_df = merged_sc_DP_df.loc[ordered_var, ordered_sc_idx]

    # ----- 2. populate the output dfs ----- #
    filled_DP_df = alleles_df.apply(
        lambda row: merged_sc_DP_df.loc[row['pos'], :] 
        if row['pos'] in merged_sc_DP_df.index 
        else pd.Series([0]*len(sc_DP_dfs_map), index = merged_sc_DP_df.columns),
        axis=1
        )
    # # for indels where the ref allele is not one single locus, fill in their DP in cells where the mutant allele is not present
    # filled_DP_df.loc[(alleles_df['ref'].str.len() != 1), :].loc[
    #     np.isnan(filled_DP_df)
    # ]
    filled_AD_df = alleles_df.apply(
        lambda row: merged_sc_AD_df.loc[row['alleles'], :] 
        if row['alleles'] in merged_sc_AD_df.index 
        else pd.Series([0]*len(sc_AD_dfs_map), index = merged_sc_AD_df.columns),
        axis=1
        )
    filled_DP_df.fillna(0, inplace=True)
    filled_AD_df.fillna(0, inplace=True)

    # ----- 3. write the output dfs ----- #
    filled_DP_df.astype(int).to_csv(output_dir / f"{sample_name}.mpileup.DP.merged.csv")
    filled_AD_df.astype(int).to_csv(output_dir / f"{sample_name}.mpileup.AD.merged.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type=str, required=True, help='sample name')
    parser.add_argument('--candidate_alleles_df', type=str, help='Dataframe containing candidate alleles in the format: CHR:POS:REF/ALT', required=True)
    parser.add_argument('--input_dir', type=str, help='Directory containing AD and DP files', required=True)
    parser.add_argument('--sc_bams_dir', type=str, help='Directory containing single cell bams; for sanity check', default=None)
    parser.add_argument('--output_dir', type=str, help='Directory to write output files', required=True)


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)

