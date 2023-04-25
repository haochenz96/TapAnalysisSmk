# @HZ 06/28/2022
# annotate a variant list with matched bulk info
# integrated into main Snakemake workflow

# @HZ 12/18/2022 --- V2
# mandatory output to {SAMPLE_NAME}.annotated_SNV.csv
# read base-quality data from command line

import argparse, sys
import pandas as pd
import re
import numpy as np
from pathlib import Path
import yaml
import mosaic.io as mio
from tea.parse import read_vcf_to_df
from tea.format import *
from tea.utils import get_simple_timestamp
from IPython import embed
import json
from typing import Tuple
import logging
from tea.format import CONDENSED_SNV_FORMAT, check_matrix_format


def main(args):
    # ----- io -----
    SAMPLE_NAME = args.sample_name
    # --- inputs
    H5 = args.input_h5
    BULK_INFO_YAML= args.bulk_info_yaml

    OUTPUT_DIR = args.output_dir
    if OUTPUT_DIR == None:
        # assume we this is in the Snakemake workflow
        OUTPUT_DIR = Path(f'{SAMPLE_NAME}/OUTPUTS_from_m2_f')
        OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    else:
        OUTPUT_DIR = Path(OUTPUT_DIR)
        OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

    # --- outputs
    BULK_ANNOTATED_H5 = args.output_bulk_annotated_h5
    if BULK_ANNOTATED_H5 is None:
        BULK_ANNOTATED_H5 = OUTPUT_DIR / f'{SAMPLE_NAME}_DNA_CNV_m2_f.bulk_annotated.h5'
    PONV1_FILTERED_H5 = args.output_ponv1_filtered_h5
    if PONV1_FILTERED_H5 is None:
        PONV1_FILTERED_H5 = OUTPUT_DIR / f'{SAMPLE_NAME}_DNA_CNV_m2_f.ponv1_filtered.h5'

    simple_timestamp = get_simple_timestamp()

    ######### -- log -- ###########
    LOG_FILE = args.log_file
    if LOG_FILE is None:
        logging.basicConfig(
            stream=sys.stdout,
            level=logging.INFO, 
            filemode='a',
            force=True,
            format='%(asctime)s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            )
    else:
        logging.basicConfig(
            filename = LOG_FILE,
            level = logging.INFO,
            filemode='a',
            force=True,
            format='%(asctime)s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            )
    ###############################

    logging.info(f'--- {simple_timestamp} ---')
    # (0) ----- preprocess single-cell SNV matrix
    sample_obj = mio.load(H5, SAMPLE_NAME, raw = False)
    # default_genotyping
    sample_obj.dna.genotype_variants(
        min_dp = 8,
        min_alt_read = 3
    )
    VAR_DF = pd.DataFrame(index=sample_obj.dna.ids())
    # embed()
    # add sc_mut_prev info
    VAR_DF['Tapestri_result-sc_mut_prev'] = sample_obj.dna.get_attribute(
        'mut_filtered', constraint='row'
    ).sum(axis=0).values
    # embed()

    # @HZ 03/27/2023: deprecated
    # try:
    #     TLOD_DF = sample_obj.dna.get_attribute('TLOD', constraint='row+col').T.astype(float)
    #     VAR_DF['Tapestri_result-sc_max_TLOD'] = TLOD_DF.max(axis=1, skipna=True)
    #     VAR_DF['Tapestri_result-sc_mean_TLOD'] = TLOD_DF.mean(axis=1, skipna=True)
    # except:
    #     logging.warning(f'no TLOD info in {H5}')
    #     VAR_DF['Tapestri_result-sc_max_TLOD'] = np.nan
    #     VAR_DF['Tapestri_result-sc_mean_TLOD'] = np.nan

    # (1) ----- read in bulk info -----
    with open(BULK_INFO_YAML, "r") as f:
        bulk_info = yaml.safe_load(f)

    case_name = bulk_info['case_name']
    # sample_name = bulk_info["sample_name"]

    bulk_dfs = {}
    # -- a. matched bulk tumor; Broad WES PoN
    for ds in ['matched_bulk_tumor', 'broad_wes_pon']:
        if bulk_info[ds] is None:
            logging.warning(f"[WARNING] --- {ds} is not given. Skipping.")
            bulk_dfs[ds] = None
            continue
        if not isinstance(bulk_info[ds], str):
            logging.warning(f"[ERROR] --- {ds} is not given as a single path!")
        if bulk_info[ds].endswith(('txt', 'maf')):
            bulk_dfs[ds] = pd.read_csv(bulk_info[ds], sep="\t")
        elif bulk_info[ds].endswith(('vcf', 'vcf.gz')):
            bulk_dfs[ds] = read_vcf_to_df(bulk_info[ds])
        else:
            logging.warning(f'[ERROR] --- {bulk_info[ds]} is not given as a TXT or MAF file!')

        # get condensed variant format
        if bulk_dfs[ds].index.name not in ('Tapestri_format', 'condensed_format'):
            if ('Tapestri_format' or 'condensed_format') not in bulk_dfs[ds].columns:
                # embed()
                raise ValueError(f'[ERROR] --- {ds} does not have Tapestri_format or condensed_format column!')
                # try:
                #     bulk_dfs[ds]['Tapestri_format'] = bulk_dfs[ds]['CHROM'] + ":" + bulk_dfs[ds]['POS'] + ":" + bulk_dfs[ds]['REF'] + "/" + bulk_dfs[ds]['ALT']
                    
                # except KeyError:
                #     embed()
                #     logging.info(f'[ERROR] --- the columns of the bulk dataframe [{ds}] are not in the correct format! Skppping.')
            else:
                for i in ['Tapestri_format', 'condensed_format']:
                    if i in bulk_dfs[ds].columns:
                        bulk_dfs[ds].index = bulk_dfs[ds][i] # <--- set Tapestri_format as the index column
                    break
        logging.info(f'[INFO] --- successfully loaded reference data {ds}')
        try:
            bulk_dfs[ds] = bulk_dfs[ds][(bulk_dfs[ds]['t_depth'] > 0) & (bulk_dfs[ds]['t_alt_count'] > 0)] # filter out vars with zero depth/alt_read
            bulk_dfs[ds]['AF'] = bulk_dfs[ds]['t_alt_count'] / bulk_dfs[ds]['t_depth']
            # bulk_dfs[ds] = bulk_dfs[ds][bulk_dfs[ds]['AF'] > 0] 
            
        except KeyError:
            logging.info(f'[INFO] --- AF cannot be calculated for {ds}')
            continue

    # -- b. matched bulk cohort 
    if bulk_info['matched_bulk_cohort'] is not None:
        if 'maf' not in bulk_info['matched_bulk_cohort'].keys():
            raise ValueError(f'[ERROR] --- for matched_bulk_cohort, {i} needs to be specified!')
        else:
            ds = 'matched_bulk_cohort'
            if 'tumor_sample_name' not in bulk_info[ds].keys():
                logging.info(f'[WARNING] --- for matched_bulk_cohort, tumor_sample_name is not specified. Not sample specific AF will be added.')
            bulk_dfs[ds] = pd.read_csv(bulk_info[ds]['maf'], sep="\t")

            # get condensed variant format
            if not ('Tapestri_format' or 'condensed_format') in bulk_dfs[ds].columns:
                try:
                    bulk_dfs[ds]['Tapestri_format'] = bulk_dfs[ds]['CHROM'] + ":" + bulk_dfs[ds]['POS'] + ":" + bulk_dfs[ds]['REF'] + "/" + bulk_dfs[ds]['ALT']
                    
                except KeyError:
                    logging.info(f'[ERROR] --- the columns of the bulk dataframe [matched_bulk_cohort] are not in the correct format! Skppping.')
            bulk_dfs[ds].index = bulk_dfs[ds]['Tapestri_format'] # <--- set Tapestri_format as the index column
            bulk_dfs[ds] = bulk_dfs[ds][(bulk_dfs[ds]['t_depth'] > 0) & (bulk_dfs[ds]['t_alt_count'] > 0)] # filter out vars with zero depth/alt_read
            bulk_dfs[ds]['AF'] = bulk_dfs[ds]['t_alt_count'] / bulk_dfs[ds]['t_depth']
            logging.info(f'[INFO] --- successfully loaded [matched_bulk_cohort]. The DataFrame has {len(bulk_dfs["matched_bulk_cohort"].index.unique().tolist())} unique variants.')
            if bulk_info[ds]['annotation_summary_method'] is None:
                annotation_summary_method = 'max' # by default, the max vaf of all samples in cohort is annotated
            else:
                annotation_summary_method = bulk_info[ds]['annotation_summary_method']
    else:
        bulk_dfs['matched_bulk_cohort'] = None
        logging.info(f'[INFO] --- matched_bulk_cohort is not given. Skipping.')

    # -- c. bulk normal        
    if bulk_info['matched_bulk_normal'] is not None:
        normal_sample_name = bulk_info['matched_bulk_normal']["normal_sample_name"]
        bulk_dfs['matched_bulk_normal'] = read_vcf_to_df(bulk_info['matched_bulk_normal']["vcf"], normal_sample_name)

        logging.info(f'[INFO] --- successfully loaded [matched_bulk_normal]. The DataFrame has {len(bulk_dfs["matched_bulk_normal"].index.unique().tolist())} unique variants.')
        # AF calculation and index setting already in the `read_vcf_to_df` function
        # bulk_dfs['matched_bulk_normal']['AF'] = bulk_dfs['matched_bulk_normal']['t_alt_count'] / bulk_dfs['matched_bulk_normal']['t_depth']
    else:
        bulk_dfs['matched_bulk_normal'] = None
        logging.info("[INFO] --- matched_bulk_normal is not given. Skipping.")

    # -- d. Tapestri PoN
    if bulk_info['tap_pon'] is not None:
        tap_pon = bulk_info['tap_pon']
        if type(tap_pon) is str:
            # @HZ 11/30/2022: superset PoN df
            logging.info(f'using SUPERSET PoN from: {tap_pon}')
            try: 
                pon_df = pd.read_csv(tap_pon, index_col=0)
                if not check_matrix_format(pon_df, CONDENSED_SNV_FORMAT):
                    raise ValueError
            except:
                logging.error('superset PoN is not in the correct format! Skipping...')
                # sys.exit(1)
            
        elif type(tap_pon) is dict:
            # @HZ 11/30/2022: 
            pon_dfs = {}
            for pon_v in tap_pon:
                pon_dfs[pon_v] = pd.read_csv(tap_pon[pon_v], index_col=0)
                # embed()
                if not check_matrix_format(pon_dfs[pon_v], CONDENSED_SNV_FORMAT):
                    logging.info(f'[ERROR] --- the index of {pon_v} is not in the correct condensed SNV format [chr]:[pos]:[ref]/[alt]. Exiting.')
                    sys.exit(1)
                # pon_dfs[pon_v]['median_sc_mut_prev'] = pon_dfs[pon_v].median(axis=1) @HZ 2022-10-31 changed to mean
                
                logging.info(f'[INFO] --- successfully loaded Tapestri PoN-{pon_v}. The DataFrame has {len(pon_dfs[pon_v].index.unique().tolist())} unique variants.')
        else:
            raise ValueError(f'[ERROR] --- tap_pon needs to be a string or a dictionary!')
    else:
        pon_dfs = None
        logging.warning('[WARNING] --- tap_pon not given')

    # @HZ 11/30/2022: 
    # -- e. Tapestri base-qual blacklist
    if args.bq_info_csv is not None:
        blacklist_var_prev_df = pd.read_csv(args.bq_info_csv, index_col=0,header=None, names=['bq_sc_prev'])
        if not check_matrix_format(blacklist_var_prev_df, CONDENSED_SNV_FORMAT):
            # exit(1)
            logging.warning("[WARNING] --- the index of base-qual blacklist df is not in the correct condensed SNV format chr[chr_number]:[pos]:[ref]/[alt]. Trying to fix by appending `chr` to it...")
            blacklist_var_prev_df.index = 'chr' + blacklist_var_prev_df.index.astype(str)
            if not check_matrix_format(blacklist_var_prev_df, CONDENSED_SNV_FORMAT):
                raise ValueError(f'[ERROR] --- the index cannot be fixed into the right format! Exiting')
                sys.exit(1)

        logging.info(f'[INFO] --- successfully loaded Tapestri base-qual blacklist. The DataFrame has {len(blacklist_var_prev_df.index.unique().tolist())} unique variants.')

    elif ('tap_base-qual_blacklist' in bulk_info) and (bulk_info['tap_base-qual_blacklist'] is not None):
        # try to fetch from bulk_info_yaml
        blacklist_var_prev_f = bulk_info['tap_base-qual_blacklist']
        blacklist_var_prev_df = pd.read_csv(blacklist_var_prev_f, index_col=0, header=None, names=['bq_sc_prev'])
        if not check_matrix_format(blacklist_var_prev_df, CONDENSED_SNV_FORMAT):
            logging.warning("[WARNING] --- the index of base-qual blacklist df is not in the correct condensed SNV format chr[chr_number]:[pos]:[ref]/[alt]. Trying to fix by appending `chr` to it...")
            blacklist_var_prev_df.index = 'chr' + blacklist_var_prev_df.index.astype(str)
            if not check_matrix_format(blacklist_var_prev_df, CONDENSED_SNV_FORMAT):
                raise ValueError(f'[ERROR] --- the index cannot be fixed into the right format! Exiting')
                sys.exit(1)
        logging.info(f'[INFO] --- successfully loaded Tapestri base-qual blacklist. The DataFrame has {len(blacklist_var_prev_df.index.unique().tolist())} unique variants.')

    else:
        logging.warning('[WARNING] --- tap_basequal_blacklist not given! It is a required input for annotation. Exiting')
        sys.exit(1)

    # 2. ----- annotate VAR_DF -----

    # save the annotation stats to a dict
    annotation_stats = {}

    # helper functions
    def fetch_var_info(var: str, ref_df: pd.DataFrame, col_name: str) -> float:
        """
        fetch the AF of a variant from a reference dataframe
        
        Args:
            var: string, the variant in the format CHROM:POS:REF/ALT
            ref_df: pd.DataFrame, the reference dataframe
            col_name: string, the column info to extract
        Returns:
            value: float, the value of the column for the variant
        """
        if var in ref_df.index:
            return ref_df.at[var, col_name]
        else:
            return np.nan
    fetch_var_info_v = np.vectorize(fetch_var_info, excluded = ['ref_df'])

    def _compare_var_dfs(df_target, df_ref, ref_data_type, ref_data_name, col_of_interest = None, summary_method = None) -> Tuple[pd.DataFrame, int, int]:
        """
        compare two dataframes and return the number of variants in each dataframe

        Args:
            df_target: pd.DataFrame, the target dataframe
            df_ref: pd.DataFrame, the reference dataframe
            ref_data_type: string, the data type of the reference dataframe. This will be added to the output dataframe column name
            ref_data_name: string, the name of the reference dataframe. This will be added to the output dataframe
            col_of_interest: string, the column of interest. If None, the presence/absence will be compared. Default: None
            summary_method: func, the summary method to use if multiple values are retrieved by df_ref.loc[var, col_of_interest]. Default: None

        Returns:
            # df_target_anntoated: pd.DataFrame, the target dataframe with the annotation
            ref_unique_var: int, the number of variants in the reference dataframe
            ref_unique_vars_covered: int, the number of variants in the reference dataframe covered by the target dataframe
        """
        ref_unique_count = len(df_ref.index.unique().tolist())
        ref_unique_count_covered = 0
        
        annotated_df = df_target.copy()
        for var_i in df_ref.index.unique().tolist():

            # get column of interest if present; else get presence/absence
            if col_of_interest is not None and col_of_interest in df_ref.columns:
                ann = df_ref.at[var_i, col_of_interest]
            else:
                # logging.info('[WARNING] --- column of interest not given/not present in df_ref column; using presence/absence.')
                ann = var_i in df_ref.index
            
            # if multiple records are present, summarize
            if type(ann) == pd.Series:
                if len(ann) > 1 and summary_method is None:
                    raise ValueError(f'[ERROR] --- variant {var_i} not unique in {ds}')
                if len(ann) > 1 and summary_method is not None:
                    ann = summary_method(ann)

            if not var_i in df_target.index:
                continue
            else:
                annotated_df.at[var_i, f'{ref_data_type}-{ref_data_name}-{col_of_interest}'] = ann
                ref_unique_count_covered += 1

        return annotated_df, ref_unique_count, ref_unique_count_covered

    # (1). bulk data
    for ds in bulk_dfs:
        if bulk_dfs[ds] is None:
            continue

        # given the knowledge that the majority of vars in VAR_DF are artifacts, we can start from the bulk dfs instead of the VAR_DF to be annotated
        if ds == 'matched_bulk_cohort':
            # first, retrieve only the tumor sample
            if bulk_info[ds]['tumor_sample_name'] is not None:
                try:
                    sample_spec_df = bulk_dfs[ds].loc[
                        bulk_dfs[ds]['Tumor_Sample_Barcode'] == bulk_info[ds]['tumor_sample_name']
                        ]
                    VAR_DF, bulk_tumor_unique_vars_count, bulk_tumor_unique_vars_count_covered = _compare_var_dfs(
                        df_target=VAR_DF, 
                        df_ref=sample_spec_df, 
                        ref_data_type='bulk',
                        ref_data_name='matched_bulk_tumor',
                        col_of_interest='AF',
                        )
                    # embed()
                    annotation_stats['bulk_tumor_unique_vars_count'] = bulk_tumor_unique_vars_count
                    annotation_stats['bulk_tumor_unique_vars_count_covered_by_Tapestri'] = bulk_tumor_unique_vars_count_covered
                    logging.info(f"[INFO] --- finished annotating matched_bulk_tumor from matched_bulk_cohort")
                except KeyError:
                    logging.info(f'[ERROR] --- Failed to retrieve {bulk_info["matched_bulk_cohort"]["tumor_sample_name"]}! Tumor sample-specific bulk AF not added')
            
            # second, get cohort-wide AF
            if annotation_summary_method == 'max':
                summary_method = np.max
            elif annotation_summary_method == 'min':
                summary_method = np.min
            elif annotation_summary_method == 'mean':
                summary_method = np.mean
            else:
                raise ValueError(f'[ERROR] --- Unknown annotation summary method: {annotation_summary_method}')

            VAR_DF, bulk_cohort_unique_vars_count, bulk_cohort_unique_vars_count_covered = _compare_var_dfs(
                df_target=VAR_DF,
                df_ref=bulk_dfs[ds],
                ref_data_type='bulk',
                ref_data_name=ds,
                col_of_interest='AF',
                summary_method = summary_method
                )
            # embed()
            annotation_stats['bulk_cohort_unique_vars_count'] = bulk_cohort_unique_vars_count
            annotation_stats['bulk_cohort_unique_vars_count_covered_by_Tapestri'] = bulk_cohort_unique_vars_count_covered
            logging.info(f"[INFO] --- finished annotating matched_bulk_cohort") 
        else:
            # logging.info(f'[INFO] --- working on [{ds}]')
            # each variant record should be unique in matched_bulk_tumor and broad_wes_pon
            VAR_DF, bulk_unique_var_count, bulk_unique_var_count_covered = _compare_var_dfs(
                df_target=VAR_DF,
                df_ref=bulk_dfs[ds],
                ref_data_type='bulk',
                ref_data_name=ds,
                col_of_interest='AF',
                )
            if ds == 'matched_bulk_tumor' or 'matched_bulk_normal':
                annotation_stats[f'{ds}_unique_vars_count'] = bulk_unique_var_count
                annotation_stats[f'{ds}_unique_vars_count_covered_by_Tapestri'] = bulk_unique_var_count_covered
            logging.info(f"[INFO] --- finished annotating {ds}") 
    
    # (2). Tapestri PoN
    if pon_df is not None:
        VAR_DF["PoN-superset-8-normals-occurence"] = VAR_DF.index.map(lambda x: pon_df.at[x, '8-normals-occurence'] if x in pon_df.index else np.nan).astype(float)
        VAR_DF["PoN-superset-8-normals-median_sc_prev"] = VAR_DF.index.map(lambda x: pon_df.at[x, '8-normals-median_sc_prev'] if x in pon_df.index else np.nan).astype(float)
        logging.info("[INFO] --- finished annotating SUPERSET PoN")
        # annotation_stats[f'{pon_vc}_unique_vars_count'] = len(pon_df.index.unique().tolist())
        annotation_stats['num_SNVs_present_in_ANY_normal_sample'] = int((VAR_DF["PoN-superset-8-normals-occurence"] > 0).sum())

    elif pon_dfs is not None:
        for pon_vc in pon_dfs:
            # embed()
            # VAR_DF[f'median_sc_mut_prev-{pon_vc}'] = VAR_DF.index.map(lambda x: pon_dfs[pon_vc].at[x, 'median_sc_mut_prev'] if x in pon_dfs[pon_vc].index else np.nan)
            VAR_DF[f'PoN-{pon_vc}-mean_sc_mut_prev'] = VAR_DF.index.map(lambda x: pon_dfs[pon_vc].at[x, 'mean_num_cells_in_8_normals'] if x in pon_dfs[pon_vc].index else np.nan).astype(float)
            logging.info(f"[INFO] --- finished annotating {pon_vc}")
            annotation_stats[f'{pon_vc}_unique_vars_count'] = len(pon_dfs[pon_vc].index.unique().tolist())
            annotation_stats[f'{pon_vc}_unique_vars_count_covered_by_Tapestri'] = int((VAR_DF[f'PoN-{pon_vc}-mean_sc_mut_prev'] > 0).sum())
            # embed()
    else:
        logging.warning('[WARNING] --- No PoN data given! The PoN filtered H5 would be the same as unfiltered')
    
    # (3). Tapestri base-qual blacklist (sample-specific)
    if blacklist_var_prev_df is not None:
        VAR_DF["blacklist-base_qual-sc_prev"] = VAR_DF.index.map(lambda x: blacklist_var_prev_df.at[x, 'bq_sc_prev'] if x in blacklist_var_prev_df.index else np.nan).astype(float)
        logging.info("[INFO] --- finished annotating base-qual blacklist")
        # annotation_stats[f'{pon_vc}_unique_vars_count'] = len(pon_df.index.unique().tolist())
        annotation_stats['num_SNVs_flagged_as_base-qual_in>3_cells'] = int((VAR_DF["blacklist-base_qual-sc_prev"] >= 3).sum())

    # 3. ----- add info to H5 -----
    for col_name in VAR_DF.columns:
        sample_obj.dna.add_col_attr(col_name, VAR_DF[col_name].values)

    # 4. ----- write output -----
    # VAR_DF
    if args.output_annotated_snv_csv is not None:
        VAR_DF.to_csv(args.output_annotated_snv_csv, index=True, header=True)
        logging.info(f'[INFO] --- Writing annotated SNV file to {args.output_annotated_snv_csv}')
    elif args.output_dir is not None:
        VAR_DF.to_csv(OUTPUT_DIR / f'{SAMPLE_NAME}.bulk_annotated_SNV.csv', index=True, header=True)
        logging.info(f'[INFO] --- Writing annotated SNV file to {OUTPUT_DIR / f"{SAMPLE_NAME}.bulk_annotated_SNV.csv"}')
    else:
        logging.warning('[WARNING] --- Please provide output SNV file name or output directory.')
        sys.exit(1)
       
    # optional outputs # V2
    BULK_ANNOTATED_H5 = args.output_bulk_annotated_h5
    if BULK_ANNOTATED_H5 is not None:
        mio.save(sample_obj, BULK_ANNOTATED_H5)
    else:
        logging.warning('[WARNING] --- Not writing bulk annotated H5.')
        # BULK_ANNOTATED_H5 = OUTPUT_DIR / f'{SAMPLE_NAME}_DNA_CNV_m2_f.bulk_annotated.h5'
    PONV1_FILTERED_H5 = args.output_ponv1_filtered_h5
    if PONV1_FILTERED_H5 is not None:
        try:
            voi = sample_obj.dna.ids()[
                np.isnan(sample_obj.dna.col_attrs['mean_sc_mut_prev-tap_pon_v1']) | 
                ~np.isnan(sample_obj.dna.col_attrs['AF-matched_bulk_normal'])
            ]
            sample_obj.dna = sample_obj.dna[sample_obj.dna.barcodes(), voi] # filter out PoN v1 variants
            mio.save(sample_obj, PONV1_FILTERED_H5)
        except KeyError:
            logging.info('[ERROR] --- No PoN v1 data given! Skipping writing PoN filtered H5.')
    else:
        logging.warning('[WARNING] --- Not writing PoN v1 filtered H5.')
    #     PONV1_FILTERED_H5 = OUTPUT_DIR / f'{SAMPLE_NAME}_DNA_CNV_m2_f.ponv1_filtered.h5'
    

    # # TLOD metrics data
    # VAR_DF.to_csv(TLOD_METRICS_TSV, sep="\t", index=True)
    # VAR_DF.loc[voi, :].to_csv(TLOD_METRICS_PONV1_FILTERED_TSV, sep="\t", index=True)

    # JSON file
    with open(OUTPUT_DIR / f"{SAMPLE_NAME}-annotation_stats.json", "w") as out_f:
        json.dump(annotation_stats, out_f)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type = str, help='sample name')
    parser.add_argument('--input_h5', type = str, help='input_h5')
    parser.add_argument('--bq_info_csv', type = str, help='base-quality info output by the custom pipeline.', default=None)
    parser.add_argument('--bulk_info_yaml', type = str, help='bulk_info_yaml')
    parser.add_argument('--output_dir', type = str, help='output_dir to write the output dataframe.', default=None)
    parser.add_argument('--output_annotated_snv_csv', type = str, help='optional output: annotated_snv_csv')
    parser.add_argument('--output_bulk_annotated_h5', type = str, help='optional output: bulk_annotated_h5', default=None)
    parser.add_argument('--output_ponv1_filtered_h5', type = str, help='optional output: ponv1_filtered_h5')
    parser.add_argument('--log_file', type = str, default=None, help='log file')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)