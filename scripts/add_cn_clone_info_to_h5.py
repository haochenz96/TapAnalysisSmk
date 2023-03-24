# assemble CN-clone info into 

import mosaic.io as mio
import pandas as pd
import numpy as np
from pathlib import Path
import re
import os, sys
import argparse
import plotly.express as px
import json
from tea.format import TAPESTRI_BARCODE_FORMAT
from IPython import embed


def main(args):

    # ===== IO =====
    cohort_name = args.cohort_name
    sample_name = args.sample_name
    raw_h5 = args.raw_h5
    cn_clone_cell_assignments_csv = args.cn_clone_cell_assignments_csv
    cn_clone_profiles_csv = args.cn_clone_profiles_csv
    output_h5 = args.output_h5
    # ===============

    # 1 ----- read h5 into mosaic.assay object -----
    sample_obj = mio.load(raw_h5, name = sample_name, raw=False)

    # 2 ----- standard genotyping -----
    sample_obj.dna.genotype_variants(
        min_dp = 8,
        min_alt_read = 3,
        assign_low_conf_genotype = True,
        )
    voi = sample_obj.dna.ids()[
        (sample_obj.dna.get_attribute('mut_filtered').sum(axis=0) >= 3)
        ] # at least confidently called in 3 cells

    sample_obj.dna = sample_obj.dna[
            sample_obj.dna.barcodes(),
            voi
            ]

    # 3 ----- add CN clone info -----
    cn_clone_cell_assignments_df = pd.read_csv(
        cn_clone_cell_assignments_csv, 
        index_col = 0,
        header = 0
        ) # ncells x 3 [sample, cell_barcode, clone_id]
    if not sample_name in cn_clone_cell_assignments_df.index:
        raise ValueError(f'[ERROR] {sample_name} not found in {cn_clone_cell_assignments_csv}')
    try:
        cn_clone_cell_assignment_df_ordered = cn_clone_cell_assignments_df.loc[sample_name].set_index('cell_barcode').loc[sample_obj.cnv.barcodes()]
    except KeyError:
        raise KeyError(f'[ERROR] some barcodes in {sample_name} are not found in {cn_clone_cell_assignments_csv}')
    if not 'clone_id' in cn_clone_cell_assignments_df.columns:
        try: # try to add clone_id column
            cn_clone_cell_assignments_df['clone_id'] = cn_clone_cell_assignments_df['cluster_id']
        except:
            raise ValueError(f'[ERROR] `clone_id`/`cluster_id` column not found in {cn_clone_profiles_csv}')
    
    cn_clone_profiles_df = pd.read_csv(
        cn_clone_profiles_csv,
        index_col = 0,
        header = 0,
    ) # nclones x namplicons

    # sanity check to make sure the number of clones in the cn_clone_profiles_df matches the number of clones in the cell_assignment_df
    for clone_i in cn_clone_cell_assignment_df_ordered.clone_id.unique():
        assert clone_i in cn_clone_profiles_df.index, f'[ERROR] clone {clone_i} in cn_clone_cell_assignment_df_ordered is not found in cn_clone_profiles_df'

    sc_amplicon_ploidy_df = cn_clone_profiles_df.loc[cn_clone_cell_assignment_df_ordered.clone_id.values]
    ploidy_layer = sc_amplicon_ploidy_df.reindex(
                index = sample_obj.cnv.barcodes(),
                columns = sample_obj.cnv.ids()
            ).fillna(0) # make sure the barcodes and amplicons are in the same order as in the input H5 
            # reindexing fills in non-existing amplicons with NaN
            # we then fill them with 0's to avoid plotting error later
    sample_obj.cnv.add_layer(f'ploidy-NB-EM-homdel', ploidy_layer.values)
    print('-'*50)
    print(f'[INFO] added ploidy layer ploidy-NB_EM to {sample_name} input H5')

    cn_clone_labels_ordered = cn_clone_cell_assignment_df_ordered['clone_id']
    assert cn_clone_labels_ordered.isna().any() == False, f'[ERROR] some barcodes are not assigned to any clone'
    sample_obj.cnv.add_row_attr(f'clone_id-NB-EM-homdel', cn_clone_labels_ordered.values)
    
    unique_cn_clone_labels = cn_clone_profiles_df.index.unique()
    cn_clone_palette = dict(zip(unique_cn_clone_labels, np.array(px.colors.qualitative.Set3)[unique_cn_clone_labels]))
    sample_obj.cnv.row_attrs['label'] = cn_clone_labels_ordered.values
    sample_obj.cnv.set_palette(cn_clone_palette)

    # copy to DNA layer as well
    sample_obj.dna.add_row_attr('clone_id-NB-EM-homdel', cn_clone_labels_ordered.values)
    sample_obj.dna.row_attrs['label'] = cn_clone_labels_ordered.values
    sample_obj.dna.set_palette(cn_clone_palette)

    print(f'[INFO] added single-cell CN clone-id assignment clone_id-NB_EM to {sample_name} input H5')
    
    # 4 ----- save -----
    Path(output_h5).parent.mkdir(parents=True, exist_ok=True)
    mio.save(sample_obj, output_h5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort_name', type=str, help='cohort name')
    parser.add_argument('--sample_name', type=str, help='sample name')
    parser.add_argument('--raw_h5', type=str, help='path to raw H5 output by TapVarCall piepline')
    parser.add_argument('--cn_clone_cell_assignments_csv', type=str, help='dataframe containing cell assignments', )
    parser.add_argument('--cn_clone_profiles_csv', type=str, help='CN clone profiles')
    parser.add_argument('--output_h5', type=str, help='output h5 path')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)