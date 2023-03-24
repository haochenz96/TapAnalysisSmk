import argparse
from pathlib import Path
import re
import os, sys
import yaml

import mosaic.io as mio
from tea.parse import * 
from tea.cravat import NONFUNC_SO
from tea.format import isNaN
from collections import Counter


from tea.plots import plot_var_sc_mut_prev_histogram

def main(args):
    
    # ===== IO =====
    wd = Path(args.wd)
    # cohort_name = args.cohort_name
    # sample_name = args.sample_name
    # cn_clone_added_h5 = args.cn_clone_added_h5  
    # cravat_output = args.cravat_output
    analysis_config_yaml = args.analysis_config
    with open(analysis_config_yaml, 'r') as f:
        analysis_config = yaml.safe_load(f)

    cohort_name = analysis_config['cohort_name']
    sample_names = analysis_config['sample_names']
    h5s = analysis_config['h5s']
    cravat_outputs = analysis_config['cravat_outputs']
    cn_assignment_df = analysis_config['cn_assignment_df']

    # ===============

    sample_objs = {}
    cravat_dfs = {}

    for sample_i in sample_names:
        try: 
            h5 = h5s[sample_i]
            sample_objs[sample_i] = mio.load(h5)
        except KeyError:
            print(f'[ERROR] {sample_i} not found in h5s.')
        try:
            cravat_txt = cravat_outputs[sample_i]
            cravat_dfs[sample_i] = pd.read_csv(
                cravat_txt, sep = '\t', index_col= 0, header = [0,1]
            )
        except KeyError:
            print(f'[ERROR] {sample_i} not found in cravat_outputs.')
        sample_objs[sample_i].dna.genotype_variants(
            min_dp = 8,
            min_alt_read = 3,
            assign_low_conf_genotype = True,
        )

    # ====== get a union of all samples SNVs and annotations =====
    #@HZ hardcoded params
    bq_prev_threshold = 0.005
    func_only = False
    # definition of a functional variant:
    # 1. not in NONFUNC_SO
    # or 
    # 2. not defined as `Benign` or `Benign/Likely benign` or `Likely benign` in cravat_df[("ClinVar", "Clinical Significance")]

    voi_union = set()
    voi_count_union = Counter()
    ann_map_union = {}
    ann_df_union = pd.DataFrame()
    for sample_i in sample_names:
        num_cells = sample_objs[sample_i].dna.shape[0]
        mask = isNaN(cravat_dfs[sample_i][('bulk_comparison', 'bulk-broad_wes_pon-AF')]) & \
        ~(cravat_dfs[sample_i][('PoN_comparison','PoN-superset-8-normals-occurence')] >= 4) & \
        ~(cravat_dfs[sample_i][('blacklist_comparison', 'blacklist-base_qual-sc_prev')] >= bq_prev_threshold*num_cells) & \
        ~(cravat_dfs[sample_i][('blacklist_comparison', 'blacklist-base_qual-sc_prev')] >= cravat_dfs[sample_i][('Tapestri_result', 'sc_mut_prev')]) 
        # (cravat_dfs[sample_i][('Tapestri_result', 'sc_mut_prev')] >= mut_prev_threshold * num_cells)
        if func_only:
            mask = mask & ~cravat_dfs[sample_i][('Variant Annotation', 'Sequence Ontology')].isin(NONFUNC_SO)

        voi = cravat_dfs[sample_i].index[mask]
        voi_mut_prev = Counter(cravat_dfs[sample_i].loc[voi, ('Tapestri_result', 'sc_mut_prev')].to_dict())

        # print(f'{sample_i}: {len(voi)} / {len(cravat_dfs[sample_i])} SNVs are kept (min prev: {mut_prev_threshold})')
        print(f'{sample_i}: {len(voi)} / {len(cravat_dfs[sample_i])} SNVs are kept')
        
        ann = cravat_dfs[sample_i].loc[voi, :].index.map(
            lambda x: 
            cravat_dfs[sample_i].loc[x, ('Variant Annotation', 'Gene')] + ' ' + cravat_dfs[sample_i].loc[x, ('Variant Annotation', 'Protein Change')] if not isNaN(cravat_dfs[sample_i].loc[x, ('Variant Annotation','Protein Change')])
            else cravat_dfs[sample_i].loc[x, ('Variant Annotation','Gene')] + ' ' + cravat_dfs[sample_i].loc[x, ('Variant Annotation','Sequence Ontology')]
        )
        ann_map = dict(zip(voi, ann))
        # ===== get a union of all samples SNVs =====
        voi_union = voi_union.union(set(voi.tolist()))
        voi_count_union.update(voi_mut_prev)
        ann_map_union.update(ann_map)
        ann_df_union = pd.concat([ann_df_union, cravat_dfs[sample_i].loc[voi, :]], axis = 0)
        ann_df_union = ann_df_union.loc[~ann_df_union.index.duplicated(keep = 'first'), :]
    
    # ====== set up for plotting ======
    # mkdir, get plotting params
    out_dir =  (wd / 'snv_mut_prev_distributions')
    out_dir.mkdir(exist_ok=True, parents=True)

    total_number_of_cells = sum([sample_objs[sample_i].dna.shape[0] for sample_i in sample_names])

    var_mut_prev_df_for_plot = pd.DataFrame(index = ann_df_union.index)
    var_mut_prev_df_for_plot['sc_mut_prev'] = ann_df_union.index.map(lambda x: voi_count_union[x])
    var_mut_prev_df_for_plot['functional?'] = ~(ann_df_union[('Variant Annotation','Sequence Ontology')].isin(NONFUNC_SO) | 
        ann_df_union[('ClinVar', 'Clinical Significance')].isin(['Benign', 'Benign/Likely benign', 'Likely benign']))
    var_mut_prev_df_for_plot['functional?'] = var_mut_prev_df_for_plot['functional?'].astype(str)
    var_mut_prev_df_for_plot['var_name_annotated'] = ann_df_union.index.map(ann_map_union.get)

    fig = plot_var_sc_mut_prev_histogram(
        df = var_mut_prev_df_for_plot,
        sample_name = cohort_name,
        split_by = 'functional?',
        color_map = {'True': 'red', 'False': 'blue'},
        num_bins = 100,
        )

    fig.update_layout(font_family="Arial",height = 600, width = 1000,)
    fig.update_yaxes(
        title_font=dict(size = 12, color='black'),
    )
    fig.write_image(
        str(out_dir / f'{cohort_name}_snv_mut_prev_distributions.pdf'),
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--wd', type=str, help='working directory')
    parser.add_argument('--analysis_config', type=str, help='analysis_config.yaml')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)