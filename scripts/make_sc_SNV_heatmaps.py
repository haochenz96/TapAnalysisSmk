import argparse
from pathlib import Path
import re
import os, sys
import yaml

import mosaic.io as mio
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from tea.cravat import NONFUNC_SO
from tea.format import isNaN
from tea.plots import plot_snv_clone
from datetime import datetime
from collections import Counter

from tea.cravat import get_technical_artifact_mask
from tea.format import CONDENSED_SNV_FORMAT, check_matrix_format
from tea.utils import get_simple_timestamp
timestamp = get_simple_timestamp()
from IPython import embed

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
    cn_assignment_f = analysis_config['cn_assignment_df']

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

    for sample_i in sample_names:
        # genotype
        if not "NGT" in sample_objs[sample_i].dna.layers:
            sample_objs[sample_i].dna.genotype_variants(
                min_dp = 8,
                min_alt_read = 3,
                assign_low_conf_genotype = True,
                )

    # ====== load cn_assignment_df ======
    cn_err = 0
    if cn_assignment_f is None:
        cn_err = 1
        print('[WARNING] No CN clone assignment file provided. All cells in the SNV heatmap will be in one cluster...')
    else:
        cn_assignment_df = pd.read_csv(cn_assignment_f, index_col = 0)
        print(f'[INFO] Loaded CN clone assignment file {cn_assignment_f}.')
        if not 'clone_id' in cn_assignment_df.columns:
            try: # try to add clone_id column
                cn_assignment_df['clone_id'] = cn_assignment_df['cluster_id']
            except:
                raise ValueError(f'[ERROR] `clone_id`/`cluster_id` column not found in CN clone assignment file!')


        unique_cluster_ids_sorted = np.sort(np.unique(cn_assignment_df['clone_id'].astype(int)))
        unique_cluster_ids_sorted_named = [ f"CN_clone-{clone_i}" for clone_i in unique_cluster_ids_sorted ]
        # embed()
        cn_clone_palette = dict(zip(unique_cluster_ids_sorted_named, np.array(px.colors.qualitative.Set3)[unique_cluster_ids_sorted]))

        for sample_i in sample_names:
            # add cn_clone info
            if not sample_i in cn_assignment_df.index:
                raise ValueError(f"{sample_i} not in cn_assignment_df.index")
            cn_assignment_dict = cn_assignment_df.loc[sample_i,:].set_index('cell_barcode').to_dict()['clone_id']

            try:
                sample_objs[sample_i].dna.row_attrs['label'] = np.array(list(map(lambda x: f"CN_clone-{int(cn_assignment_dict[x])}", sample_objs[sample_i].dna.barcodes())))
                sample_objs[sample_i].dna.set_palette(cn_clone_palette)
            except KeyError:
                print (f'[ERROR] {sample_i} not found in cn_assignment_df.index!')
                cn_err = 1

            num_cells = sample_objs[sample_i].dna.shape[0]
            print(f'[INFO] {sample_i} has {num_cells} cells.')
            if args.write_cn_clone_added_h5s and cn_err == 0:
                out_dir = wd / 'cn_clone_added_h5s'
                out_dir.mkdir(exist_ok=True, parents=True)
                try:
                    mio.save(sample_objs[sample_i], str(out_dir / f"{sample_i}_cn_clone_added.h5"))
                except FileExistsError:
                    print(f'[WARNING] {out_dir / f"{sample_i}_cn_clone_added.h5"} already exists! Overwriting...')
                    try: 
                        # delete and try again
                        os.remove(str(out_dir / f"{sample_i}_cn_clone_added.h5"))
                        mio.save(sample_objs[sample_i], str(out_dir / f"{sample_i}_cn_clone_added.h5"))
                    except:
                        print(f'[ERROR] {out_dir / f"{sample_i}_cn_clone_added.h5"} already exists and cannot be overwritten!')

    # ====== set up for plotting ======
    # mkdir, get plotting params
    (wd / 'sc_heatmaps').mkdir(exist_ok=True, parents=True)
    
    ####################################
    # ### ----- SNV selection ----- ####
    ####################################

    snv_selection_params = analysis_config['snv_selection_params']
    
    # 1. mutational prevalence
    mut_prev_threshold = snv_selection_params['mut_prev_threshold']
    if not type(mut_prev_threshold) is list:
        mut_prev_threshold = [mut_prev_threshold]
    
    # 2. technical artifact filters
    bq_prev_threshold = snv_selection_params['bq_prev_threshold']
    if bq_prev_threshold is not None and type(bq_prev_threshold) is not float: # single-value
        raise ValueError(f"bq_prev_threshold must be float, not {type(bq_prev_threshold)}")

    normals_occurences = snv_selection_params['normals_occurences'] if 'normals_occurences' in snv_selection_params else 3 # <--- default to 3
    
    if 'ado_threshold' in snv_selection_params and snv_selection_params['ado_threshold'] is not None:
        print(f"[INFO] filtering out SNVs with ADO > {snv_selection_params['ado_threshold']} in ANY sample...")
        ado_threshold = snv_selection_params['ado_threshold']
    else:
        ado_threshold = None

    # 3. functional SNVs
    topic = snv_selection_params['topic']
    if not type(topic) is str: # single-value
        raise ValueError(f"topic must be str, not {type(topic)}")
    try: 
        func_only = snv_selection_params['func_only']
        func_only = bool(func_only)
    except KeyError:
        func_only = False
    except TypeError:
        func_only = False

    # 4. germline SNPs
    germline_attrs = {}
    for germline_attr_i in ['select_hom_germline_snps_af', 'rescue_1000genome_af']:
        if not germline_attr_i in snv_selection_params:
            germline_attrs[germline_attr_i] = False
        else:
            germline_attrs[germline_attr_i] = snv_selection_params[germline_attr_i]

    # 5. plotting params
    attribute = snv_selection_params['attribute']
    if not type(attribute) is list:
        attribute = [attribute]

    # whitelist snvs
    if 'whitelist_snvs' in snv_selection_params and snv_selection_params['whitelist_snvs']:
        whitelist_snvs = set(snv_selection_params['whitelist_snvs'])
    else:
        whitelist_snvs = []
    
    # blacklist snvs
    if 'blacklist_snvs' in snv_selection_params and snv_selection_params['blacklist_snvs']:
        if type(snv_selection_params['blacklist_snvs']) is list:
            blacklist_snvs = set(snv_selection_params['blacklist_snvs'])
        else:
            try:
                blacklist_snv_df = pd.read_csv(snv_selection_params['blacklist_snvs'], index_col=0)
                if not check_matrix_format(blacklist_snv_df, CONDENSED_SNV_FORMAT):
                    raise ValueError(f"[ERROR] blacklist_snvs file not in the correct format (index column needs to be in condensed SNV format).")
                blacklist_snvs = set(blacklist_snv_df.index)
            except:
                raise ValueError(f"[ERROR] blacklist_snvs should either be a list of SNVs or a path to a CSV file whose index is the list of SNVs.")
    else:
        blacklist_snvs = []

    for mut_prev_i in mut_prev_threshold:
        print(f"""
              ===== 
              [INFO] mut_prev_i = {mut_prev_i} 
              =====
              """)
        
        (wd / 'sc_heatmaps' / topic / f"mut_prev={mut_prev_i}").mkdir(exist_ok=True, parents=True)

        # ----- @HZ 07/17/2023 filter T>C artifacts that are rare -----
        if 'filter_TtoC_artifact' in snv_selection_params:
            try:
                filter_TtoC_artifact = snv_selection_params['filter_TtoC_artifact']['filter'] 
                filter_TtoC_artifact_lower_thres = snv_selection_params['filter_TtoC_artifact']['lower_thres']
                filter_TtoC_artifact_upper_thres = snv_selection_params['filter_TtoC_artifact']['upper_thres']
            except KeyError:
                raise ValueError(f"[ERROR] filter_TtoC_artifact must have keys ['filter', 'lower_thres', 'upper_thres']")
            else:
                if filter_TtoC_artifact_lower_thres >= filter_TtoC_artifact_upper_thres:
                    raise ValueError(f"[ERROR] filter_TtoC_artifact_lower_thres must be strictly smaller than filter_TtoC_artifact_upper_thres")
        else:
            if mut_prev_i < 0.01:
                print(f"[INFO] mut_prev_i is lower than default upper_thres (0.01) for T>C filter. The filter will be applied.")
                filter_TtoC_artifact = True
                filter_TtoC_artifact_lower_thres = mut_prev_i
                filter_TtoC_artifact_upper_thres = 0.01 
            else:
                print(f"[WARNING] mut_prev_i is higher than default upper_thres (0.01) for T>C filter. The filter will not be applied.")
                filter_TtoC_artifact = False
        # -----------------------------------------------------------------

        voi_union = set()
        voi_count_union = {}
        ann_map_union = {}
        bulk_germline_vars = set()
        bulk_somatic_vars = set()
        TtoC_artifact_blacklist = set()
        ado_blacklist = set()

        # ====== for each sample_i, filter SNVs and get a union set for plotting  ======
        # embed()
        for sample_i in sample_names:
            num_cells = sample_objs[sample_i].dna.shape[0]
            
            # this masks on PoN's (rescuing 1000genome SNVs) and bq_prev_threshold
            mask = get_technical_artifact_mask(cravat_dfs[sample_i], num_cells = num_cells, bq_prev_threshold = bq_prev_threshold, normals_pon_occurence=normals_occurences, rescue_1000genome_af = germline_attrs['rescue_1000genome_af'], filter_broad_wes_pon = False)
            
            # # filters on mut_prev_threshold
            # mask = mask & (cravat_dfs[sample_i][('Tapestri_result', 'sc_mut_prev')] >= mut_prev_i * num_cells)
            
            # filter on T>C artifacts that are rare
            if filter_TtoC_artifact:
                tc_mask = (cravat_dfs[sample_i].index.str.endswith('T/C')) & (cravat_dfs[sample_i][('Tapestri_result', 'sc_mut_prev')] >= filter_TtoC_artifact_lower_thres * num_cells) & (cravat_dfs[sample_i][('Tapestri_result', 'sc_mut_prev')] <= filter_TtoC_artifact_upper_thres * num_cells)
                print(f"[INFO] {sample_i} has {tc_mask.sum()} T>C artifacts that are rare. They will be filtered out.")
                mask = mask & ~tc_mask
                TtoC_artifact_blacklist = TtoC_artifact_blacklist.union(
                    set(cravat_dfs[sample_i].index[tc_mask].tolist())
                )

            # filters on functional SNVs
            if func_only:
                mask = mask & ~cravat_dfs[sample_i][('Variant Annotation', 'Sequence Ontology')].isin(NONFUNC_SO)

            voi = cravat_dfs[sample_i].index[mask].tolist()

            # filters on mut_prev_threshold
            prev_filtered_vars = sample_objs[sample_i].dna.ids()[
                sample_objs[sample_i].dna.get_attribute("mut_filtered", constraint="row").sum(axis=0) >= (mut_prev_i * num_cells)
            ]
            # take intersection
            voi = [ v for v in voi if v in prev_filtered_vars ]

            # blacklist SNVs that have too high ADO in any sample
            if ado_threshold is not None:
                ado_high_vars = sample_objs[sample_i].dna.ids()[
                    (sample_objs[sample_i].dna.get_attribute('NGT',constraint='row') == 3).sum(axis=0) > (ado_threshold*num_cells)
                ]
                
                voi = [ v for v in voi if v not in ado_high_vars ]
                print(f"[INFO] In sample {sample_i}, filtered {len([v for v in voi if v in ado_high_vars])} SNVs due to high ADO")
                ado_blacklist = ado_blacklist.union(set(ado_high_vars))

            voi = list(set(voi).union(whitelist_snvs))
            voi_mut_prev = Counter(cravat_dfs[sample_i].loc[voi, ('Tapestri_result', 'sc_mut_prev')].to_dict())

            print(f'{sample_i}: {len(voi)} / {len(cravat_dfs[sample_i])} SNVs are kept (min prev: {mut_prev_i})')
            ann = cravat_dfs[sample_i].loc[voi, :].index.map(
                lambda x: 
                cravat_dfs[sample_i].loc[x, ('Variant Annotation', 'Gene')] + ' ' + cravat_dfs[sample_i].loc[x, ('Variant Annotation', 'Protein Change')] if not isNaN(cravat_dfs[sample_i].loc[x, ('Variant Annotation','Protein Change')])
                else cravat_dfs[sample_i].loc[x, ('Variant Annotation','Gene')] + ' ' + cravat_dfs[sample_i].loc[x, ('Variant Annotation','Sequence Ontology')]
            )
            ann_map = dict(zip(voi, ann))

            # ====== save bulk annotated vars ======

            # select SNVs detected in matched bulk normal
            # embed()
            # sys.exit()
            try:
                bulk_normal_vars = cravat_dfs[sample_i].index[(cravat_dfs[sample_i][('bulk_comparison', 'bulk-matched_bulk_normal-AF')] > 0)]
            except KeyError:
                print(f'bulk normal annotation not found in CRAVAT DF for {sample_i}')
            else:
                if len(bulk_normal_vars) == 0:
                    print(f'[WARNING] No bulk normal SNVs detected in {sample_i}')
                bulk_germline_vars.update(bulk_normal_vars)

            # select SNVs detected in matched bulk cohort
            try:
                bulk_cohort_vars = cravat_dfs[sample_i].index[(cravat_dfs[sample_i][('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] > 0)]
            except KeyError:
                print(f'bulk tumor annotation not found in CRAVAT DF for {sample_i}')
            else:
                if len(bulk_cohort_vars) == 0:
                    print(f'[WARNING] No bulk cohort SNVs detected in {sample_i}')
                bulk_somatic_vars.update(bulk_cohort_vars)

            # try:
            #     bulk_annotation_df = cravat_dfs[sample_i].loc[voi, [
            #         ('bulk_comparison', 'bulk-matched_bulk_normal-AF'), 
            #         ('bulk_comparison', 'bulk-matched_bulk_cohort-AF'),
            #         ]] > 0
            #     for var_i in bulk_annotation_df.index:
            #         if bulk_annotation_df.loc[var_i, ('bulk_comparison', 'bulk-matched_bulk_normal-AF')] == True:
            #             bulk_germline_vars.add(var_i)
            #         elif bulk_annotation_df.loc[var_i, ('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] == True:
            #             bulk_somatic_vars.add(var_i)
            # except KeyError:
            #     print(f'bulk annotation not found in CRAVAT DF for {sample_i}')

            # ===== get a union of all samples SNVs =====
            voi_union = voi_union.union(set(voi))
            voi_count_union.update(voi_mut_prev)
            ann_map_union.update(ann_map)

        # select homozygous germline SNVs, if specified
        if germline_attrs['select_hom_germline_snps_af']:
            alt_combined = pd.concat(
                [ sample_objs[sample_i].dna.get_attribute('alt_read_count', constraint='row') for sample_i in sample_names ], axis=0
            )
            dp_combined = pd.concat(
                [ sample_objs[sample_i].dna.get_attribute('DP', constraint='row') for sample_i in sample_names ], axis=0 
            )
            overall_af = alt_combined.sum(axis=0) / dp_combined.sum(axis=0)
            germline_hom_snps_from_tapestri = [x for x in voi_union if overall_af[x] > germline_attrs['select_hom_germline_snps_af']]
            print(f"[INFO] selected {len(germline_hom_snps_from_tapestri)} homozygous germline SNVs (AF > {germline_attrs['select_hom_germline_snps_af']})")
        else:
            germline_hom_snps_from_tapestri = []

        # remove SNVs that are blacklisted
        print(f"[DEBUG] {len(voi_union)} SNVs before blacklist filtering")
        voi_union = voi_union.difference(TtoC_artifact_blacklist)
        print(f"[DEBUG] {len(voi_union)} SNVs after TtoC blacklist filtering")
        voi_union = voi_union.difference(ado_blacklist)
        print(f"[DEBUG] {len(voi_union)} SNVs after ADO blacklist filtering")
        voi_union = voi_union.difference(blacklist_snvs)
        print(f"[DEBUG] {len(voi_union)} SNVs after manual blacklist filtering")

        # embed()
        voi_sorted = sorted(voi_union, key=voi_count_union.get, reverse=True)
        # label germline and bulk somatic vars:
        def __annotate_snvs(var, germline_hom_vars, other_bulk_germline_vars, bulk_somatic_vars):
            if var in germline_hom_vars:
                return "germline_HOM"
            elif var in other_bulk_germline_vars:
                return "germline_HET"
            elif var in bulk_somatic_vars:
                return "bulk_somatic"
            else:
                return "NA"
            
        with open(wd / 'sc_heatmaps' / topic / f"mut_prev={mut_prev_i}" / f'{topic}-voi.txt', 'w') as f:
            voi_df = pd.DataFrame.from_dict(voi_count_union, orient='index', columns=['mut_prev']).loc[voi_sorted]
            voi_df['HGVSp'] = voi_df.index.map(lambda x: ann_map_union[x])
            voi_df['annotation'] = voi_df.index.map(
                lambda x: __annotate_snvs(x, germline_hom_snps_from_tapestri, bulk_germline_vars, bulk_somatic_vars)
            )
            voi_df.index.name = 'condensed_format'
            f.write(f"# {timestamp}\n")
            f.write(f"# cohort_name: {cohort_name}")
            f.write(f"# sample_names: {sample_names}")
            f.write('# ======  [snv_selection_params] ======\n')
            f.write(f"# attribute: {attribute}\n")
            f.write(f"# mut_prev_threshold: {mut_prev_i}\n")
            f.write(f"# bq_prev_threshold: {bq_prev_threshold}\n")
            f.write(f"# topic: {topic}\n")
            f.write(f"# func_only: {func_only}\n")
            f.write(f"# normals_occurences: {normals_occurences}\n")
            f.write(f"# ado_threshold: {ado_threshold}\n")
            f.write(f"# germline_attrs: {germline_attrs}\n")
            f.write(f"# filter_TtoC_artifact: {snv_selection_params['filter_TtoC_artifact']}\n")
            f.write(f"# whitelist_snvs: {whitelist_snvs}\n")
            if 'blacklist_snvs' in snv_selection_params:
                f.write(f"# blacklist_snvs: {snv_selection_params['blacklist_snvs']}\n")
            else:
                f.write(f"# blacklist_snvs: None\n")
            f.write('# =====================================\n')

            voi_df.to_csv(f, sep='\t', index=True, header=True)
            
        if args.filter_snv_only:
            # stop after filtering
            sys.exit()

        # ====== plot ======
        # ===== highlight vars with bulk annotation ====='
        # highlight vars
        germline_hom_var_col = '#00cc66' # green
        germline_het_var_col = '#2693ff' # blue
        somatic_var_col = '#ff0000' # red

        for var_i in ann_map_union:
            # germline
            if var_i in germline_hom_snps_from_tapestri:
                ann_map_union[var_i] = f'<span style="color:{germline_hom_var_col};">' + ann_map_union[var_i] + '</span>'
            elif var_i in bulk_germline_vars:
                ann_map_union[var_i] = f'<span style="color:{germline_het_var_col};">' + ann_map_union[var_i] + '</span>'
            elif var_i in bulk_somatic_vars:
                ann_map_union[var_i] = f'<span style="color:{somatic_var_col};">' + ann_map_union[var_i] + '</span>'
            else:
                pass

        for sample_i in sample_names:
            for attribute_i in attribute:
                fig = plot_snv_clone(
                    sample_objs[sample_i],
                    sample_name=sample_i,
                    story_topic = 'high_conf_mutations-{topic}',
                    voi = voi_sorted,
                    attribute = attribute_i,
                    ann_map = ann_map_union
                )
                fig.update_layout(
                    width = max(len(voi_sorted) * 20, 600),
                )
                fig.write_image(
                    wd / 'sc_heatmaps' / topic / f"mut_prev={mut_prev_i}" / f'{sample_i}-{topic}-{attribute_i}.pdf',
                )

                # fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--wd', type=str, help='working directory')
    parser.add_argument('--analysis_config', type=str, help='analysis_config.yaml')
    parser.add_argument('--filter_snv_only', type=bool, help='only filter and write SNVs; not proceed into making heatmaps.', default=False)
#     parser.add_argument('--combined_heatmap', type=bool, help='make a combined heatmap for all samples', default=False)
    parser.add_argument('--write_cn_clone_added_h5s', type=bool, help='write_cn_clone_added_h5s', default=False)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)