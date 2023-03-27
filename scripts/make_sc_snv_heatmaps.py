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

    # ====== load cn_assignment_df ======
    cn_assignment_df = pd.read_csv(cn_assignment_df, index_col = 0)
    if not 'clone_id' in cn_assignment_df.columns:
        try: # try to add clone_id column
            cn_assignment_df['clone_id'] = cn_assignment_df['cluster_id']
        except:
            raise ValueError(f'[ERROR] `clone_id`/`cluster_id` column not found in CN clone assignment file!')


    unique_cluster_ids_sorted = np.sort(np.unique(cn_assignment_df['clone_id']))
    unique_cluster_ids_sorted_named = [ f"CN_clone-{clone_i}" for clone_i in unique_cluster_ids_sorted ]
    cn_clone_palette = dict(zip(unique_cluster_ids_sorted_named, np.array(px.colors.qualitative.Set3)[unique_cluster_ids_sorted]))

    for sample_i in sample_names:
        # genotype
        if not "NGT" in sample_objs[sample_i].dna.layers:
            sample_objs[sample_i].dna.genotype_variants(
                min_dp = 8,
                min_alt_read = 3,
                assign_low_conf_genotype = True,
                )

        # add cn_clone info
        if not sample_i in cn_assignment_df.index:
            raise ValueError(f"{sample_i} not in cn_assignment_df.index")
        cn_assignment_dict = cn_assignment_df.loc[sample_i,:].set_index('cell_barcode').to_dict()['clone_id']

        sample_objs[sample_i].dna.row_attrs['label'] = np.array(list(map(lambda x: f"CN_clone-{cn_assignment_dict[x]}", sample_objs[sample_i].dna.barcodes())))
        sample_objs[sample_i].dna.set_palette(cn_clone_palette)

        num_cells = sample_objs[sample_i].dna.shape[0]
        print(f'[INFO] {sample_i} has {num_cells} cells.')
        if args.write_cn_clone_added_h5s:
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

    sc_heatmaps_params = analysis_config['sc_heatmaps_params']

    bq_prev_threshold = sc_heatmaps_params['bq_prev_threshold']
    if not type(bq_prev_threshold) is float: # single-value
        raise ValueError(f"bq_prev_threshold must be float, not {type(bq_prev_threshold)}")
    # embed()
    mut_prev_threshold = sc_heatmaps_params['mut_prev_threshold']
    if not type(mut_prev_threshold) is list:
        mut_prev_threshold = [mut_prev_threshold]

    topic = sc_heatmaps_params['topic']
    if not type(topic) is str: # single-value
        raise ValueError(f"topic must be str, not {type(topic)}")
    try: 
        func_only = sc_heatmaps_params['func_only']
        func_only = bool(func_only)
    except KeyError:
        func_only = False
    except TypeError:
        func_only = False

    attribute = sc_heatmaps_params['attribute']
    if not type(attribute) is list:
        attribute = [attribute]

    for mut_prev_i in mut_prev_threshold:
        (wd / 'sc_heatmaps' / topic / f"mut_prev>={mut_prev_i}").mkdir(exist_ok=True, parents=True)
        for attribute_i in attribute:

            voi_union = set()
            voi_count_union = {}
            ann_map_union = {}

            # ====== for each sample_i, filter SNVs and get a union set for plotting  ======
            # embed()
            for sample_i in sample_names:
                num_cells = sample_objs[sample_i].dna.shape[0]
                mask = isNaN(cravat_dfs[sample_i][('bulk_comparison', 'bulk-broad_wes_pon-AF')]) & \
                ~(cravat_dfs[sample_i][('PoN_comparison','PoN-superset-8-normals-occurence')] >= 4) & \
                ~(cravat_dfs[sample_i][('blacklist_comparison', 'blacklist-base_qual-sc_prev')] >= bq_prev_threshold * num_cells) & \
                ~(cravat_dfs[sample_i][('blacklist_comparison', 'blacklist-base_qual-sc_prev')] >= cravat_dfs[sample_i][('Tapestri_result', 'sc_mut_prev')]) & \
                (cravat_dfs[sample_i][('Tapestri_result', 'sc_mut_prev')] >= mut_prev_i * num_cells)
                if func_only:
                    mask = mask & ~cravat_dfs[sample_i][('Variant Annotation', 'Sequence Ontology')].isin(NONFUNC_SO)
                voi = cravat_dfs[sample_i].index[mask]
                voi_mut_prev = Counter(cravat_dfs[sample_i].loc[voi, ('Tapestri_result', 'sc_mut_prev')].to_dict())

                print(f'{sample_i}: {len(voi)} / {len(cravat_dfs[sample_i])} SNVs are kept (min prev: {mut_prev_i})')
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

            voi_sorted = sorted(voi_count_union, key=voi_count_union.get, reverse=True)
            with open(wd / 'sc_heatmaps' / topic / f"mut_prev>={mut_prev_i}" / f'{topic}-{attribute_i}-voi.txt', 'w') as f:
                f.write(f"timestamp\n")
                f.write('=====================\n')
                f.write(f"bq_prev_threshold: {bq_prev_threshold}\n")
                f.write(f"mut_prev_threshold: {mut_prev_i}\n")
                f.write(f"topic: {topic}\n")
                f.write(f"func_only: {func_only}\n")
                for voi_i in voi_sorted:
                    f.write(f"{voi_i}\t{voi_count_union[voi_i]}\n")

            for sample_i in sample_names:

                fig = plot_snv_clone(
                    sample_objs[sample_i],
                    sample_name=sample_i,
                    story_topic = 'major_driver_snvs',
                    voi = voi_sorted,
                    attribute = attribute_i,
                    ann_map = ann_map_union
                )

                fig.write_image(
                    wd / 'sc_heatmaps' / topic / f"mut_prev>={mut_prev_i}" / f'{sample_i}-{topic}-{attribute_i}.png',
                    scale = 3
                )

                # fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--wd', type=str, help='working directory')
    parser.add_argument('--analysis_config', type=str, help='analysis_config.yaml')
    parser.add_argument('--write_cn_clone_added_h5s', type=bool, help='write_cn_clone_added_h5s', default=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)