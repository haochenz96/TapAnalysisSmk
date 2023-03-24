import argparse
from pathlib import Path
import re
import os, sys
import yaml

def main(args):
    output_dir = Path(args.wd)
    analysis_config_yaml = args.analysis_config
    with open(analysis_config_yaml, 'r') as f:
        analysis_config = yaml.safe_load(f)
    
    sample_i = 'RA17_22-27_2'

    ALT_df = sample_objs[sample_i].dna.get_attribute('alt_read_count', constraint='row')
    # ALT_df.columns = [i + '_v' for i in ALT_df.columns]
    DP_df = sample_objs[sample_i].dna.get_attribute('DP', constraint='row')
    # DP_df.columns = [i + '_t' for i in DP_df.columns]

    AF_DP_df = pd.DataFrame(index = ALT_df.index)
    for i in ALT_df.columns:
        AF_DP_df[i + '_v'] = ALT_df[i]
        AF_DP_df[i + '_t'] = DP_df[i]

    scarlet_df = pd.DataFrame(index = sample_objs[sample_i].dna.barcodes())
    scarlet_df['c'] = sample_objs[sample_i].cnv.row_attrs['clone_id-NB_EM_nclones=5']

    scarlet_df = pd.concat([scarlet_df, AF_DP_df], axis=1)
    scarlet_df.index.name = 'cell_id'

    scarlet_df.to_csv(output_dir / f'{sample_i}-scarlet-rc.csv', index = True, header = True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--wd', type=str, help='working directory')
    parser.add_argument('--prefix', type=str, help='prefix for output files', default='')
    parser.add_argument('--analysis_config', type=str, help='analysis_config.yaml')


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)