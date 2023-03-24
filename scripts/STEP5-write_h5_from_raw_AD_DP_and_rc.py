import argparse, sys, os, re, json, shutil
from pathlib import Path
import logging
import pandas as pd
import numpy as np
import mosaic.io as mio
from h5.create import create_cnv_assay
from h5.data import Assay, H5Writer
from h5.constants import (
    AF,
    ALT,
    CHROM,
    DATE_CREATED,
    DEFAULT_SAMPLE,
    DNA_ASSAY,
    DNA_READ_COUNTS_ASSAY,
    DP,
    GQ,
    ID,
    NGT,
    POS,
    QUAL,
    REF,
    RGQ,
    SAMPLE,
)
from tea.format import CONDENSED_SNV_FORMAT, check_matrix_format
from tea.utils import get_simple_timestamp
from tea.annotate import annotate_snv_amplicon_coverage

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

from IPython import embed

def compute_af(ad: np.ndarray, dp: np.ndarray) -> np.ndarray:
    """Compute allele frequency
    Args:
        ad: allele depth
        dp: read depth
    Returns:
        allele frequency
    """
    # suppress true_divide warning
    with np.errstate(invalid="ignore"):
        return np.nan_to_num(ad / dp * 100)

def create_variant_col_attrs(var_list) -> pd.DataFrame:
    """Create a variant dataframe from a list of variants
    Args:
        var_list: array-like variants in the format CHROM:POS:REF/ALT
    Returns:
        variant dataframe
    """
    var_df = pd.DataFrame(index = var_list)

    # handle b37 --> hg19 `chr` notation issue
    if not check_matrix_format(var_df, CONDENSED_SNV_FORMAT):
        logging.warning("[WARNING] --- the index of the SNV df is not in the correct condensed SNV format chr[chr_number]:[pos]:[ref]/[alt]. Trying to fix by appending `chr` to it...")
        var_df.index = 'chr' + var_df.index.astype(str)
        if not check_matrix_format(var_df, CONDENSED_SNV_FORMAT):
            raise ValueError(f'[ERROR] --- the index cannot be fixed into the right format! Exiting')
            sys.exit(1)
        
    var_df['REF'] = var_df.index.str.split(':').str[2].str.split('/').str[0]
    var_df['ALT'] = var_df.index.str.split(':').str[2].str.split('/').str[1]
    var_df['CHROM'] = var_df.index.str.split(':').str[0]
    var_df['POS'] = var_df.index.str.split(':').str[1]
    return var_df

def add_amplicon_metadata(cnv_assay, amplicons):
    ca = cnv_assay.col_attrs
    chrom = np.full((cnv_assay.shape[1],), "", dtype=object)
    start_pos = np.full((cnv_assay.shape[1],), 0, dtype=int)
    end_pos = np.full((cnv_assay.shape[1],), 0, dtype=int)

    for _, amplicon in amplicons.iterrows():
        matching_ids = ca[ID] == amplicon['amplicon']
        chrom[matching_ids] = amplicon['CHROM'].strip("chrCHR")
        start_pos[matching_ids] = amplicon['start_pos']
        end_pos[matching_ids] = amplicon['end_pos']

    cnv_assay.add_col_attr('CHROM', chrom)
    cnv_assay.add_col_attr('start_pos', start_pos)
    cnv_assay.add_col_attr('end_pos', end_pos)

def main(args):
    # ----- io ----------------------------------------
    sample_name = args.sample_name
    metadata_json = args.metadata_json
    if metadata_json is not None:
        with open(metadata_json, 'r') as f:
            metadata = json.load(f)
        metadata['date_created'] = get_simple_timestamp()
        for i in metadata:
            # resolve nested dict
            if isinstance(metadata[i], dict):
                metadata[i] = json.dumps(metadata[i])
    else:
        metadata = {
            'sample_name': np.array([sample_name]),
            'date_created': get_simple_timestamp(),
            "pass2_route": "bcf_mpileup; raw AD, DP extraction",
        }

    # --- inputs for creating the DNA assay
    AD_merged_df = pd.read_csv(args.AD_merged_df, index_col=0)
    DP_merged_df = pd.read_csv(args.DP_merged_df, index_col=0)
    AD_merged_df = AD_merged_df[~AD_merged_df.index.duplicated()]
    DP_merged_df = DP_merged_df[~DP_merged_df.index.duplicated()]
    # sanity check: the single cell indices and variant order need to be uinified
    # assert AD_merged_df.shape == DP_merged_df.shape, f"AD and DP layers have different shapes; AD: {AD_merged_df.shape}, DP: {DP_merged_df.shape}"
    assert np.all(AD_merged_df.index == DP_merged_df.index), f"AD and DP layers have different indices; AD: {AD_merged_df.index}, DP: {DP_merged_df.index}"
    panel_insert_file = args.panel_insert_file # this need to be UCSC format ('chr1' instead of '1')

    # --- inputs for creating the CNV assay
    if args.read_counts_tsv is not None: # @HZ 03/08/2022: allows for the creation of H5 without the CNV assay
        read_counts_tsv = args.read_counts_tsv
        barcode_num_map_f = args.barcode_num_map_f
        barcode_num_map_df = pd.read_csv(barcode_num_map_f, sep='\t', index_col= 1, header=0, names = ['numerical_index'])
        barcode_num_map = barcode_num_map_df['numerical_index'].to_dict()
        panel_amplicon_file = args.panel_amplicon_file
    else:
        logging.warning(f'[WARNING] --- no read_counts_tsv provided, skipping CNV assay creation')
    
    # --- outputs
    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
    output_h5 = args.output_h5
    if output_h5 is None:
        output_h5 = str(output_dir / f'{sample_name}.mpileup.h5')
    read_counts_tsv_renamed = args.read_counts_tsv_renamed
    if read_counts_tsv_renamed is None:
        read_counts_tsv_renamed = str(output_dir / f'{sample_name}.per_amplicon_read_counts.tsv')
    # --------------------------------------------------
    logging.info(f'[STEP5-mpileup_route] starting to write to output file: {output_h5}')

    ###############################
    # part 1 ----- create DNA assay
    ###############################
    dna = Assay.create(DNA_ASSAY)

    alt_read_count = AD_merged_df.T.values
    dp = DP_merged_df.T.values

    dna.add_layer('alt_read_count', alt_read_count)
    # embed()
    dna.add_layer(DP, dp)
    dna.add_layer(AF, compute_af(alt_read_count, dp))

    # add row and column attributes
    var_df = create_variant_col_attrs(AD_merged_df.index)
    dna.add_row_attr("barcode", AD_merged_df.columns.values)
    dna.add_row_attr(SAMPLE, np.tile(np.array([sample_name]), len(AD_merged_df.columns)))
    dna.add_col_attr(ID, var_df.index.values)
    for name, series in var_df.iteritems():
        dna.add_col_attr(name, series.values)
    logging.info(f'[STEP5-mpileup_route] finished creating DNA assay.')
    # add metadata
    for name, value in metadata.items():
        dna.add_metadata(name, value)
    logging.info(f'[STEP5-mpileup_route] added DNA metadata from: {metadata_json}.')

    ###############################
    # part 2 ----- create CNV assay
    ###############################
    if args.read_counts_tsv is not None:
        rc_df = pd.read_csv(read_counts_tsv, sep='\t', index_col=0)
        # rename cell barcode --> cell numerical index
        rc_df.index = rc_df.index.map(barcode_num_map)
        rc_df.sort_index(inplace=True)
        # save renamed matrix to output
        rc_df.to_csv(read_counts_tsv_renamed, sep='\t')

        cnv = Assay.create(DNA_READ_COUNTS_ASSAY)
        cnv.add_layer("read_counts", rc_df.values)

        cnv.add_row_attr('barcode', rc_df.index.values)
        cnv.add_col_attr('id', rc_df.columns.values)

        for name, value in metadata.items():
            cnv.add_metadata(name, value)

        cnv.add_row_attr(SAMPLE, np.array([sample_name] * cnv.shape[0]))
        
        # cnv.row_attrs[SAMPLE] = np.tile(np.array([sample_name]), cnv.shape[0])
        # # ^^^^ this circumvents the bug in MB's h5.create.cnv.create_cnv_assay() where the sample_name is inferrd from the CNV.METADATA
        logging.info(f'[STEP5-mpileup_route] finished creating CNV assay.')
        logging.info(f'[STEP5-mpileup_route] CNV assay row_attr `barcode` shape: {cnv.row_attrs["barcode"].shape}')
        logging.info(f'[STEP5-mpileup_route] CNV assay row_attr `sample` shape: {cnv.row_attrs[SAMPLE].shape}')
        amplicons = pd.read_csv(
                    panel_amplicon_file,
                    sep="\t",
                    lineterminator="\n",
                    names=['CHROM', 'start_pos', 'end_pos', 'amplicon'],
                    dtype={'start_pos': int, 'end_pos': int},
        )
        add_amplicon_metadata(cnv, amplicons)
        logging.info(f'[STEP5-mpileup_route] added amplicon metadata from: {panel_amplicon_file}')
    else:
        cnv = None

    
    ###############################
    # part 3 ----- add both to H5
    ###############################

    # assays = [dna, cnv]
    assays = [ assay for assay in [dna, cnv] if assay is not None ]

    if Path(output_h5).is_file():
        logging.warning(f'[STEP5-mpileup_route] output file already exists: {output_h5}; overwriting!')
        os.remove(output_h5)
    with H5Writer(output_h5) as writer:
        for assay in assays:
            writer.write(assay)
        logging.info('[STEP5-mpileup_route] finished adding assays to output H5 file.')

    # add amplicon coverage for each SNV
    sample = mio.load(output_h5, name = sample_name, raw = False)
    os.remove(output_h5)
    temp_dir = output_dir / 'temp'
    Path(temp_dir).mkdir(parents=True, exist_ok=True)
    snv_panel_bed_isec_df = annotate_snv_amplicon_coverage(
        sample, 
        insert_bed = panel_insert_file, 
        working_dir = temp_dir,
        )
    sample.dna.add_col_attr('amplicon_coverage', snv_panel_bed_isec_df.values)
    mio.save(sample, output_h5)
    shutil.rmtree(str(temp_dir))
    logging.info('[STEP5-mpileup_route] finished adding amplicon coverage info for SNV layer in H5.')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type=str, required=True, help='sample name')
    parser.add_argument('--metadata_json', type=str, default=None, help='metadata json file')
    parser.add_argument('--AD_merged_df', type=str, help='Dataframe containing sc-AD values. The index needs to be in the format: CHR:POS:REF/ALT', required=True)
    parser.add_argument('--DP_merged_df', type=str, help='Directory containing sc-DP values. The index needs to be in the format: CHR:POS', required=True)
    parser.add_argument('--read_counts_tsv', type=str, help='read counts tsv output by the MB pipeline.', default=None)
    parser.add_argument('--barcode_num_map_f', type=str, help='sc barcode to numerical index map file', default=None)
    parser.add_argument('--panel_insert_file', type=str, help='insert file', required=True)
    parser.add_argument('--panel_amplicon_file', type=str, help='amplicon file')
    parser.add_argument('--output_dir', type=str, help='output directory', required=True)
    parser.add_argument('--output_h5', type=str, help='path to write output H5.', default=None)
    parser.add_argument('--read_counts_tsv_renamed', type=str, help='path to write output read counts tsv, with cell barcodes renamed to numerical indices.', default=None)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)