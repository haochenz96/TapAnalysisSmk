# @HZ 06/20/2022
# script to query variants from SNV_df with CRAVAT
# inputs:
#   - an SNV_df file
# outputs:
#   - a file with the variants and their CRAVAT annotations

# first, install customized-mosaic and tea from Github

from pathlib import Path
import mosaic.io as mio
from tea.utils import *
from tea.plots import *
from tea.cravat import *
import requests
import base64
import json
import sys, argparse, yaml

# for testing
# SAMPLE_NAME = 'M13-1'
# H5 = '/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/Caitlin/M13-1_combined/OUTPUTS_from_m2_f/M13-1_combined_DNA_CNV_m2_f.bulk_annotated.ponv1_filtered.h5'
# __output_dir = Path('/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/Caitlin/M13-1_combined/OUTPUTS_from_m2_f/annotations')
# CRAVAT_OUTPUT = __output_dir / 'M13-1_ponv1_filtered_CRAVAT_output.txt'
# LOG_FILE = __output_dir / 'STEP6_CRAVAT_log.txt'
# OPENCRAVAT_USERNAME = 'zhangh5@mskcc.org'
# OPENCRAVAT_PASSWORD = 'Hfhd68MdTtHn5UE'
# CRAVAT_ANNOTATORS = ["clinvar", "gnomad3", "chasmplus", "chasmplus_PAAD","civic","cosmic","dbsnp","dbsnp_common","clinvar","gnomad3",'thousandgenomes','go','ndex']
# CRAVAT_INPUT_PATH = __output_dir / 'M13-1_ponv1_filtered_CRAVAT_input.txt'

def main(args):
    # ----- io -----
    SAMPLE_NAME = args.sample_name
    # --- inputs
    snv_df = pd.read_csv(args.input_snv_csv, index_col=0)
    # --- outputs
    if args.output_dir is None:
        CRAVAT_OUTPUT = args.cravat_output
        CRAVAT_INPUT_PATH = CRAVAT_OUTPUT.parent / f'{SAMPLE_NAME}_CRAVAT_input.txt'
        CRAVAT_OUTPUT_CLEANED = CRAVAT_OUTPUT.parent / f'{SAMPLE_NAME}_CRAVAT_output_cleaned.txt'
    else:
        __output_dir = Path(args.output_dir)
        __output_dir.mkdir(exist_ok=True, parents=True)
        CRAVAT_OUTPUT = __output_dir / f'{SAMPLE_NAME}_CRAVAT_output.txt'
        CRAVAT_INPUT_PATH = __output_dir / f'{SAMPLE_NAME}_CRAVAT_input.txt'
        CRAVAT_OUTPUT_CLEANED = __output_dir / f'{SAMPLE_NAME}_CRAVAT_output_cleaned.txt'
    # --- params
    ###############################
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
            # filename = LOG_FILE,
            handlers=[logging.FileHandler(LOG_FILE, 'a'), logging.StreamHandler()],
            level = logging.INFO,
            force=True,
            format='%(asctime)s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            )
    ###############################  
    cravat_settings_yaml = args.cravat_settings_yaml
    with open(cravat_settings_yaml, "r") as f:
        cravat_settings = yaml.safe_load(f)

    try:
        OPENCRAVAT_USERNAME = cravat_settings['username']
        OPENCRAVAT_PASSWORD = cravat_settings['password']
        CRAVAT_ANNOTATORS = cravat_settings['annotators']
    except KeyError:
        logging.error('Please check your cravat_settings.yaml file')
        sys.exit(1)
    
    if not CRAVAT_OUTPUT.is_file():
        # 1. establish connection to OpenCRAVAT server:
        new_session = requests.Session()

        # log in with credentials 
        cred_str = base64.b64encode(bytes(OPENCRAVAT_USERNAME + ":" + OPENCRAVAT_PASSWORD, 'utf-8')).decode('utf-8')
        reply = new_session.get(
            'https://run.opencravat.org/server/login', 
            headers={'Authorization': 'Basic ' + cred_str})
        # original cred_str ----->> base64.b64encode(b'zhangh5@mskcc.org:Hfhd68MdTtHn5UE').decode() 

        if not reply.json() == 'success':
            logging.warning('CRAVAT login failed!')
            exit(1)

        # @HZ 04/01/2023: this is the new method based on snv_df
        vars_of_interest = snv_df.index[
            (snv_df['Tapestri_result-sc_mut_prev'] >= 3)
        ].tolist()
        logging.info(f'{len(vars_of_interest)} variants written to CRVAT input')
        mut_prev_column = snv_df.loc[vars_of_interest, 'Tapestri_result-sc_mut_prev']
        write_cravat_input(vars_of_interest, CRAVAT_INPUT_PATH, SAMPLE_NAME, mut_prev_column)
        logging.info(f'{len(vars_of_interest)} variants of interest written to CRAVAT input')
    
        cravat_params = json.dumps({
            "annotators": CRAVAT_ANNOTATORS, 
            "reports": ["text"], 
            "assembly": "hg19", 
            "note": f"{SAMPLE_NAME}-analysis",
            })

        post = new_session.post(
            'https://run.opencravat.org/submit/submit', 
            files={'file_0': open(CRAVAT_INPUT_PATH)}, 
            data={'options': cravat_params
                }
        )
        
        cravat_job_id = post.json()['id']
        get_cravat_output(new_session, cravat_job_id, CRAVAT_OUTPUT)
        logging.info(f'Finished CRAVAT analysis for sample -- {SAMPLE_NAME}.')
    else:
        logging.warning('CRAVAT output file already exists! Skipping...')

    # @HZ 04/01/2023: clean and format the CRAVAT output
    # 1. read in the CRAVAT report
    with read_until_match(
        in_file = CRAVAT_OUTPUT, 
        match_pattern = (lambda line: "#CRAVAT Report" in line), 
        num_occurences = 2
        ) as f:
        
        f.seek(0) # !!!reset cursor!!!
        cravat_df = pd.read_csv(f, sep='\t', skiprows=4, header=[0,1], low_memory=False)
        logging.info(f'[INFO] {SAMPLE_NAME} CRAVAT report [ {CRAVAT_OUTPUT} ] read in successfully.')

    # 2. clean and format the CRAVAT report
    cravat_df = clean_and_format_cravat_df(cravat_df)
    logging.info(f'[INFO] {SAMPLE_NAME} CRAVAT report cleaned and formatted successfully (see https://github.com/haochenz96/tea/blob/3b9d40642869608f0fff3fe2827f03271675dff1/src/tea/cravat.py#L102).')

    # 3. add key information from SNV_DF
    for data_type in snv_df.columns:
        if 'bulk' in data_type:
            logging.info(f'[INFO] adding bulk comparisons for {data_type}')
            if type(snv_df[data_type]) == str: # True/False
                snv_df[data_type] = bool(snv_df[data_type])
            cravat_df.loc[:, ('bulk_comparison', data_type)] = snv_df[data_type]
        elif 'PoN' in data_type:
            logging.info(f'[INFO] adding PoN comparisons for {data_type}')
            if type(snv_df[data_type]) == str:
                snv_df[data_type] = bool(snv_df[data_type])
            cravat_df.loc[:, ('PoN_comparison', data_type)] = snv_df[data_type]
        elif 'blacklist' in data_type:
            logging.info(f'[INFO] adding blacklist comparisons for {data_type}')
            cravat_df.loc[:, ('blacklist_comparison', data_type)] = snv_df[data_type]
        else:
            continue
    # adjust the order of bulk & pon comparison columns
    for lv1_col_name in ['bulk_comparison', 'PoN_comparison', 'blacklist_comparison']:
        cols = cravat_df.pop(lv1_col_name)
        for i in range(cols.shape[1]):
            cravat_df.insert(i, (lv1_col_name, cols.columns[i]), cols.iloc[:,i])
    # logging.info(f'[INFO] pushed forward the columns for {ds}')

    # # insert TLOD max value into the DataFrame
    # i = np.where(cravat_df.columns == ('Tapestri_result', 'sc_mut_prev'))[0][0] + 1
    # tlod_col = cravat_df.pop(('Tapestri_result', 'sc_max_TLOD'))
    # cravat_df.insert(int(i), ('Tapestri_result', 'sc_max_TLOD'), tlod_col)

    # important: fetch variant annotation before any sorting
    gene_HGVSp = cravat_df.index.map(
        lambda x: 
        cravat_df.loc[x, ('Variant Annotation', 'Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation', 'Protein Change')] if not isNaN(cravat_df.loc[x, ('Variant Annotation','Protein Change')])
        else cravat_df.loc[x, ('Variant Annotation','Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation','Sequence Ontology')]
    )
    # add gene_HGVSp column
    cravat_df.insert(0, ('var_ann', 'Gene_HGVSp'), gene_HGVSp)

    # sort by Tapestri sc mutational prevalence
    cravat_df.sort_values(by = ('Tapestri_result', 'sc_mut_prev'), ascending=False, inplace=True)
    cravat_df.to_csv(CRAVAT_OUTPUT_CLEANED, sep='\t')
    logging.info(f'[INFO] saved cleaned CRAVAT report to {CRAVAT_OUTPUT_CLEANED}')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type = str, help='sample name')
    parser.add_argument('--input_snv_csv', type = str, help='input_snv_csv, output by STEP6-annotate_h5_with_bulk.py')
    parser.add_argument('--cravat_settings_yaml', type = str, required=True, help='bulk_info_yaml')
    parser.add_argument('--cravat_output', type = str, help='output CRAVAT report')
    parser.add_argument('--output_dir', type = str, help='output_dir; overrides the output file name if specified.')
    parser.add_argument('--log_file', type = str, default=None, help='log file')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)