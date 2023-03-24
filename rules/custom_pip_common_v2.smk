# @HZ 03/09/2022
# prepare all inputs for running custom single-cell variant calling pipeline
# @HZ 10/01/2022
# V2: allow for custom barcode file

from snakemake import load_configfile
from pathlib import Path
import pysam
import pandas as pd
import sys, os

genome_version = config['reference_info']['genome_version']

# ----- specify directories -----
cohort_name = config['sample_info']['cohort_name']
top_dir = Path(config['top_dir'])
if config['working_dir'] == "default":
    print(f"[INFO] working_dir is set to default: {config['top_dir']}/{cohort_name}")
    working_dir = top_dir / cohort_name
else:
    print(f"[INFO] working_dir specified to be: {config['working_dir']}")
    working_dir = Path(config['working_dir'])

print(f'[INFO] --- working directory ---- {working_dir}')

# ----- get current date and time
eastern = tz.timezone('US/Eastern') # <------ uses US eastern time by default
now = datetime.now(tz.utc).astimezone(eastern) 
timestamp = now.strftime("%a %m %d %H:%M:%S %Z %Y")
date_simple = now.strftime("%Y-%m-%d")

# ----- copy the Snakemake and config files over
config_copy = working_dir / f'run_config--{date_simple}.yaml'
snake_file_copy = working_dir / f'snakefile--{date_simple}'
if not config_copy.is_file():
    args = sys.argv
    try:
        config_file_path = args[args.index("--configfile") + 1]
        shutil.copyfile(config_file_path, config_copy)
        print("[INFO] --- copied over the config file.")
    except:
        print("[WARNING] --- Failed to copy over the config file.")
if not snake_file_copy.is_file():
    try:
        snake_file_path = 'workflow/Snakefile'
        shutil.copyfile(snake_file_path, snake_file_copy)
        print("[INFO] --- copied over the snake file.")
    except:
        print("[WARNING] --- Failed to copy over the snake file.")

# ----- get sample and reference info -----
    # @HZ 10/01/2022
    # Note: for this version, the sample_info in the YAML file needs to be in the format:
    # -- cohort_name:
    # -- sample_bams:
    # -- -- sample1: /path/to/sample1.bam
    # -- -- sample2: /path/to/sample2.bam
    # -- -- ...
    # -- sample_barcode_maps:
    # -- -- sample1: /path/to/sample1_barcode_map.tsv
    # -- -- sample2: /path/to/sample2_barcode_map.tsv
    # -- -- ...

# ----- Fetch all single-cell barcodes from each sample's combined BAM ----- 
def fetch_sample_barcode_map(sample_i, bam_file = None, wd = None):
    '''
    Try to fetch the sample barcode maps for valid sc-BAMs to use
    
    '''
    if bam_file is None:
        if wd is None:
            raise(ValueError('Must provide either a BAM file or a working directory to infer the BAM file!'))
        else:
            bam_file = wd / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam'
    
    # write the barcode map file to this directory
    bars_map_file = wd / sample_i / 'reference' / f'{sample_i}.barcode_map.txt'
    if not bars_map_file.is_file():
        print(f'[WARNING] --- {sample_i} barcode map file not found. Fetching barcode map from BAM file.')
        (wd / sample_i / 'reference').mkdir(exist_ok=True, parents=True)
        bars_map = {}
        with pysam.AlignmentFile(bam_file, "rb") as cells_bam:
            # @HZ 04/19/2022 create numerical index for each barcode
            cell_count = 0
            for i in cells_bam.header['RG']:
                cell_count += 1
                bars_map[f'cell_{cell_count}'] = i['SM'] # <-------------------- numerical index to cell barcode map
                #bars.append(i['SM'])

            bars_map_df = pd.DataFrame.from_dict(bars_map, orient='index', columns = ['cell_barcode'])
            bars_map_df.index.name = 'num_index'
            bars_map_df.to_csv(str(bars_map_file), sep='\t', index=True)
        bars_map = bars_map_df.to_dict()['cell_barcode']
    else:
        bars_map_df = pd.read_csv(bars_map_file, sep='\t', index_col=0)
        # sanity check barcode map:
        with pysam.AlignmentFile(bam_file, "rb") as cells_bam:
            barcodes_from_BAM = pd.DataFrame(cells_bam.header['RG'])['SM']

            if not bars_map_df['cell_barcode'].isin(barcodes_from_BAM).all():
                #raise KeyError(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
                sys.exit(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
            elif bars_map_df.shape[0] != barcodes_from_BAM.shape[0]:
                sys.exit(f'[ERROR] --- {sample_i} barcode map file length does not match number of barcodes from BAM! Abort.')
            else:
                print(f'[INFO] --- {sample_i} barcode map file is OK.')
                
    bars_map = bars_map_df.to_dict()['cell_barcode']
    return bars_map

sample_barcode_maps = {} # cohort-wide, each sample has a barcode map
if 'sample_names' in config['sample_info']:
    # V1 method: just infer the BAM path and use all barcodes
    # this also overrides all the manually input BAMs and barcode_maps
    sample_names = config['sample_info']['sample_names']
    for sample_i in sample_names:
        sample_barcode_maps[sample_i] = fetch_sample_barcode_map(
            sample_i,
            config['sample_info']['sample_bams'][sample_i],
            working_dir,
            )
else:
    # V2 method
    assert 'sample_bams' in config['sample_info'].keys(), "Must specify either sample_names or sample_bams in the config file."
    assert 'sample_rc_tsvs' in config['sample_info'].keys(), "Must specify either sample_rc_tsvs in the config file."
    sample_names = config['sample_info']['sample_bams'].keys()
    for sample_i in sample_names:
        (working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam').mkdir(parents=True, exist_ok=True)
        (working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'tsv').mkdir(parents=True, exist_ok=True)
        if not os.path.lexists(str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam')):
            os.symlink(
                config['sample_info']['sample_bams'][sample_i], 
                str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'bam' / f'{sample_i}.tube1.cells.bam')
                )
        if not os.path.lexists(str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'tsv' / f'{sample_i}.tube1.barcode.cell.distribution.tsv')):
            os.symlink(
                config['sample_info']['sample_rc_tsvs'][sample_i], 
                str(working_dir / sample_i / 'tap_pipeline_output' / 'results' / 'tsv' / f'{sample_i}.tube1.barcode.cell.distribution.tsv')
                )

    if config['sample_info']['sample_barcode_maps'] is not None:
        # sanity check:
        assert sample_names == config['sample_info']['sample_barcode_maps'].keys(), '[WARNING] Sample names in sample_bams and sample_barcodes do not match!'
        sample_barcode_maps = config['sample_info']['sample_barcode_maps']

        for sample_i, map_i in sample_barcode_maps.items():
            if map_i is None:
                print(f'[WARNING] No barcode map provided for sample {sample_i}! Will use all barcodes in provided BAM file.')
                sample_barcode_maps[sample_i] = fetch_sample_barcode_map(
                    sample_i,
                    config['sample_info']['sample_bams'][sample_i], 
                    working_dir
                    )
            else:
                bars_map_file_copy = working_dir / sample_i / 'reference' / f'{sample_i}.barcode_map.txt'
                (working_dir / sample_i / 'reference').mkdir(exist_ok=True, parents=True)
                shutil.copyfile(map_i, bars_map_file_copy)
                bars_map_df = pd.read_csv(map_i, sep='\t', index_col=0)
                # sanity check barcode map:
                with pysam.AlignmentFile(
                    config['sample_info']['sample_bams'][sample_i],
                    "rb"
                    ) as cells_bam:
                    barcodes_from_BAM = pd.DataFrame(cells_bam.header['RG'])['SM']
                    if not bars_map_df['cell_barcode'].isin(barcodes_from_BAM).all():
                        #raise KeyError(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
                        sys.exit(f'[ERROR] --- {sample_i} barcode map file does not contain all the barcodes from BAM file! Abort.')
                    # elif bars_map_df.shape[0] != barcodes_from_BAM.shape[0]:
                    #     sys.exit(f'[ERROR] --- {sample_i} barcode map file length does not match number of barcodes from BAM! Abort.')
                    else:
                        print(f'[INFO] --- {sample_i} barcode map file is OK.')
                        
                bars_map = bars_map_df.to_dict()['cell_barcode']
                sample_barcode_maps[sample_i] = bars_map
    else:
        print('[INFO] No barcode map provided for any sample! Will use all barcodes in provided BAM file.')
        for sample_i in sample_names:
            sample_barcode_maps[sample_i] = fetch_sample_barcode_map(sample_i, None, working_dir)


# sample_barcode_maps = fetch_sample_barcode_maps(sample_names, working_dir)
print('[INFO] barcode reading done!')

# (1) get step1 sc_bams output
def get_step1_sc_bams(sample_names, sample_barcode_maps):
    # for the main Snakefile as target output
    # gets all sc_bams for ALL sample
    out = []
    for sample_i in sample_names:
        for cell_num_index in list(sample_barcode_maps[sample_i].keys()):
            out.append(f"{sample_i}/1-sc_bams/{sample_i}_{cell_num_index}.bam")
    return out

def get_step1_sc_bams_by_sample(wildcards):
    # for Snakemake rule use
    # gets all sc_bams for A GIVEN sample
    sample_i = wildcards.sample_name
    out = []
    for cell_num_index in list(sample_barcode_maps[sample_i].keys()):
        out.append(f"{sample_i}/1-sc_bams/{sample_i}_{cell_num_index}.bam")
    return out

# (2) get the single-cell "filter_added.vcf.gz" outputs
def get_step2_filter_added_mutect2_vcfs(sample_names, sample_barcode_maps):
    # for the main Snakefile as TARGET OUTPUT
    out = []
    for sample_i in sample_names:
        for cell_num_index in sample_barcode_maps[sample_i].keys():
            out.append(f"{sample_i}/2-sc_mutect2_call/m2_sc_vcfs_filter_added/{sample_i}_{cell_num_index}_somatic_m2_filter_added.vcf.gz")
    return out

def get_step2_m2_vcfs_by_sample(wildcards):
    # for Snakemake rule use only
    # gets all step2 mutect2 output vcf for A GIVEN sample
    # as input for step2 bcftools filter
    out = []
    for cell_num_index in sample_barcode_maps[wildcards.sample_name].keys():
        out.append(f"{wildcards.sample_name}/2-sc_mutect2_call/m2_sc_vcfs_filter_added/{wildcards.sample_name}_{cell_num_index}_somatic_m2_filter_added.vcf.gz")
    return out

def get_step2_strelka2_vcfs(sample_names, sample_barcode_maps):
    # for the main Snakefile as target output
    out = []
    for sample_i in sample_names:
        for cell_num_index in sample_barcode_maps[sample_i].keys():
            #"{sample_name}/sc_strelka2_call/{sample_name}_{cell_num_index}/results/variants/variants.vcf.gz",
            out.append(f"{sample_i}/sc_strelka2_call/{sample_i}_{cell_num_index}/results/variants/variants.vcf.gz")
    return out

###########
# ### STEP3 ----- single-cell Mutect2 filtering; merge single cell VCFs
###########
def get_step2_m2_filtered_vcfs_by_sample(wildcards):
    # for Snakemake rule use only
    # gets all step2 bcftools-filtered vcf for A GIVEN sample
    # as input for step3 bcftools merge
    out = []
    for cell_barcode in sample_barcode_maps[wildcards.sample_name]:
        out.append(f"{wildcards.sample_name}/3-bcf_filter/filtered_vcfs/{wildcards.sample_name}_{cell_barcode}_hard_filtered.vcf.gz")
    return out


###########
# ### STEP4 ----- single-cell genotying of SNVs in q_vcf from above
###########

# def get_step4_outputs(sample_names, sample_barcode_maps):
#     # for the main Snakefile as target output
#     out = []
#     for sample_i in sample_names:
#         for cell_barcode in sample_barcode_maps[sample_i]:
#             out.append(f"{sample_i}/mutect2_sc_f_pass2/m2_sc_vcfs_filter_added/{sample_i}_{cell_barcode}_somatic_m2_filter_added.vcf.gz")
#     return out

def get_step4_m2_f_vcfs(wildcards):
    out = []
    for cell_barcode in sample_barcode_maps[wildcards.sample_name]:
        out.append(f"{wildcards.sample_name}/4-sc_mutect_f_call/m2_sc_vcfs/{wildcards.sample_name}_{cell_barcode}_somatic_m2.vcf.gz")
    return out

def get_step4_mpileup_vcfs(wildcards):
    out = []
    for cell_num_index in sample_barcode_maps[wildcards.sample_name].keys():
        out.append(f"{wildcards.sample_name}/4-bcf_genotyping/sc_mpileup_vcfs/{wildcards.sample_name}_{cell_num_index}_raw_counts_added.vcf.gz")
    return out