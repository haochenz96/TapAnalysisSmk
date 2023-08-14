from pathlib import Path
import mosaic.io as mio
import pandas as pd
from tea.format import isNaN

pipeline_output_dir = Path('/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom')
case_name = 'RA17_22'
sample_names = ['RA17_22-04_1', 'RA17_22-06_2', 'RA17_22-11_1', 'RA17_22-27_2', 'RA17_22-32_1', 'RA17_22-33_2', 'RA17_22-35_1', 'RA17_22-39_6', 'RA17_22-42_2']

output_dir = pipeline_output_dir / case_name / 'cw_genotype'
output_dir.mkdir(exist_ok=True, parents=True)
voi_sets = []

for sample_i in sample_names:

    # Get H5 and CRAVAT files
    h5 = pipeline_output_dir / case_name / sample_i / 'OUTPUTS_from_mpileup' / f"{sample_i}.mpileup.h5"
    cravat_f = pipeline_output_dir / case_name / sample_i / 'OUTPUTS_from_mpileup' / 'annotations' / f"{sample_i}_CRAVAT_output_cleaned.txt"

    sample_obj = mio.load(h5)
    cravat_df = pd.read_csv(cravat_f, sep='\t', index_col=0, header=[0,1])

    bulk_mask =  (cravat_df[('bulk_comparison', 'bulk-matched_bulk_normal-AF')] > 0) | (cravat_df[('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] > 0)

    voi = cravat_df[bulk_mask].index.tolist()
    voi_sets.append(set(voi))

# Get intersection of all sets
voi_intersect = set.intersection(*voi_sets)

# write to output
with open(output_dir / f"{case_name}-bulk_snvs.txt", 'w') as f:
    for snv in voi_intersect:
        f.write(f"{snv}\n")