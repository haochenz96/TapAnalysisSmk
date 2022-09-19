from pathlib import Path
import sys
import glob
import re

# ----- fetch variables from config -----
sample_names = config['input']['sample_sc_bams'].keys()
cohort_name = config['input']['cohort_name']
sc_bams_dir_map = config['input']['sample_sc_bams']

workdir: Path(config['run_info']['work_dir']) / cohort_name

# Get each sample's cell indexes from glob'ing the input scBAMs dir
sample_to_sc_idx_map = {}
for sample_i in sample_names:
    sample_to_sc_idx_map[sample_i] = [
        re.findall('cell_\d+', x)[0] for x in glob.glob(f"{sc_bams_dir_map[sample_i]}/*.bam")
        ]
    # print(sample_to_sc_idx_map[sample_i])
    print(f"for {sample_i}, {len(sample_to_sc_idx_map[sample_i])} single cells' BAM files found!")

# ----- utility functions -----
def get_sc_bam(wildcards):
    return f"{sc_bams_dir_map[wildcards.sample_name]}/{wildcards.sample_name}_{wildcards.cell_num_index}.bam"

def get_sc_mpileup_vcfs(wildcards):
    out = []
    # cell_num_index is inferred here from the map
    for cell_num_index in sample_to_sc_idx_map[wildcards.sample_name]:
    #"{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}/{sample_name}_{cell_num_index}_raw_counts.vcf.gz",
        out.append(f"{wildcards.sample_name}/sc_mpileup_vcfs/{wildcards.sample_name}_{cell_num_index}_raw_counts_added.vcf.gz")
    print(f"for {wildcards.sample_name}, {len(out)} single cells' VCFs found!")
    # for i in range(30):
    #     print(out[i])
    return out

# ===== target ======
rule all:
    input:
        merged_mpileup_vcf = expand(
            "{sample_name}/combined_vcf/{sample_name}-genotyped_combined.vcf.gz",
            sample_name = sample_names,
        ),
        merged_mpileup_vcf_AD_py = expand(
            "{sample_name}/{sample_name}-genotyped_combined_AD_for_py.txt",
            sample_name = sample_names
        ),
# ===================

include: "single-cell_mpileup.smk"
include: "format_combined_vcf.smk"
