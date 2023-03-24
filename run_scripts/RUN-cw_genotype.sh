#!/bin/bash 
#BSUB -J RA17_22-cw_genotype-more_samples
#BSUB -n 1                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=16]    # expected resorce consumption for memory
#BSUB -W 5:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/Tapestri_batch2/analysis/cohort-wide_genotype/RA17_22/RA17_22-cohort_genotype-more_samples.stdout  # standard output directory
#BSUB -eo /home/zhangh5/work/Tapestri_batch2/analysis/cohort-wide_genotype/RA17_22/RA17_22-cohort_genotype-more_samples.stderr # standard error output file

# if [ -f ~/.bashrc ] ; then
#     . ~/.bashrc
# fi

# conda activate snakemake 
# # # # export SNAKEMAKE_OUTPUT_CACHE=/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/snakemake-cache/ # <--- define cache location
# # # export XDG_CACHE_HOME=/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/snakemake-cache/ # <--- define cache location

# snakemake \
#     -s /home/zhangh5/work/Tapestri_analysis/tap_analysis_SNAKEMAKE/rules/main.smk \
#     --configfile /home/zhangh5/work/Tapestri_batch2/analysis/cohort-wide_genotype/RA17_22/for_SCARLET/input/RA17_22-cw_genotype.yaml \
#     --conda-prefix /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/conda \
#     --profile lsf \
#     --quiet rules \
#     --rerun-incomplete 

conda activate mosaic-custom

CASE_NAME=RA17_22
SAMPLE_NAMES=(RA17_22-32_1 RA17_22-33_2 RA17_22-35_1 RA17_22-39_6 RA17_22-42_2 RA17_22-04_1 RA17_22-06_2 RA17_22-11_1 RA17_22-27_2)

ADD_AD_DP_TO_H5_SCRIPT=/home/zhangh5/work/Tapestri_project/TapVarCallSmk/workflow/scripts/STEP5-write_h5_from_raw_AD_DP_and_rc.py

for sample_i in ${SAMPLE_NAMES[@]}; do
    echo "Processing sample: ${sample_i}"
    python ${ADD_AD_DP_TO_H5_SCRIPT} \
        --sample_name ${sample_i} \
        --AD_merged_df /home/zhangh5/work/Tapestri_batch2/analysis/cohort-wide_genotype/${CASE_NAME}/for_SCARLET/output-03-23-23/${CASE_NAME}/${sample_i}/${sample_i}.mpileup.AD.merged.csv \
        --DP_merged_df /home/zhangh5/work/Tapestri_batch2/analysis/cohort-wide_genotype/${CASE_NAME}/for_SCARLET/output-03-23-23/${CASE_NAME}/${sample_i}/${sample_i}.mpileup.DP.merged.csv \
        --panel_insert_file /home/zhangh5/work/Tapestri_batch2/panel_3359/panel_3359_hg19-ucsc/3359.bed \
        --output_dir /home/zhangh5/work/Tapestri_batch2/analysis/cohort-wide_genotype/RA17_22/for_SCARLET/output-03-23-23/cw_genotype-H5
done