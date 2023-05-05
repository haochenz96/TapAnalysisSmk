#!/bin/bash 
#BSUB -J SNV-CNV-combined_analysis
#BSUB -n 4                # number of core(tasks) for parallel jobs          
#BSUB -sla IACOBUZC
#BSUB -R rusage[mem=36]    # expected resorce consumption for memory
#BSUB -W 4:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/SNV-CNV-combined_analysis.stdout
#BSUB -eo /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/SNV-CNV-combined_analysis.stderr

if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate mosaic-custom

# COHORT_NAME=(RA15_06 RA15_16 RA17_13 M04 M07 M13 TP6 RA21_17) # RA16_08 RA16_29 
# COHORT_NAME=RA15_06
# COHORT_NAME=(M04 M07 M13)
# COHORT_NAME=RA15_16
COHORT_NAME=RA17_13
# # SAMPLE_NAME=RA17_22-04_1
# COHORT_NAME=(M11 M12)
# COHORT_NAME=RA21_17
# COHORT_NAME=TP12
# COHORT_NAME=(M07) # M11 M12 M13)


for COHORT_NAME_i in ${COHORT_NAME[@]}
do
    echo "Processing ${COHORT_NAME_i} ..."
    python /home/zhangh5/work/Tapestri_analysis/TapAnalysisSmk/scripts/make_sc_snv_heatmaps.py \
        --wd /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/${COHORT_NAME_i} \
        --analysis_config /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/${COHORT_NAME_i}/${COHORT_NAME_i}-analysis_config.yaml \
        2>&1 | tee /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/${COHORT_NAME_i}/${COHORT_NAME_i}-make_snv_cnv_heatmaps.log

    # python /home/zhangh5/work/Tapestri_analysis/TapAnalysisSmk/scripts/make_snv_mut_prev_histogram.py \
    #     --wd /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/${COHORT_NAME_i} \
    #     --analysis_config /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/${COHORT_NAME_i}/${COHORT_NAME_i}-analysis_config.yaml \
    #     2>&1 | tee /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/${COHORT_NAME_i}/${COHORT_NAME_i}-make_snv_mut_prev_histogram.log
done