# run in mosaic-custom
CASE_NAME=RA17_13
SAMPLE_NAMES=(RA17_13-34_1 RA17_13-44_1 RA17_13-50_1 RA17_13-55_1 RA17_13-56_1 RA17_13-67_1)

for SAMPLE_i in ${SAMPLE_NAMES[@]}; do
    echo processing: ${SAMPLE_i}
    # run orthogonal annotation
    # python /home/zhangh5/work/Tapestri_analysis/tap_analysis_SNAKEMAKE/scripts/STEP6-annotate_h5_with_bulk.py \
    #     --sample_name ${SAMPLE_i} \
    #     --input_h5 /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/${CASE_NAME}/${SAMPLE_i}/OUTPUTS_from_mpileup/${SAMPLE_i}.mpileup.h5 \
    #     --bq_info_csv /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/${CASE_NAME}/${SAMPLE_i}/3-bcf_filter/merged_bq_info/${SAMPLE_i}-bq_merged.sc_prev.csv \
    #     --bulk_info_yaml /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/${CASE_NAME}/${SAMPLE_i}/reference/${SAMPLE_i}_matched_bulk_info.yaml \
    #     --output_dir /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/${CASE_NAME}/${SAMPLE_i}/OUTPUTS_from_mpileup 

    # CRAVAT annotation
    python /home/zhangh5/work/Tapestri_analysis/tap_analysis_SNAKEMAKE/scripts/STEP6-query_snv_df_with_CRAVAT.py \
        --sample_name ${SAMPLE_i} \
        --input_snv_csv /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/${CASE_NAME}/${SAMPLE_i}/OUTPUTS_from_mpileup/${SAMPLE_i}.bulk_annotated_SNV.csv \
        --cravat_settings_yaml /home/zhangh5/work/Tapestri_analysis/tap_analysis_SNAKEMAKE/configs/cravat_settings.yaml \
        --output_dir /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/${CASE_NAME}/${SAMPLE_i}/OUTPUTS_from_mpileup/annotations \
        --log_file /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/${CASE_NAME}/${SAMPLE_i}/OUTPUTS_from_mpileup/annotations/${SAMPLE_i}.annotate_cravat.log
done