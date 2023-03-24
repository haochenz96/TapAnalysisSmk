from pathlib import Path

cohort_name = config['sample_info']['cohort_name']
sample_names = config['sample_info']['sample_names']
wd = '/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom'
wd = f"{wd}/{cohort_name}"
print(f"[INFO] wd set to: {wd}")

# snakemake workdir
workdir: wd

rule all:
    input:
        expand("{SAMPLE_i}/OUTPUTS_from_mpileup/{SAMPLE_i}.bulk_annotated_SNV.csv", SAMPLE_i=sample_names),
        expand("{SAMPLE_i}/OUTPUTS_from_mpileup/annotations/{SAMPLE_i}_CRAVAT_output.cleaned.txt", SAMPLE_i=sample_names),


rule annotate_with_bulk:
    input: 
        h5 = "{SAMPLE_i}/OUTPUTS_from_mpileup/{SAMPLE_i}.mpileup.h5",
        bq_info_csv = "{SAMPLE_i}/3-bcf_filter/merged_bq_info/{SAMPLE_i}-bq_merged.sc_prev.csv",
        bulk_info_yaml = "{SAMPLE_i}/reference/{SAMPLE_i}_matched_bulk_info.yaml",
    output: 
        annotated_snv_csv = "{SAMPLE_i}/OUTPUTS_from_mpileup/{SAMPLE_i}.bulk_annotated_SNV.csv",
    params:
        sample_name = "{SAMPLE_i}",
        # h5 = "{SAMPLE_i}/OUTPUTS_from_mpileup/{SAMPLE_i}.mpileup.h5",
        # bulk_info_yaml = "{SAMPLE_i}/reference/{SAMPLE_i}_matched_bulk_info.yaml",
        # output_annotated_snv_csv = "{SAMPLE_i}/OUTPUTS_from_mpileup/{SAMPLE_i}.bulk_annotated_SNV.csv"
    log:
        "{SAMPLE_i}/OUTPUTS_from_mpileup/annotations/{SAMPLE_i}.bulk-annotate.log"
    conda: 
    # hardcoded
        '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/envs/mosaic-custom.yaml'
    resources:
        mem_mb = 8000,
        n_threads = 4,
        time_min = 49,
    shell:
        """
        python /home/zhangh5/work/Tapestri_analysis/tap_analysis_SNAKEMAKE/scripts/STEP6-annotate_h5_with_bulk.py \
            --sample_name {params.sample_name} \
            --input_h5 {input.h5} \
            --bq_info_csv {input.bq_info_csv} \
            --bulk_info_yaml {input.bulk_info_yaml} \
            --output_annotated_snv_csv {output.annotated_snv_csv} \
            --log_file {log}
        """
rule annotate_with_CRAVAT:
    input:
        annotated_snv_csv = "{SAMPLE_i}/OUTPUTS_from_mpileup/{SAMPLE_i}.bulk_annotated_SNV.csv",
        cravat_settings_yaml = "/home/zhangh5/work/Tapestri_analysis/tap_analysis_SNAKEMAKE/configs/cravat_settings.yaml",
    output:
        cravat_input_df = protected("{SAMPLE_i}/OUTPUTS_from_mpileup/annotations/{SAMPLE_i}_CRAVAT_input.txt"),
        cravat_output_df = protected("{SAMPLE_i}/OUTPUTS_from_mpileup/annotations/{SAMPLE_i}_CRAVAT_output.txt"), 
        cravat_output_df_cleaned = "{SAMPLE_i}/OUTPUTS_from_mpileup/annotations/{SAMPLE_i}_CRAVAT_output.cleaned.txt",
    params:
        sample_name = "{SAMPLE_i}",
        output_dir = "{SAMPLE_i}/OUTPUTS_from_mpileup/annotations" # overrides in command line input for query Python script below
    conda: 
    # hardcoded
        '/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/workflow/envs/mosaic-custom.yaml'
    log: 
        "{SAMPLE_i}/OUTPUTS_from_mpileup/annotations/{SAMPLE_i}.CRAVAT-annotate.log"
    resources:
        mem_mb = 8000,
        n_threads = 4,
        time_min = 49,
    shell:
        """
        python /home/zhangh5/work/Tapestri_analysis/tap_analysis_SNAKEMAKE/scripts/STEP6-query_snv_df_with_CRAVAT.py \
            --sample_name {params.sample_name} \
            --input_snv_csv {input.annotated_snv_csv} \
            --cravat_settings_yaml {input.cravat_settings_yaml} \
            --cravat_output {output.cravat_output_df} \
            > {log} 2>&1
        """