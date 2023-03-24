rule h5_from_mpileup_raw_data_and_rc:
    # write SNV and CNV matrices.
    input:
        SC_MPILEUP_AD_LAYER_MERGED = "{sample_name}/{sample_name}.mpileup.AD.merged.csv",
        SC_MPILEUP_DP_LAYER_MERGED = "{sample_name}/{sample_name}.mpileup.DP.merged.csv",
        barcode_num_map_f =  "{sample_name}/reference/{sample_name}.barcode_map.txt",
        read_counts_tsv = "{sample_name}/tap_pipeline_output/results/tsv/{sample_name}.tube1.barcode.cell.distribution.tsv",
        panel_insert_file = config['reference_info']['panel_insert_file_UCSC'], # this need to be UCSC format ('chr1' instead of '1')
        panel_amplicons_file = config['reference_info']['panel_amplicon_file'],
    output:
        read_count_tsv_renamed = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}.per_amplicon_read_counts.tsv",
        output_h5 = "{sample_name}/OUTPUTS_from_mpileup/{sample_name}.mpileup.h5",
    params:
        output_dir = "{sample_name}/OUTPUTS_from_mpileup",
        write_h5_script = config['python_scripts']['WRITE_H5_FROM_MPILEUP_RAW_DATA_AND_RC'],
    log: 
        std = "{sample_name}/logs/{sample_name}.STEP5-write_h5_from_raw_data.log",
        err = "{sample_name}/logs/{sample_name}.STEP5-write_h5_from_raw_data.err",
    conda:
        "../envs/mosaic-custom.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = 59
    shell:
        """
        python {params.write_h5_script} \
            --sample_name {wildcards.sample_name} \
            --AD_merged_df {input.SC_MPILEUP_AD_LAYER_MERGED} \
            --DP_merged_df {input.SC_MPILEUP_DP_LAYER_MERGED} \
            --panel_insert_file {input.panel_insert_file} \
            --output_dir {params.output_dir} \
            1> {log.std} 2> {log.err}
        """