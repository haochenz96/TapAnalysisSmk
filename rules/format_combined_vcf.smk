rule get_AD_from_combined_vcf:
    input:
        merged_mpileup_vcf = "{sample_name}/combined_vcf/{sample_name}-genotyped_combined.vcf.gz",
    output:
        merged_mpileup_vcf_AD = temp("{sample_name}/combined_vcf/{sample_name}-genotyped_combined_AD.txt"),
    log: "logs/{sample_name}.bcf_merge.log",
    conda:
        config['envs']['bcftools']
    shell:
        """
        bcftools query \
            -H \
            -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
            {input.merged_mpileup_vcf} > {output.merged_mpileup_vcf_AD}
        """
rule format_AD_for_py:
    input:
        merged_mpileup_vcf_AD = "{sample_name}/combined_vcf/{sample_name}-genotyped_combined_AD.txt",
    output:
        merged_mpileup_vcf_AD_py = "{sample_name}/{sample_name}-genotyped_combined_AD_for_py.txt",
    params:
        FORMAT_SCRIPT=config['scripts']['FORMAT_AD_MATRIX_FOR_PY'],
    log: 
        std = "logs/{sample_name}.bcf_merge.log",
        err = "logs/{sample_name}.bcf_merge.err",
    # run in native Snakemake environment, which should have pandas installed
    shell:
        """
        python {params.FORMAT_SCRIPT} \
            --sample_name {wildcards.sample_name} \
            --AD_TXT {input.merged_mpileup_vcf_AD} \
            --AD_TXT_FORMATTED_FOR_PY {output.merged_mpileup_vcf_AD_py} \
            1> {log.std} 2> {log.err}
        """
