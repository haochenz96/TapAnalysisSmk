# ----- rules -----
rule generate_candidate_allele:
    # generate a candidate allele from Q_VCF
    # @HZ 08/28/2022: need to merge multiallelic sites for the next step (mpileup+call) to work
    input:
        Q_VCF = config['input']['Q_VCF'],
    output: 
        CANDIDATE_ALLELE = config['input']['Q_alleles_txt'],
        CANDIDATE_ALLELE_for_py = config['input']['Q_alleles_txt_for_py'],
    threads: lambda wildcards, attempt: attempt * 2,
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 59,
    conda:
        config['envs']['bcftools']
    params:
        Q_VCF_multiallelic = "input/candidate_alleles.multiallelic.vcf.gz",
        REF_GENOME = config['reference_info']['reference_genome'],
        PANEL_BED = config['reference_info']['panel_insert_file'],
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    log: 
        "logs/bcf_preprocess_input.log"
    shell:
        """
        TIMESTAMP=[`date`]
        echo $TIMESTAMP >> {log} 
        echo '--- generating candidate alleles from filtered, combined, normed Mutect2 VCF.' >> {log} 

        # @HZ 09/18/2022 get a condensed form for Python use later
        {params.BCFTOOLS_EXEC} query \
            -f '%CHROM:%POS:%REF/%ALT\n' \
            {input.Q_VCF} > \
            {output.CANDIDATE_ALLELE_for_py}

        # @HZ 08/28/2022: merge multiallelic sites 
        {params.BCFTOOLS_EXEC} norm \
            -Oz \
            -o {params.Q_VCF_multiallelic} \
            --multiallelics + \
            {input.Q_VCF}

        {params.BCFTOOLS_EXEC} query \
            -f'%CHROM\t%POS\t%REF,%ALT\n' \
            {params.Q_VCF_multiallelic} | \
        bgzip -c > {output.CANDIDATE_ALLELE} && \
        tabix -s1 -b2 -e2 {output.CANDIDATE_ALLELE} && \
        echo '--- finished writing candidate alleles file.' >> {log}
        
        echo '--- next, bcftools mpileup with:' >> {log} && \
        echo '--- --- BED FILE: {params.PANEL_BED}' >> {log} && \
        echo '--- --- GENOME FILE: {params.REF_GENOME}' >> {log} 
        """

rule prepare_header:
    # prepara a header file for next step
    input:
        CANDIDATE_ALLELE = config['input']['Q_alleles_txt'],
    output:
        HEADER_FILE = "input/candidate_alleles.header.txt", # this is for the next step
    shell:
        """
        # write header for next step
        echo -e '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Total allelic depths (high-quality bases)">' >> {output.HEADER_FILE}
        """

rule sc_mpileup:
    # scattered by single cell
    # mpileup and call variants 
    input:
        SC_BAM = get_sc_bam, # retrieve each sc_bam given the wildcards: sample_name and cell_num_index
        CANDIDATE_ALLELE = config['input']['Q_alleles_txt'],
        HEADER_FILE = "input/candidate_alleles.header.txt",
    output:
        SC_RAW_COUNTS_TXT = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_raw_counts.txt.gz",
        SC_MPILEUP_VCF = temp("{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.vcf.gz"),
        SC_MPILEUP_VCF_raw_counts = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_raw_counts_added.vcf.gz",
    params:
        REF_GENOME = config['reference_info']['reference_genome'],
        PANEL_BED = config['reference_info']['panel_insert_file'],
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    threads: 2
    resources: 
        mem_mb = lambda wildcards, attempt: min(8000, attempt * 4000),
        time_min = lambda wildcards, attempt: attempt * 119,
    conda:
        config['envs']['bcftools']
    shell:
        """
        {params.BCFTOOLS_EXEC} mpileup \
            -Ou \
            -R {input.CANDIDATE_ALLELE} \
            -f {params.REF_GENOME} \
            --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
            --max-depth 100000 \
            --max-idepth 100000 \
            {input.SC_BAM} | \
        {params.BCFTOOLS_EXEC} call \
            --keep-alts \
            -C alleles \
            -T {input.CANDIDATE_ALLELE} \
            --multiallelic-caller \
            -Ou | \
        {params.BCFTOOLS_EXEC} norm \
            -Ou \
            --multiallelics - | \
        {params.BCFTOOLS_EXEC} sort \
            -Oz \
            -o {output.SC_MPILEUP_VCF} && \
        tabix {output.SC_MPILEUP_VCF}

        {params.BCFTOOLS_EXEC} query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t[ %AD]\n' \
            {output.SC_MPILEUP_VCF} | \
        bgzip -c > {output.SC_RAW_COUNTS_TXT} && \
        tabix -s1 -b2 -e2 {output.SC_RAW_COUNTS_TXT}

        SAMPLE_NAME=$({params.BCFTOOLS_EXEC} query -l {output.SC_MPILEUP_VCF})

        {params.BCFTOOLS_EXEC} annotate \
            -s $SAMPLE_NAME \
            -a {output.SC_RAW_COUNTS_TXT} \
            -h {input.HEADER_FILE} \
            -Oz \
            -o {output.SC_MPILEUP_VCF_raw_counts} \
            -c CHROM,POS,REF,ALT,FORMAT/AD \
            {output.SC_MPILEUP_VCF} && \
        tabix {output.SC_MPILEUP_VCF_raw_counts}

        # clean up:
        rm {output.SC_RAW_COUNTS_TXT}.tbi
        rm {output.SC_MPILEUP_VCF}.tbi

        """

rule norm_and_merge_sc_mpileup:
    # scatter by sample
    # (1) for each single cell VCF, split multiallelic sites and index (OBSOLETE WITH RAW COUNTS EXTRACTION); 
    # (2) merge filtered single cell VCFs; no filter based on mutational prevalence across all cells; create statistics.
    input:
        sc_mpileup_vcfs = get_sc_mpileup_vcfs,
    output:
        merged_mpileup_vcf = "{sample_name}/combined_vcf/{sample_name}-genotyped_combined.vcf.gz",
        merged_mpileup_vcf_stats = "{sample_name}/combined_vcf/{sample_name}-genotyped_combined.vcf.gz.stats",
    params:
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    log: 
        "logs/{sample_name}.bcf_merge.log"
    conda:
        config['envs']['bcftools']
    threads: lambda wildcards, attempt: attempt * 4,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 119,
    shell:
        """
        {params.BCFTOOLS_EXEC} merge \
            --merge none \
            {input.sc_mpileup_vcfs} | \
        {params.BCFTOOLS_EXEC} sort | \
        bgzip -c > {output.merged_mpileup_vcf} && \
        tabix {output.merged_mpileup_vcf}

        {params.BCFTOOLS_EXEC} stats {output.merged_mpileup_vcf} > {output.merged_mpileup_vcf_stats}

        """

