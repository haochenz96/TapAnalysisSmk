if Path(config['input']['Q_alleles_txt']).exists():
    print(f"[INFO] -- candidate alleles TXT file exists: {config['input']['Q_alleles_txt']}. It will be directly used for mpileup+call. Please make sure it has multiallelic sites merged.")

else:
    rule generate_candidate_alleles_for_bcf:
        # generate a candidate allele from Q_VCF
        # @HZ 08/28/2022: need to merge multiallelic sites for the next step (mpileup+call) to work
        input:
            Q_VCF = config['input']['Q_VCF'],
        output: 
            CANDIDATE_ALLELE = config['input']['Q_alleles_txt'],
            CANDIDATE_ALLELE_for_py = "input/candidate_alleles.multiallelic.for_py.csv",   
        threads: lambda wildcards, attempt: attempt * 2,
        resources: 
            mem_mb = lambda wildcards, attempt: attempt * 2000,
            time_min = lambda wildcards, attempt: attempt * 59,
        conda:
            config['envs']['bcftools'],
        params:
            Q_VCF_atomized = "input/candidate_alleles.atomized.vcf.gz",
            Q_VCF_multiallelic = "input/candidate_alleles.multiallelic.vcf.gz",
            REF_GENOME = config['reference_info']['reference_genome'],
            PANEL_BED = config['reference_info']['panel_insert_file'],
            BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
        log: 
            "logs/mpileup_genotype.log"
        shell:
            """
            TIMESTAMP=[`date`]
            echo $TIMESTAMP >> {log} 
            echo '--- generating candidate alleles from filtered, combined, normed Mutect2 VCF.' >> {log} 

            # @HZ 03/21/2023 atomize MNPs
            {params.BCFTOOLS_EXEC} norm \
                -Oz \
                -o {params.Q_VCF_atomized} \
                --atomize \
                {input.Q_VCF}

            # @HZ 09/18/2022 get a condensed form for Python use later
            {params.BCFTOOLS_EXEC} query \
                -f '%CHROM:%POS:%REF/%ALT\n' \
                {params.Q_VCF_atomized} > \
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

# rule generate_candidate_alleles_for_py:
#     input:
#         CANDIDATE_ALLELE = config['input']['Q_alleles_txt'],
#     output:
#         CANDIDATE_ALLELE_for_py = "input/candidate_alleles.multiallelic.for_py.csv",        
#     run:
#         import pandas as pd
#         alleles_df = pd.read_csv(input.CANDIDATE_ALLELE, sep='\t', names = ['CHR', 'POS', 'ALLELES'])
#         alleles_df['REF'] = alleles_df['ALLELES'].str.split(',').str[0]
#         alleles_df['ALT'] = alleles_df['ALLELES'].str.split(',').str[1]
#         output_df = pd.DataFrame(index = alleles_df['CHR'].astype(str) + ':' + alleles_df['POS'].astype(str) + ':' + alleles_df['REF'] + '/' + alleles_df['ALT'])
#         if not output_df.index.str.startswith('chr').any():
#             output_df.index = 'chr' + output_df.index
#         output_df.to_csv(output.CANDIDATE_ALLELE_for_py, header=False, index=True)

rule sc_mpileup:
    # scattered by single cell
    # mpileup and call variants 
    input:
        SC_BAM = get_sc_bam, # retrieve EACH sc_bam given the wildcards: sample_name and cell_num_index
        CANDIDATE_ALLELE = config['input']['Q_alleles_txt'],
        # HEADER_FILE = "input/candidate_alleles.header.txt",
    output:
        # SC_RAW_COUNTS_TXT = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_raw_counts.txt.gz",
        SC_MPILEUP_VCF = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.normed.vcf.gz",
        SC_MPILEUP_VCF_TBI = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.normed.vcf.gz.tbi",
        # SC_MPILEUP_VCF_raw_counts = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_raw_counts_added.vcf.gz",
        SC_MPILEUP_AD_LAYER = "{sample_name}/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.AD.csv",
        SC_MPILEUP_DP_LAYER = "{sample_name}/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.DP.csv",
    params:
        REF_GENOME = config['reference_info']['reference_genome'],
        # PANEL_BED = config['reference_info']['panel_insert_file'],
        BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
    threads: 2
    group: "sc_mpileup"
    resources: 
        mem_mb = lambda wildcards, attempt: min(8000, attempt * 4000),
        time_min = lambda wildcards, attempt: attempt * 119,
    conda:
        config['envs']['bcftools'],
    shell:
        """
        # --- 1. mpileup ---
        # @HZ 09/19/2022: need to add `-f [REF_GENOME]` to the norm command for indels to be normalized
        # otherwise the merging step won't consider them as they are not present in the alleles file

        # @HZ 03/08/2023: added sort step because the candidate alleles file could be unsorted
        {params.BCFTOOLS_EXEC} mpileup \
            {input.SC_BAM} \
            -R {input.CANDIDATE_ALLELE} \
            -f {params.REF_GENOME} \
            --annotate FORMAT/AD,FORMAT/DP,INFO/AD \
            --max-depth 100000 \
            --max-idepth 100000 \
            -Ou | \
        {params.BCFTOOLS_EXEC} norm \
        	-m- \
            -f {params.REF_GENOME} \
        	-Ou | \
        {params.BCFTOOLS_EXEC} sort \
            -Oz \
            > {output.SC_MPILEUP_VCF} && \
        tabix {output.SC_MPILEUP_VCF}

        # --- 2. extract AD, DP ---
        {params.BCFTOOLS_EXEC} query \
        	{output.SC_MPILEUP_VCF} \
        	-f '%CHROM:%POS:%REF/%ALT,%AD{{1}}\n' \
        	-i'ALT!="<*>"' \
        	> {output.SC_MPILEUP_AD_LAYER}
        {params.BCFTOOLS_EXEC} query \
            {output.SC_MPILEUP_VCF} \
            -f '%CHROM:%POS:%REF,%DP\n' \
            > {output.SC_MPILEUP_DP_LAYER}
        """

# rule extract_sc_mpileup_raw:
#     input:
#         SC_MPILEUP_VCF = "{sample_name}/sc_mpileup_vcfs/{sample_name}_{cell_num_index}_mpileup.normed.vcf.gz",
#     output:
#         SC_MPILEUP_DP_LAYER = "{sample_name}/sc_mpileup_raw_data/{sample_name}_{cell_num_index}_mpileup.normed.DP.csv",
#     params:
#         BCFTOOLS_EXEC = config['bcftools']['BCFTOOLS_EXEC'],
#     shell:
#         """
#         {params.BCFTOOLS_EXEC} query \
#             {input.SC_MPILEUP_VCF} \
#             -f '%CHROM:%POS:%REF,%DP\n' \
#             > {output.SC_MPILEUP_DP_LAYER}
#         """

rule merge_raw_layers:
    # scatter by sample
    # merge raw counts layers
    input:
        SC_MPILEUP_AD_LAYER = get_sc_mpileup_AD,
        SC_MPILEUP_DP_LAYER = get_sc_mpileup_DP,
        CANDIDATE_ALLELE_for_py = "input/candidate_alleles.multiallelic.for_py.csv",
    output:
        SC_MPILEUP_AD_LAYER_MERGED = "{sample_name}/{sample_name}.mpileup.AD.merged.csv",
        SC_MPILEUP_DP_LAYER_MERGED = "{sample_name}/{sample_name}.mpileup.DP.merged.csv",
    params:
        input_dir = "{sample_name}/sc_mpileup_raw_data",
        output_dir = "{sample_name}",
        MERGE_SCRIPT = config['scripts']['MERGE_AD_DP_MATRIX'],
    log:
        std = "{sample_name}/logs/mpileup_genotype.log",
        err = "{sample_name}/logs/mpileup_genotype.err",
    threads: lambda wildcards, attempt: attempt * 4,
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time_min = lambda wildcards, attempt: attempt * 59,
    shell:
        """
        python {params.MERGE_SCRIPT} \
            --sample_name {wildcards.sample_name} \
            --input_dir {params.input_dir} \
            --candidate_alleles_df {input.CANDIDATE_ALLELE_for_py} \
            --output_dir {params.output_dir} \
            1> {log.std} 2> {log.err}
        """