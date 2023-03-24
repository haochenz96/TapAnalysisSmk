INPUT_DIR=/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas
OUT_DIR=/home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/normal_pancreas/RA18_18-11_1/4-bcf_genotyping/test
BCFTOOLS_EXEC=/home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/resources/bcftools-1.15.1/bcftools


cell_idx=(1011 1 212)

for i in ${cell_idx[@]}
do
	echo "processing cell: ${i}"
	# # --- 1. mpileup ---
	# ${BCFTOOLS_EXEC} mpileup \
	# 	-R ${INPUT_DIR}/RA18_18-11_1/4-bcf_genotyping/RA18_18-11_1-candidate_alleles.multiallelic.tsv.gz \
	# 	-f /home/zhangh5/work/commons/Tapestri_ref/hg19-b37/hg19-b37.fasta \
	# 	--annotate FORMAT/AD,FORMAT/DP,INFO/AD \
	# 	--max-depth 100000 \
	# 	--max-idepth 100000 \
	# 	${INPUT_DIR}/RA18_18-11_1/1-sc_bams/RA18_18-11_1_cell_${i}.bam \
	# 	-Ou | \
	# ${BCFTOOLS_EXEC} norm \
	# 	-m- \
	# 	-f /home/zhangh5/work/commons/Tapestri_ref/hg19-b37/hg19-b37.fasta \
	# 	-Oz > \
	# 	${OUT_DIR}/RA18_18-11_1_cell_${i}.mpileup.normed.vcf.gz
	
	# --- 2. extract AD, DP ---
	${BCFTOOLS_EXEC} query \
		${OUT_DIR}/RA18_18-11_1_cell_${i}.mpileup.normed.vcf.gz \
		-f '%CHROM:%POS:%REF/%ALT\t%AD{1}\n' \
		-i'ALT!="<*>"' \
		> ${OUT_DIR}/RA18_18-11_1_cell_${i}.mpileup.normed.AD.txt
	${BCFTOOLS_EXEC} query \
		${OUT_DIR}/RA18_18-11_1_cell_${i}.mpileup.normed.vcf.gz \
		-f '%CHROM:%POS\t%DP\n' \
		> ${OUT_DIR}/RA18_18-11_1_cell_${i}.mpileup.normed.DP.txt

done
# /home/zhangh5/work/Tapestri_project/tap_custom_SNAKEMAKE/resources/bcftools-1.15.1/bcftools call \
# 	-C alleles \
# 	-T ${INPUT_DIR}/RA18_18-11_1/4-bcf_genotyping/RA18_18-11_1-candidate_alleles.tsv.gz \
# 	--multiallelic-caller \
# 	-Oz > \
#     ${OUT_DIR}/RA18_18-11_4_cell_212.genotyped.no_alts.vcf.gz