#!/bin/bash -e

for i in 16011 16031 17004 17005 17008 17011 17012 17017 17028 17029 17030 18001;do

	vcftools --vcf /GATK_orig_results/$i.recalibrate_variants.vcf --minDP 3 --max-missing 0.7 --keep sample.list --recode --stdout |awk '!/^#/{print $1"\t"$2"\t"$2}' > GATK_germline_0.7/$i.recalibrate_variants_0.7.bed
	       
done



