#!/bin/bash
for j in `cat sample.list`;do
	for i in `find 1* -name "$j*fpfilter*.vcf"`;do

	vcftools --vcf $i --bed /data/exome_bed/xgen_exome_target.bed --recode --stdout > ${i%%vcf}target.vcf

	done
done
