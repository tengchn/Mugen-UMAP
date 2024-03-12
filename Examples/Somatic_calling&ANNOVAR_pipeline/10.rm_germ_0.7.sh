#!/bin/bash -e

for i in `ls SomSNVs_annovar/*.txt`;do
	name=`echo $(basename $i)|cut -d "." -f1`
	patient=${name::5}
	bedtools intersect -a $i -b GATK_germline_0.7/$patient.recalibrate_variants_0.7.bed -wa -u -header > ${i%%txt}germline_0.7.txt

done

