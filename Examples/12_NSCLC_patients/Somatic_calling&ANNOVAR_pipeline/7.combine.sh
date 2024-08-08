#!/bin/bash -e
mkdir -p SomSNVs
for j in `cat sample.list`;do
	i=`find 1* -name "$j*fpfilter.target.vcf"`
	patient=`echo $i|cut -d "/" -f1`
	name=`basename $i|cut -d "." -f1`
	bgzip -@ 10 $i
	bcftools index $i.gz

	bgzip -@ 10 $patient/6.fpfilter_LOH/$name.fpfilter.LOH.target.vcf
	bcftools index $patient/6.fpfilter_LOH/$name.fpfilter.LOH.target.vcf.gz

	bcftools concat -a $i.gz $patient/6.fpfilter_LOH/$name.fpfilter.LOH.target.vcf.gz -D > SomSNVs/$name.fpfilter.target.SomSNVs.vcf

done
