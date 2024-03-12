#!/bin/bash -e
ref=/data/hs37d5/hs37d5.fa
wk_dir=/data/picard_bam
for patient_number in 17017 17029 17005 17008 17012 17017 17029 16011 16031 18001 17011 17028;do

mkdir $patient_number
cd $patient_number

normal_bam=`ls $wk_dir/$patient_number/normal|grep "bam$"`

    for tumor_bam in `ls $wk_dir/$patient_number/|grep "bam$"`;do
     
      g=`echo $tumor_bam|cut -d "." -f 1`

      samtools mpileup -f $ref -q 20 -Q 20 -B $wk_dir/$patient_number/normal/$normal_bam $wk_dir/$patient_number/$tumor_bam| varscan somatic /dev/stdin $patient_number/$g.varscan.total.output --mpileup 1 --min-coverage 10 --min-coverage-normal 10 --min-coverage-tumor 10 --min-var-freq 0.1 --min-freq-for-hom 0.85 --normal-purity 1.0 --tumor-purity 1.00 --p-value 0.99 --somatic-p-value 0.05 --strand-filter 0 --output-vcf 1 &
      
    done

mkdir snp && mv *.snp.* snp
mkdir indel && mv *.indel.* indel

cd /data/varscan_result 

done 
