#!/bin/bash

for patient_number in  18001 17004 16011 16031 17005 17008 17011 17012 17017 17028 17029 17030

do 

{ 
  snp_dir=./$patient_number/snp/
  indel_dir=./$patient_number/indel/
  mkdir ./$patient_number/sF
  sF_dir=./$patient_number/sF

  for h in `ls $snp_dir|grep "snp.vcf"`

  do

  {
      g=`echo $h|cut -d "." -f 1`
      m=`ls $indel_dir|grep "$g"`
      varscan somaticFilter $snp_dir/$h  --min-coverage 10 --min-strands2 2 --min-var-freq 0.1 --p-value  0.05 --indel-file $indel_dir/$m --output-file $sF_dir/$g.sF.snv.vcf 
   }&

  done

wait

}

done
