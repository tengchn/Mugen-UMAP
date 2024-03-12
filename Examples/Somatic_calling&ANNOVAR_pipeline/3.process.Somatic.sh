#!/bin/bash

for patient_number in 17004 16011 16031 17005 17008 17011 17012 17017 17028 17029 17030 18001;do 

  sF_dir=./$patient_number/sF/
  mkdir ./$patient_number/sF_pS/
  sF_pS_dir=./$patient_number/sF_pS/
 
  for h in `ls $sF_dir|grep "snv.vcf"`

  do

  {
      varscan processSomatic $sF_dir/$h --min-tumor-freq 0.1 --max-normal-freq 0.02  --p-value 0.05 
   }&

  done

  wait

  mkdir $sF_pS_dir/hc_somatic $sF_pS_dir/hc_germline $sF_pS_dir/hc_LOH $sF_pS_dir/left

  mv $sF_dir/*.Somatic.hc.vcf $sF_pS_dir/hc_somatic

  mv $sF_dir/*.Germline.hc.vcf $sF_pS_dir/hc_germline 

  mv $sF_dir/*.LOH.hc.vcf $sF_pS_dir/hc_LOH

  mv $sF_dir/*LOH.vcf $sF_dir/*Somatic.vcf $sF_dir/*Germline.vcf $sF_pS_dir/left

done
