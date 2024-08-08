#!/bin/bash -e

bam_directory=/data/picard_bam/

for i in `cat sample.list`;do
  patient_number=${i::5}
  
  sF_pS_dir=$patient_number/sF_pS/
  rmrp=$patient_number/sF_pS_rmrp
  rmrp_LOH=$patient_number/sF_pS_rmrp_LOH
  rmrp_Germline=$patient_number/sF_pS_rmrp_Germline

  bam-readcount-0.8.0/build/bin/bam-readcount -q 20 -b 20 -f /data/hs37d5/hs37d5.fa -w 1 -l $rmrp/$i.rmrp.bed $bam_directory/$i.markdup.realigned.bam > $patient_number/5.bamreadcount/$i.readcount
  bam-readcount-0.8.0/build/bin/bam-readcount -q 20 -b 20 -f /data/hs37d5/hs37d5.fa -w 1 -l $rmrp_LOH/$i.sF.snv.LOH.hc.vcf.rmrp.LOH.bed $bam_directory/$i.markdup.realigned.bam > $patient_number/5.bamreadcount_LOH/$i.LOH.readcount
  bam-readcount-0.8.0/build/bin/bam-readcount -q 20 -b 20 -f /data/hs37d5/hs37d5.fa -w 1 -l $rmrp_Germline/$i.sF.snv.Germline.hc.vcf.rmrp.Germline.bed $bam_directory/$i.markdup.realigned.bam > $patient_number/5.bamreadcount_Germline/$i.Germline.readcount

  varscan fpfilter $sF_pS_dir/hc_somatic/$i.*.Somatic.hc.vcf $patient_number/5.bamreadcount/$i.readcount --output-file $patient_number/6.fpfilter/$i.fpfilter.vcf --dream3-settings 1
  varscan fpfilter $sF_pS_dir/hc_LOH/$i.*.LOH.hc.vcf $patient_number/5.bamreadcount_LOH/$i.LOH.readcount --output-file $patient_number/6.fpfilter_LOH/$i.fpfilter.LOH.vcf --dream3-settings 1
  varscan fpfilter $sF_pS_dir/hc_germline/$i.*.Germline.hc.vcf $patient_number/5.bamreadcount_Germline/$i.Germline.readcount --output-file $patient_number/6.fpfilter_Germline/$i.fpfilter.Germline.vcf --dream3-settings 1

done









