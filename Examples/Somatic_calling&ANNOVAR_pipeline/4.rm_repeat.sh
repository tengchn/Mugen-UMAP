#!/bin/bash -e

repeat_region=/data/list/repeat.blacklist.bed

for i in `cat sample.list`;do
   patient_number=${i::5}
   sF_pS_dir=$patient_number/sF_pS
   mkdir ./$patient_number/sF_pS_bedfile/ ./$patient_number/sF_pS_rmrp
   sF_pS_bed_dir=$patient_number/sF_pS_bedfile
   rmrp=$patient_number/sF_pS_rmrp
   rmrp_LOH=$patient_number/sF_pS_rmrp_LOH
   rmrp_Germline=$patient_number/sF_pS_rmrp_Germline
   mkdir -p $rmrp

    awk '!/^#/{print $1"\t"$2"\t"$2}' $sF_pS_dir/hc_somatic/$i.sF.snv.Somatic.hc.vcf > $sF_pS_bed_dir/$i.sF_pS.bed
    awk '!/^#/{print $1"\t"$2"\t"$2}' $sF_pS_dir/hc_LOH/$i.sF.snv.LOH.hc.vcf > $sF_pS_bed_dir/$i.sF.snv.LOH.hc.vcf.LOH.bed
    awk '!/^#/{print $1"\t"$2"\t"$2}' $sF_pS_dir/hc_germline/$i.sF.snv.Germline.hc.vcf > $sF_pS_bed_dir/$i.sF.snv.Germline.hc.vcf.Germline.bed

    bedtools intersect -v -a $sF_pS_bed_dir/$i.sF_pS.bed -b $repeat_region -wa > $rmrp/$i.rmrp.bed
    bedtools intersect -v -a $sF_pS_bed_dir/$i.sF.snv.LOH.hc.vcf.LOH.bed -b $repeat_region -wa > $rmrp_LOH/$i.sF.snv.LOH.hc.vcf.rmrp.LOH.bed
    bedtools intersect -v -a $sF_pS_bed_dir/$i.sF.snv.Germline.hc.vcf.Germline.bed -b $repeat_region -wa > $rmrp_Germline/$i.sF.snv.Germline.hc.vcf.rmrp.Germline.bed
    
done
