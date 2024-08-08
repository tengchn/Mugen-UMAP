#!/bin/bash -e

for j in `cat sample.list`;do
	i=SomSNVs/$j.*.vcf
	/tools/annovar/table_annovar.pl $i /tools/annovar/humandb/ -buildver hg19 -out SomSNVs_annovar/$j -remove -protocol refGene,cytoBand,dbnsfp35c,clinvar_20200316,cosmic97_coding,avsnp150,gnomad211_exome,exac03,intervar_20180118,revel -operation g,r,f,f,f,f,f,f,f,f -nastring . -vcfinput --thread 16 2>&1 | tee SomSNVs_annovar/$j.annovar.log

done





