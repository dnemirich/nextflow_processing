#!bin/bash


less data/gencode.v31.chr_patch_hapl_scaff.annotation.gtf | grep exon | grep chr | awk '{print $1"\t"$4"\t"$5}' > data/exons.coords.bed

bcftools view  -m2 -M2  data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf  | bgzip -c \
    > data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.biallelic.vcf.gz 

bedtools intersect -header -wa -v -a data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.biallelic.vcf.gz \
    -b data/exons.coords.bed | bgzip -c > data/ALL.wgs.biallelic.exons.vcf.gz

bcftools view --min-af 0.02 data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.biallelic.vcf.gz | bgzip -c > \
     data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.biallelic.af.vcf.gz 

bedtools intersect -header -wa -v -a data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.biallelic.af.vcf.gz \
    -b data/exons.coords.bed | \
    bgzip -c > data/ALL.wgs.biallelic.exons.af.vcf.gz

tabix -p vcf data/ALL.wgs.biallelic.exons.af.vcf.gz


