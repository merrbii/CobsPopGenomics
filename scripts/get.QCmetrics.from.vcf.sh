#!/bin/bash
# This script will extract important metrics that will help you define your cutoffs.
# By Errbii M (added 18.08.21)

#get all variables you need..

echo "Type the name of the vcf file you are working on now, followed by [ENTER]:"

read vcf

echo "Type the window size here, followed by [ENTER]:"
read windsize

echo "Type the window step here, followed by [ENTER]:"
read windstep


base=$(basename $vcf .vcf.gz)

echo "$base"

mkdir $(echo "$base"_QC)

out=$(echo "$base"_QC)

#no_snps=$(bcftools view $vcf|grep -v "#"|wc -l)

#echo "This file contains $no_snps"

echo "calculating allele frequency and count..."

#freq
$(vcftools --gzvcf $vcf --freq2 --out ./$out/$base.freq)

#count
$(vcftools --gzvcf $vcf --counts2 --out ./$out/$base.count)


echo "calculating site mean depth..."

$(vcftools --gzvcf $vcf --site-mean-depth --out ./$out/$base.site.depth)

echo "calculating the missingness on a per-individual basis...."

$(vcftools --gzvcf $vcf --missing-indv --out ./$out/$base.missing.indv)


echo "calculating the missingness on a per-site basis..."

$(vcftools --gzvcf $vcf --missing-site --out ./$out/$base.missing.site)


echo "calculating the number and density of SNPs in bins of 10kb..."

$(vcftools --gzvcf $vcf --SNPdensity 10000 --out ./$out/$base.density)

echo "done with calculating 5 metrics! Now got to R and visualize the data"
