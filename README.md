# CobsPopGenomics
# Bioinformatic pipeline from analyzing genomic data of _Cardiocondyla obscurior_
##### By Mohammed Errbii "Simo"

This markdown document gives a detailed description of the bioinformatic steps followed to analyze WGS of pools and single individuals of _Cardiocondyla obscurior_. Please have a look at our [paper](https://onlinelibrary.wiley.com/doi/10.1111/mec.16099) for further details on our research hypotheses and major findings.

### I- Mapping, variant calling and filtering
#### 1) Filter raw reads for BQ > 20 and minimum length > 40bp by using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```bash

java -jar trimmomatic-0.38.jar PE \
reads_R1.fastq.gz \
reads_R2.fastq.gz \
reads_R1_paired.fastq.gz \
reads_R1_unpaired.fastq.gz \
reads_R2_paired.fastq.gz \
reads_R2_unpaired.fastq.gz \
SLIDINGWINDOW:4:20 \
MINLEN:40
```

#### 2) map trimmed reads with [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)

```bash

bwa mem \
-t 16 \
Cobs2.1.reference.fa \
reads_R1_paired.fastq.gz \
reads_R2_paired.fastq.gz > alignment.sam
```

#### 3) pre-processing of the alignments using [Picard](https://broadinstitute.github.io/picard/) and [SAMtools](http://samtools.sourceforge.net/)

```bash
# a- clean sam files
java -jar picard.jar CleanSam I=alignment.sam O=alignment.C.sam

# b- covert to bam format
samtools view -S -b alignment.C.sam > alignment.C.bam

# c- sort bam file (sort by default assume sorting by coordinates):
samtools sort alignment.C.bam > alignment.CS.bam

# d- fix mate information if needed
java -jar picard.jar FixMateInformation \
I=alignment.CS.bam \
O=alignment.CSF.bam \
ADD_MATE_CIGAR=true ASSUME_SORTED=true

# e- fix read group information
java -jar picard.jar AddOrReplaceReadGroups \
I=alignment.CSF.bam \
O=alignment.CSFR.bam \
RGID=id \
RGLB=library \
RGPL=illumina \
RGPU=lane \
RGSM=sample

# f- locates and tags duplicate reads in bam file
java -jar picard.jar MarkDuplicates \
I=alignment.CSFR.bam \
O=alignment.CSFRD.bam \
M=alignment.CSFRD.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
ASSUME_SORTED=true

# g- index bam file
samtools index alignment.CSFRD.bam
```


#### 4) Variant calling using [GATK](https://gatk.broadinstitute.org/hc/en-us)

```bash

# a- generate gVCFs
gatk HaplotypeCaller \
-R Cobs2.1.reference.fa \
-I alignment.CSFRD.bam \
-O alignment.CSFRD.g.vcf.gz \
-ERC GVCF -ploidy 2

# b- combine gVCFs
gatk CombineGVCFs \
-R Cobs2.1.reference.fa \
-V input.list \
-O all.g.vcf.gz

# c- genotype the combined gVCF
gatk GenotypeGVCFs \
-R Cobs2.1.reference.fa \
-V all.g.vcf.gz \
-O al.vcf.gz
```

#### 5) Apply hard filters using [GATK](https://gatk.broadinstitute.org/hc/en-us)

```bash
# a- it is better to visualize the variants attributes and decide for thresholds.
bcftools query -f '%INFO/DP\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\n' all.vcf.gz > info_density.tsv
# visualize each attribute in R (i.e., make density plots)

# b- get SNPs only
gatk SelectVariants -V all.vcf.gz --select-type-to-include SNP -O all_snps.vcf.gz

# c- apply hard filters:
gatk VariantFiltration \
-V all_snps.vcf.gz  \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -5.0" --filter-name "MQRankSum-5" \
-filter "ReadPosRankSum < -5.0" --filter-name "ReadPosRankSum-5" \
-O all_snps_filtered.vcf.gz

# d- extract variants that have passed the applied filters
gatk SelectVariants -V all_snps_filtered.vcf.gz -select 'vc.isNotFiltered()' -O all_snpsPASS.vcf.gz
```
#### 6) More filtering using [VCFtools](https://vcftools.github.io/index.html)

Similar to above, generate and visualize annotations associated with each variant to decide for cutoffs. Use this [script](./scripts/get.QCmetrics.from.vcf.sh) to generate a set of metrics that you can inspect using [R](https://www.r-project.org/). Now use VCFtools to apply filters.

```bash
vcftools \
--gzvcf all_snpsPASS.vcf.gz \
--recode --recode-INFO-all \
--max-missing X \
--minDP Y \
--maxDP YY \
--maf Z \
--min-alleles 2 \
--max-alleles 2 \
--out all_snpsPASS_filtered
```

### II- Population structure analyses

#### 1) Convert vcf to plink format using [VCFtools](https://vcftools.github.io/index.html) and [PLINK](https://www.cog-genomics.org/plink/)

```bash
# generate chromosomes map:
grep -v "^#" all_snpsPASS_filtered.vcf|cut -f1 | uniq| sort -V|awk '{print $0"\t"$0}' > chrom.map

# vcf >> plink
vcftools \
--vcf all_snpsPASS_filtered.vcf \
--plink \
--chrom-map chrom.map\
--out all_snpsPASS

# generate bed file
plink \
--file all_snpsPASS \
--make-bed \
--aec \
--out all_snpsPASS
```
#### 2) Linkage Disequilibrium (LD) pruning using [PLINK](https://www.cog-genomics.org/plink/)

```bash
# first generate a list of position to keep/remove
plink \
--bfile all_snpsPASS \
--aec \
--indep-pairwise 1 kb 1 0.2 \
--out all_snpsPASS.prune
# outputs two files: all_snpsPASS.prune.in (to keep) and all_snpsPASS.prune.out (to remove)

# then
plink \
--bfile all_snpsPASS \
--exclude all_snpsPASS.prune.out
--aec \
--out all_snpsPASS.LDpruned \
--make-bed
```


#### 3) PCA using [PLINK](https://www.cog-genomics.org/plink/) and admixture [ADMIXTURE](http://dalexander.github.io/admixture/)

```bash
# PCA
plink \
--bfile all_snpsPASS.LDpruned \
--aec \
--out all_snpsPASS.LDpruned \
--pca 4 # to limit the analysis to the first 4 eigenvectors

# ADMIXTURE
$ for K in `seq 2 4`; do admixture --cv all_snpsPASS.LDpruned.bed $K | tee log${K}.out ; done
# to inspect CV error values: grep -h CV log*.out
```
### III- Demographic population history analysis
#### 1) Phasing using [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

```bash
# extract positions with no missing data
vcftools --gzvcf all_snpsPASS_filtered.vcf.gz \
--max-missing 1 \
--recode --recode-INFO-all \
--out all_snpsPASS_nomissing.vcf

# compress and index
bgzip all_snpsPASS_nomissing.vcf

tabix -p vcf all_snpsPASS_nomissing.vcf.gz

# get call in each scaffold individually
for i in {1..30}; do bcftools view -r scaffold_$i all_snpsPASS_nomissing.vcf.gz -Oz -o scaffold_$i.vcf.gz; done

# perform phasing (3 steps)
for i in {1..30}; do shapeit -check --input-vcf scaffold_$i.vcf.gz --output-log scaffold_$i.alignments;done

for i in {1..30}; do shapeit -phase -V scaffold_$i.vcf.gz -O scaffold_$i.phased.haps.gz scaffold_$i.phased.samples --output-log scaffold_$i.main -T 10;done

for i in {1..30}; do shapeit -convert --input-haps scaffold_$i.phased.haps.gz scaffold_$i.phased.samples --output-vcf scaffold_$i.Onlyphased.vcf --output-log scaffold_$i.convert -T 10;done
```

#### 2) Demographic population history analysis [MSMC2](https://github.com/stschiff/msmc2)
To run MSMC2, check ([Schiffels & Wang, 2020](https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_7)) for a very nice and detailed description of how to estimate effective population size (_N_<sub>e</sub>) and the relative cross-coalescence rate (rCCR). The [GitHub](https://github.com/stschiff/msmc2) repository has also some nice resources.

### IV- Population genomic analyses

#### 1) Computing nucleotide diversity (&#960;) and Tajima's _D_ from pool-seq data using [PoPoolation](https://sourceforge.net/projects/popoolation/)



```bash
# generate mpileup files
samtools mpileup pool.bam > pool.mpileup

# exclude InDels
perl popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input pool.mpileup --output pool.InDels.gtf

perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl  --input pool.mpileup --output pool.indelsfree.pileup --gtf pool.InDels.gtf

# calculate nucleotide diversity
perl popoolation_1.2.2/Variance-sliding.pl \
--input pool.indelsfree.pileup \
--output pool.D \
--snp-output pool.D \
--measure pi \
--window-size 100000 \
--step-size 10000 \
--min-qual 20 \
--pool-size 16 \ # 30 for itabuna
--min-covered-fraction 0 \
--min-count 2 \
--min-coverage 4 \
--max-coverage 79 \ # 75 for itabuna
--fastq-type sanger

# calculate Tajima's D
perl popoolation_1.2.2/Variance-sliding.pl \
--input pool.indelsfree.pileup \
--output pool.D \
--snp-output pool.D \
--measure D \ #
--window-size 100000 \
--step-size 10000 \
--min-qual 20 \
--pool-size 16 \ # 30 for itabuna
--min-covered-fraction 0 \
--min-count 2 \
--min-coverage 4 \
--max-coverage 79 \ # 75 for itabuna
--fastq-type sanger
```

#### 2) Computing genetic differentiation (_F_<sub>ST</sub>) from pool-seq data by [PoPoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/):

```bash
# construct the mpileup:
samtools mpileup -B pool_tenerife.bam pool_itabuna.bam > tenerife_itabuna.mpileup

# generate the sync file:
java -ea -Xmx12g -jar popoolation2_1201/mpileup2sync.jar \
--input tenerife_itabuna.mpileup \
--output tenerife_itabuna.sync \
--fastq-type sanger \
--min-qual 20 --threads 20

# calculate Fst
perl popoolation2_1201/fst-sliding.pl \
--input tenerife_itabuna.sync \
--output tenerife_itabuna_100kb10kb.fst \
--min-count 4 \
--min-coverage 20 \
--max-coverage 79,75 \
--min-covered-fraction 0 \
--window-size 100000 \
--step-size 10000 \
--pool-size 16:30 \
--suppress-noninformative
```

### V- Repeat quantification and TE insertions identification

#### 1) Running [dnaPipeTE](https://github.com/clemgoub/dnaPipeTE)

```bash
python3 ./dnaPipeTE.py \
-input tenerife.R1.fastq.gz \
-output /tenerife.0.1x.4it.RB.Dfam.lib \
-cpu 30 \
-genome_size 193051228 \
-genome_coverage 0.1 \
-sample_number 4 \
-RM_lib RB.Dfam.combined.lib

# for itabuna, same command using itabuna.R1.fastq.gz instead
```
* **_processing the reads_per_component_and_annotation_**
Estimate the proportion of the genome that each repeat makes

```bash

# get entries of classified repeats
awk 'FS=" " {if ($7>=0) print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' reads_per_component_and_annotation |tr "/" "\t"|awk '{if (!$8) {print  $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$6"\t"$7"\t"$8} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 }}' > tenerife.0.1x.4it.RB.Dfam.classified.txt

awk 'FS=" " {if ($7>=0) print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' reads_per_component_and_annotation |tr "/" "\t"|awk '{if (!$8) {print  $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$6"\t"$7"\t"$8} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 }}' > itabuna.0.1x.4it.RB.Dfam.classified.txt


# get entries of unclassified repeats
awk 'FS=" " {if (!$7) print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' reads_per_component_and_annotation |tr "/" "\t"|cut -f1,2,3 > tenerife.0.1x.4it.RB.Dfam.unclassified.txt

awk 'FS=" " {if (!$7) print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' reads_per_component_and_annotation |tr "/" "\t"|cut -f1,2,3 > itabuna.0.1x.4it.RB.Dfam.unclassified.txt

```
Then go to R:
```R
#get the no of reads from the blast sample in the Annotation folder
Tenerife <- 19173155
Itabuna <- 19173092

#load classified data
annoTenerife <- read.table("tenerife.0.1x.4it.RB.Dfam.classified.txt", fill = T)
colnames(annoTenerife) <- c("contig","no_reads","aligned_bases","contig_length","ann", "Order", "Family","RM_hit_coverage")

annoItabuna <- read.table("itabuna.0.1x.4it.RB.Dfam.classified.txt", fill = T)
colnames(annoItabuna) <- c("contig","no_reads","aligned_bases","contig_length","ann", "Order", "Family","RM_hit_coverage")

#load unclassified data
unclasTenerife <- read.table("tenerife.0.1x.4it.RB.Dfam.unclassified.txt", fill = T)
colnames(unclasTenerife) <- c("contig","no_reads","aligned_bases")

unclasItabuna <- read.table("itabuna.0.1x.4it.RB.Dfam.unclassified.txt", fill = T)
colnames(unclasItabuna) <- c("contig","no_reads","aligned_bases")


#Calculate proportions for classified
annoTenerife$proportions <- annoTenerife$aligned_bases/Tenerife
annoItabuna$proportions <- annoItabuna$aligned_bases/Itabuna

#Calculate proportions for unclassified
unclasTenerife$proportions <- unclasTenerife$aligned_bases/Tenerife
unclasItabuna$proportions <- unclasItabuna$aligned_bases/Itabuna

#only most abundant repeat families:
df <- data.frame(Family=c("Gypsy","Simple_repeat","Maverick", "R1","R2","Helitron", "Pao", "Low_complexity",
                          "TcMar-Tc","Copia", "CMC-Chapaev-3", "hAT-Blackjack","RTE-X", "Unclassified"),
                 Tenerife=c(sum(annoTenerife$proportions[annoTenerife$Family=="Gypsy"])*100, sum(annoTenerife$proportions[annoTenerife$Family=="Simple_repeat"])*100,
                            sum(annoTenerife$proportions[annoTenerife$Family=="Maverick"])*100, sum(annoTenerife$proportions[annoTenerife$Family=="R1"])*100,sum(annoTenerife$proportions[annoTenerife$Family=="R2"])*100,
                            sum(annoTenerife$proportions[annoTenerife$Family=="Helitron"])*100, sum(annoTenerife$proportions[annoTenerife$Family=="Pao"])*100,
                            sum(annoTenerife$proportions[annoTenerife$Family=="Low_complexity"])*100, sum(annoTenerife$proportions[annoTenerife$Family=="TcMar-Tc1" | annoTenerife$Family=="TcMar-Tc4"])*100,
                            sum(annoTenerife$proportions[annoTenerife$Family=="Copia"])*100, sum(annoTenerife$proportions[annoTenerife$Family=="CMC-Chapaev-3"])*100,
                            sum(annoTenerife$proportions[annoTenerife$Family=="hAT-Blackjack"])*100, sum(annoTenerife$proportions[annoTenerife$Family=="RTE-X"])*100,
                            sum(unclasTenerife$proportions)*100),
                 Itabuna=c(sum(annoItabuna$proportions[annoItabuna$Family=="Gypsy"])*100, sum(annoItabuna$proportions[annoItabuna$Family=="Simple_repeat"])*100,
                           sum(annoItabuna$proportions[annoItabuna$Family=="Maverick"])*100, sum(annoItabuna$proportions[annoItabuna$Family=="R1"])*100,sum(annoItabuna$proportions[annoItabuna$Family=="R2"])*100,
                           sum(annoItabuna$proportions[annoItabuna$Family=="Helitron"])*100, sum(annoItabuna$proportions[annoItabuna$Family=="Pao"])*100,
                           sum(annoItabuna$proportions[annoItabuna$Family=="Low_complexity"])*100, sum(annoItabuna$proportions[annoItabuna$Family=="TcMar-Tc1" | annoItabuna$Family=="TcMar-Tc4"])*100,
                           sum(annoItabuna$proportions[annoItabuna$Family=="Copia"])*100, sum(annoItabuna$proportions[annoItabuna$Family=="CMC-Chapaev-3"])*100,
                           sum(annoItabuna$proportions[annoItabuna$Family=="hAT-Blackjack"])*100, sum(annoItabuna$proportions[annoItabuna$Family=="RTE-X"])*100,
                           sum(unclasItabuna$proportions)*100)



df[, 2:3] <- round(df[, 2:3], digits = 2)
df
```




#### 2) Running [PoPoolatioTE2](https://sourceforge.net/projects/popoolation-te2/)
* **_generate a TE-merged genome and generate TE hierarchy_**
```bash
# mask ref genome
perl RepeatMasker/RepeatMasker \
-gccalc \
-s \
-cutoff 200 \
-no_is \
-nolow \
-norna \
-gff \
-u \
-pa 30 \
-lib TE.lib \ # TE library
Cobs2.1.reference.fa
```

```bash
# merge TE lib and ref genome
cat Cobs2.1.reference.fa.masked TE.lib > Cobs2.1.temergedref.fa
```

```bash
# generate TE hierarchy

cat TE.lib|grep '^>'|perl -pe 's/>//' |cut -f1 -d" "|tr "#" "\t"|tr "/" "\t" > te-hierarchy.txt
awk '{print $1"\t"$3"\t" $2}' te-hierarchy.txt > tmp.te && mv tmp.te te-hierarchy.txt
```

* **_Map reads to the TE-combined-reference_**

```bash
# Index the TE-masked ref genome
bwa index Cobs2.1.temergedref.fa

# map reads with bwa bwasw in batch
bwa bwasw -t 20 Cobs2.1.temergedref.fa pool.R1.fastq.gz >pool.R1.sam
```


* **_Restore paired-end information with PoPoolationTE2 se2pe_**

```bash

java -jar popoolationte2/popte2.jar se2pe \
--fastq1 pool_R1_paired.fastq.gz \
--fastq2 pool_R2_paired.fastq.gz \
--bam1 pool.R1.sam \
--bam2 pool.R2.sam \
--sort \
--output population1.sorted.bam
```

* **_Generate the ppileup file_**

```bash

java -jar popoolationte2/popte2.jar ppileup \
--bam tenerife.sort.bam \
--bam itabuna.sort.bam \
--hier te-hierarchy.txt \
--map-qual 15 \
--output tenerife.itabuna.ppileup.RM.gz

```

* **_Subsample the ppileup file_**

Subsampling an equal physical coverage is required to correct for insert size and read coverage differences between samples (see [Kofler et al. 2016](https://doi.org/10.1093/molbev/msw137)).

```bash

java -jar popoolationte2/popte2.jar subsamplePpileup \
--ppileup tenerife.itabuna.ppileup.RM.gz \
--target-coverage 20 \
--output ss.tenerife.itabuna.ppileup.RM.gz

```

* **_Run PoPoolationTE2_**

```bash

# first
java -jar popoolationte2/popte2.jar identifySignatures \
--ppileup ss.tenerife.itabuna.ppileup.RM.gz \
--mode separate \
--output ss.tenerife.itabuna.ppileup.RM.signatures --min-count 2

# then
java -jar popoolationte2/popte2.jar frequency \
--ppileup ss.tenerife.itabuna.ppileup.RM.gz \
--signature ss.tenerife.itabuna.ppileup.RM.signatures \
--output ss.tenerife.itabuna.ppileup.RM.freqsig

# and
java -jar popoolationte2/popte2.jar filterSignatures \
--input ss.tenerife.itabuna.ppileup.RM.freqsig \
--output ss.tenerife.itabuna.ppileup.RM.filtered.freqsig \
--max-otherte-count 2 \
--max-structvar-count 2

# finally
java -jar popoolationte2/popte2.jar pairupSignatures \
--signature ss.tenerife.itabuna.ppileup.RM.filtered.freqsig \
--ref-genome Cobs.alpha.2.1.temergedref.fasta \
--hier te-hierarchy.txt \
--min-distance -200 \
--max-distance 300 \
--output ss.tenerife.itabuna.ppileup.RM.filtered.teinsertions
```
* **_Processing the ss.tenerife.itabuna.ppileup.RM.filtered.teinsertions_**

```bash
# add -/+150 bp to the TE location calculated by PoPoolatioTE2
awk '{print $1 "\t" $2 "\t" $3-150 "\t" $3+150 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' ss.tenerife.itabuna.ppileup.RM.filtered.teinsertions > start-end.teinsertions

# get TE insertions identified in the 30 largest scaffold
awk '{gsub(/scaffold/,""); print}' start-end.teinsertions |awk '{if ($2 <= 30) print $0}'|sort -V -k2,2 -k3,3|awk '{print $1 "\tscaffold" $2 "\t" $3 "\t" $4 "\t" $5"\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11}' > start-end.teinsertions.sorted

# exctract insertions found in each population
awk '{if ($1 == 1) print $0}' start-end.teinsertions.sorted > tenerife.tmp.bed
awk '{if ($1 == 2) print $0}' start-end.teinsertions.sorted > itabuna.tmp.bed

# add ids (population+number) to make manual curation slightly easier
awk -F'\t' 'BEGIN {OFS = FS} {id++}{print $2,$3,$4,"tenerife"id,$6,$7,$8,$9,$10}' tenerife.tmp.bed > tmp && mv tmp tenerife.tmp.bed
awk -F'\t' 'BEGIN {OFS = FS} {id++}{print $2,$3,$4,"itabuna"id,$6,$7,$8,$9,$10}' itabuna.tmp.bed > tmp && mv tmp itabuna.tmp.bed

#Get shared entries
bedtools intersect -a tenerife.tmp.bed -b itabuna.tmp.bed -names  itabuna -sorted -wao -f 0.5 -r |awk '{if ($19>=150 && $6==$15) print $0}' > te.insertions.shared.between.tenerife.itabuna.txt

#Get unique entries
bedtools intersect -a tenerife.tmp.bed -b itabuna.tmp.bed -names  itabuna -sorted -wao -f 0.5 -r |awk '{if ($19==0 || $6!=$15) print $0}' > te.insertions.unique.tenerife.txt
bedtools intersect -a itabuna.tmp.bed -b tenerife.tmp.bed -names  tenerife -sorted -wao -f 0.5 -r |awk '{if ($19==0 || $6!=$15) print $0}' > te.insertions.unique.itabuna.txt

# manually check for double entries to get the final set of TEs:
```

