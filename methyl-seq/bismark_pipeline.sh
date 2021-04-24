#!/bin/bash

################################ FASTQC ################################
# saves concatenated, trimmed files in a new directory
# run from analysis directory
# save large files that don't need to be backed up in /dmpi/analysis/
# When running this file, use the command nohup rna-seq_test.sh >> rna-seq_test.out 2>&1 & tail -f rna-seq_test.out to keep a run log of the pipeline

################################ Setting up File Structure ################################

# Change the batch based on whatever the name is of the raw file directory
# Directories removed prior to uploading to GitHub

CURRENT_BATCH="2019-10-24_order_5987"
RAWS=".../$CURRENT_BATCH"
WORKDIR="/.../sureselect"
TRIM="$WORKDIR/trim/$CURRENT_BATCH"
CAT="$WORKDIR/cat/$CURRENT_BATCH"
ALIGN="$WORKDIR/alignment/$CURRENT_BATCH"

# set up directories if needed
mkdir -p "$WORKDIR/FASTQC_orig/$CURRENT_BATCH"
mkdir -p "$WORKDIR/FASTQC_trim/$CURRENT_BATCH"
mkdir -p "$WORKDIR/trim/$CURRENT_BATCH"
mkdir -p "$WORKDIR/tempcat/$CURRENT_BATCH"
mkdir -p "$WORKDIR/alignment/$CURRENT_BATCH"

################################ FastQC ################################

# Get a list of all of the samples in the raw file directory
cd $RAWS
ls *.gz > raw_sample_list.txt
mv raw_sample_list.txt $WORKDIR
cd $WORKDIR


cd $RAWS
for f in *.gz; do
	fastqc $f -outdir $WORKDIR/FASTQC_orig/$CURRENT_BATCH;
done
cd $WORKDIR


################################ cutadapt ################################ 

# In the sureselect mouse guide the authors use the following parameters
# "Agilent SureSelectXT Methyl-Seq Applications with Low-Input DNA and Smaller Capture Libraries"
# It looks like their pipeline did not filter by quality at this step 

### NOT USING FLEXBAR, just for reference ####
flexbar -r IMR90_1_R1_PF.fastq -p # -r = read 1 input file, -p = read 2 input file
IMR90_1_R2_PF.fastq -t IMR90_1 -a # -a = adapters that need to be trimmed as a fasta file, -t = output file name
Truseq_adapters.fa -f i1.8 -n 10 -ao 6 # -f = fasta output, -n = number of threads to use, -ao = minimum overlap to remove adapter (default 3)
-m 21 -at 2 -ae RIGHT -u 2 -j -am 3 -ai # -m = minimum read length, -at = type of adapter removal, -u = maximum number of uncalled bases for each read, -j = length list, -am = alignment match score, -ai = adapter mismatch score
-3 -ag -20 # ag = alignment gap score 

### FILTERING READS ###
# FASTQC didn't show any adapter contamination

# Get the unique identifiers for the paired read of each sample 
grep -oP '^(?:[^_]*_){3}' raw_sample_list.txt > cut_raw_sample_list.txt
sort -u -o sample_list.txt cut_raw_sample_list.txt

# Run cutadapt on each sample for both paired ends
samples=`cat sample_list.txt`
for i in $samples
do
	echo "running cutadapt for"'$i'
	cutadapt  --pair-filter=any -m 21 -q 15,10 -o "$TRIM"'/'"$i"'R1.trim.fastq.gz' -p "$TRIM"'/'"$i"'R2.trim.fastq.gz' "$RAWS"'/'"$i"'R1_001.fastq.gz' "$RAWS"'/'"$i"'R2_001.fastq.gz'
done

################################ post-cutadapt FastQC ################################ 

# Final fastqc
cd $TRIM
for f in *.gz; do
	fastqc $f -outdir $WORKDIR/FASTQC_trim/$CURRENT_BATCH;
done
cd $WORKDIR

################################ Concatenating together the lanes ################################ 

# Make another list of samples to cat together the different lanes
cd $WORKDIR
grep -oP '^(?:[^_]*_){2}' raw_sample_list.txt > cat_sample_temp
sort -u -o cat_sample_list cat_sample_temp
rm -f cat_sample_temp

# Concatenate all of the lanes on the flowcell for each sample 

cd $WORKDIR
for i in `cat cat_sample_list`
do	
	echo 'concatenating for sample'"$i"
	cd $TRIM
	cat "$i"'L001_R1.trim.fastq.gz' "$i"'L001_R1.trim.fastq.gz' "$i"'L001_R1.trim.fastq.gz' "$i"'L001_R1.trim.fastq.gz' > "$i"'R1.cat.fastq.gz'
	cat "$i"'L001_R2.trim.fastq.gz' "$i"'L001_R2.trim.fastq.gz' "$i"'L001_R2.trim.fastq.gz' "$i"'L001_R2.trim.fastq.gz' > "$i"'R2.cat.fastq.gz'
	
done

mv *.cat* $CAT

################################ Alignment Using Bismarck and Tophat2 ################################
# https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf 


# Downloading the grcm38 fla file to create the index file
cd $ALIGN 
mkdir index
wget -P grcm38 https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
nohup bowtie2-build Mus_musculus.GRCm38.dna.primary_assembly.fa grcm38 >> bowtie2_build.out 2>&1 &
tail -f bowtie2_build.out
rm -f grcm38.fa

################################ Bismark Genome Preparation ################################

cd $ALIGN
mkdir -p methyl_genome_build
cd "ALIGN/methyl_genome_build"
nohup bismark_genome_preparation --bowtie2 --verbose "/dmpi/analysis/SGregory/steve/SureSelect/qc_check/alignment/2019-10-24_order_5987/index/grcm38/" >> bismark.build.out 2>&1 &
tail -f bismark.build.out

# THIS WILL NOT BE RUN, provided in the guidance by Agilent #
bismark -q -un --ambiguous -L 32 -X 500 /genomes/indices/Homo_sapiens/Bismark_bt2_hg19 -1 IMR90_1_R1.fastq -2 IMR90_1_R2.fastq
	## -q = input is fastq, 
	## -un = Write all reads that could not be aligned to the file _unmapped_reads.txt in the output directory
	## --ambiguous = Write all reads which produce more than one valid alignment with the same number of 
		# lowest mismatches or other reads that fail to align uniquely to _ambiguous_reads.txt
	## -L = Sets the length of the seed substrings to align during multiseed alignment. Smaller values 
		# make alignment slower but more senstive. Bowtie2 defaults to 20 (wonder why they went with bigger here)
	## -X = The maximum insert size for valid paired-end alignments. E.g. if -X 100 is specified and 
		# a paired-end alignment consists of two 20-bp alignments in the proper orientation with a 
		# 60-bp gap between them, that alignment is considered valid

nohup bismark -q -un --ambiguous --bowtie2 -L 32 -X 500 "/dmpi/analysis/SGregory/steve/SureSelect/qc_check/alignment/2019-10-24_order_5987/index/grcm38/" -1 "/dmpi/analysis/SGregory/steve/SureSelect/qc_check/cat/2019-10-24_order_5987/f9-h-13_S3_R1.cat.fastq.gz" -2 "/dmpi/analysis/SGregory/steve/SureSelect/qc_check/cat/2019-10-24_order_5987/f9-h-13_S3_R2.cat.fastq.gz" >> bismark_align.out 2>&1 &
tail -f bismark_align.out


nohup bismark_methylation_extractor -s --comprehensive g6-h-5_S2_R1.cat.fastq.gz_bismark_bt2_pe.sam >> methlation_extractor_g6-h-5_S2_.out 2>&1 &
tail -f methlation_extractor_g6-h-5_S2_.out

# Generating summary reports of the alignment
bismark2report # run in directory containing the results

# to run this
nohup cutadapt.sh >> cutadapt.sh.out 2>> cutadapt.sh.error.out & 
tail -f cutadapt.sh.out
