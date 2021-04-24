#!/bin/bash

################################ FASTQC ################################
# saves concatenated, trimmed files in a new directory
# run from analysis directory
# save large files that don't need to be backed up in **redacted**
# When running this file, use the command nohup rna-seq_test.sh >> rna-seq_test.out 2>&1 & tail -f rna-seq_test.out to keep a run log of the pipeline

################################ Setting up File Structure ################################

# Change the batch based on whatever the name is of the raw file directory
# directories removed prior to uploading to GitHub
set CURRENT_BATCH = "2019-05-08_order_5638"
set RAWS = "/.../$CURRENT_BATCH"
set WORKDIR = "/.../RNA-seq"
set TRIM = "$WORKDIR/trim/$CURRENT_BATCH"
set CAT = "$WORKDIR/tempcat/$CURRENT_BATCH"

# set up directories if needed
if [! -d "$WORKDIR/FASTQC_orig/$CURRENT_BATCH"] then
  mkdir $WORKDIR/FASTQC_orig
  mkdir $WORKDIR/FASTQC_orig/$CURRENT_BATCH
  echo "Creating directory $WORKDIR/FASTQC_orig/$CURRENT_BATCH"
fi

if (! -d "$WORKDIR/FASTQC_trim" ) then
  mkdir $WORKDIR/FASTQC_trim
  mkdir $WORKDIR/FASTQC_trim/$CURRENT_BATCH
  echo "Creating directory $WORKDIR/FASTQC_trim/$CURRENT_BATCH"
fi

if (! -d "$TRIM") then
  mkdir $WORKDIR/trim/
  mkdir $WORKDIR/trim/$CURRENT_BATCH
  echo "Creating directory $WORKDIR/trim/$CURRENT_BATCH"
fi

if (! -d "$CAT" ) then
  mkdir $WORKDIR/tempcat
  mkdir $WORKDIR/tempcat/$CURRENT_BATCH
  echo "Creating directory $WORKDIR/tempcat/$CURRENT_BATCH"
fi


################################ FastQC ################################

# Get a list of all of the samples in the raw file directory
cd $RAWS
ls *.gz > raw_sample_list.txt
mv raw_sample_list.txt $WORKDIR
cd $WORKDIR

# Remove the extra file name information to get a list of samples, then remove the duplicates
sed 's/\_.*//' raw_sample_list.txt > cut_raw_sample_list.txt
sort -u cut_raw_sample_list.txt > raw_sample_cut_dup.txt

# Initial FastQC on every individual fastq file
for sample in $(cat raw_sample_cut_dup.txt) 
do
    if (! -e /$WORKDIR/FASTQC_orig/$CURRENT_BATCH/{$sample}_R1_fastqc.html) then
    echo "initial fastqc for sample {$sample}_R1"
    fastqc $RAWS/{$sample}.fastq.gz --outdir $WORKDIR/FASTQC_orig/$CURRENT_BATCH
    else
    echo "an initial FASTQC report is already available for {$sample}_R1, moving to the next"
    fi
    if (! -e /$WORKDIR/FASTQC_orig/$CURRENT_BATCH/{$sample}_R2_fastqc.html) then
    echo "initial fastqc for {$sample}_R2"
    fastqc $RAWS/{$sample}.fastq.gz --outdir $WORKDIR/FASTQC_orig/$CURRENT_BATCH
    else
    echo "an initial FASTQC report is already available for {$sample}, moving to the next"
    fi
done

################################ cutadapt ################################ 

#Cutting out adapters, filtering based on quality
for sample in $(cat raw_sample_cut_dup.txt)
do
    echo "running cutadapt for sample {$sample}"
        if (! -e $TRIM/{$sample}_R1.trim.fastq.gz) then
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -e 0.1 -o $TRIM/{$sample}_R1.trim.fastq.gz -p $TRIM/{$sample}_R2.trim.fastq.gz $RAWS/{$sample}_R1.fastq.gz $RAWS/{$sample}_R2.fastq.gz
        else
        echo "Cutadapt has already been performed for the sample $sample, moving onto the next."
    fi
echo "cutadapt complete for sample {$sample}"
end

################################ FastQC_r2 ################################

for sample in $(cat raw_sample_cut_dup.txt)
    if (! -e /$WORKDIR/FASTQC_trim/$CURRENT_BATCH/{$sample}_R1.trim.fastqc.html) then
        echo "post-trimming fastqc for {$sample}_R1"
        fastqc $TRIM/{$sample}_R1.trim.fastq.gz --outdir $WORKDIR/FASTQC_trim/$CURRENT_BATCH
        else
        echo "a post-trimming FastQC report is already available for sample {$sample}_R1, moving to the next"
    fi
    if (! -e /$WORKDIR/FASTQC_trim/$CURRENT_BATCH/{$sample}_R2.trim.fastqc.html) then
        echo "post-trimming fastqc for {$sample}_R2"
        fastqc $TRIM/{$sample}_R2.trim.fastq.gz --outdir $WORKDIR/FASTQC_trim/$CURRENT_BATCH
        else
        echo "a post-trimming FastQC report is already available for sample {$sample}_R2, moving to the next"
    fi
end

################################ Concatenation ################################
# it appears this is unnecessary for aligning with STAR, which takes both reads as inputs  

# for sample (`cat raw_sample_cut_dup.txt`)
    # if (! -e $CAT/{$sample}_cat.fastq.gz) then
        # echo "concatenating the trimmed fastq files for sample {$sample}" 
        # cat $TRIM/{$sample}_R1.trim.fastq.gz $TRIM/{$sample}_R2.trim.fastq.gz > $CAT/{$sample}_cat.fastq.gz
        # else
        # echo "The trimmed fastq files have already been concatenated for the sample {$sample}, moving onto the next."
    # fi
# end


nohup fastqc *.gz >> fastqc_post_trim.out 2>&1 & tail -f fastqc_post_trim.out


################################ Alignment using STAR ################################
for i in $(cat raw_sample_cut_dup.txt)
do
STAR --runThreadN 20 \
--genomeDir /dmpi/analysis/SGregory/madi/index \
--readFilesIn /dmpi/analysis/SGregory/steve/unc_mouse_collaboration/RNA-seq/trim/"$i"_R1_trimmed.fastq.gz /dmpi/analysis/SGregory/steve/unc_mouse_collaboration/RNA-seq/trim/"$i"_R2_trimmed.fastq.gz \
--runMode alignReads --outFileNamePrefix /dmpi/analysis/SGregory/steve/unc_mouse_collaboration/RNA-seq/STAR_output/"$i"_ \
--genomeLoad LoadAndKeep --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat
done

STAR --genomeDir /dmpi/analysis/SGregory/madi/index --genomeLoad Remove

for i in $(cat raw_sample_cut_dup.txt)
do
featureCounts -T 20 -p -s 2 -a /dmpi/analysis/SGregory/madi/Mus_musculus.GRCm38.96.chr.gtf \
/.../RNA-seq/STAR_output/"$i"_Aligned.out.bam -o /.../RNA-seq/counts/"$i"_counts.txt

done 
