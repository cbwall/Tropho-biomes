#!/bin/bash
#PBS -q hotel
#PBS -N cutadapt
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:30:00
#PBS -o /projects/ps-shurinlab/users/cbwall/Yos_water_16S/output/outs/cutadapt_out
#PBS -e /projects/ps-shurinlab/users/cbwall/Yos_water_16S/output/errs/cutadapt_err
#PBS -V
#PBS -M cbwall@ucsd.edu
#PBS -m ae
#PBS -A shurin-group


# activate conda and cutadapt, if there are problems the culprit is probably one of these two lines
source ~/.bashrc
conda activate cutadaptenv


# change directory to easily access raw fastq sequence files
cd /projects/ps-shurinlab/Yos_water_DNA/16S_Yos_data

#remove any index files
# rm index*


# unzip sequence files, if necessary
gunzip *.fastq.gz


# create output file in personal folder
# mkdir /projects/ps-shurinlab/users/cbwall/Yos_water_16S/data/trimmed
# mkdir /projects/ps-shurinlab/users/cbwall/Yos_water_16S/data/filtered


# create text file with list of all sample names
# assumes samples are in folder data with month.year after, and have format of...
# sampleid_S###_L001_R1_001.fastq
ls /projects/ps-shurinlab/Yos_water_DNA/16S_Yos_data/ | sed 's/_L[[:digit:]]\+_.*//g' > /projects/ps-shurinlab/users/cbwall/Yos_water_16S/data/samples.txt


# --- back to the code...
# code written specifically for samples that only have the 16S adapter sequence on the 5' end 
# (this should be standard for our library runs)

### re: -g and -G
### this is 'regular 5' adapter' that is not anchored
### 5’ adapter is a piece of DNA ligated to the 5’ end of the DNA fragment of interest.
### For this type of adapter to be found, the adapter sequence needs to either appear in full
### somewhere within the read (internal match) or at the start (5’ end) of it
### In all cases, the adapter itself and the sequence preceding it is removed.

# -g ADAPTER for forward primer/adapter
# -G adapter for reverse primer/adapter
# -o for output of forward read
# -p for output of reverse read
# -m for minimum and -M for maximum read length

# works for samples with the file names in format: SAMPLEID_SAMPLENUMBER_L001_R(1/2)_001.fastq, 
# may need to adjust

# also, cutadapt output stats (“cutadapt_primer_trimming_stats.txt”), assess how things went

#### --- back to the code...


# change to directory where 'data' and 'output' folder is for the loop
cd /projects/ps-shurinlab/users/cbwall/Yos_water_16S/


#run the loop


for SAMPLEID in $(cat data/samples.txt);
do
        echo "On sample: $SAMPLEID"
	cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT \
	-o data/trimmed/${SAMPLEID}_R1_trimmed.fastq \
	-p data/trimmed/${SAMPLEID}_R2_trimmed.fastq \
	-m 1 \
	/projects/ps-shurinlab/Yos_water_DNA/16S_Yos_data/${SAMPLEID}_L001_R1_001.fastq \
	/projects/ps-shurinlab/Yos_water_DNA/16S_Yos_data/${SAMPLEID}_L001_R2_001.fastq \
		>> output/cutadapt_primer_trimming_stats16S.txt 2>&1
done



# re-zip sequence files
gzip /projects/ps-shurinlab/Yos_water_DNA/16S_Yos_data/*.fastq
gzip /projects/ps-shurinlab/users/cbwall/Yos_water_16S/data/trimmed/*.fastq

# final output is a set of trimmed fastq files in a "trimmed" folder
