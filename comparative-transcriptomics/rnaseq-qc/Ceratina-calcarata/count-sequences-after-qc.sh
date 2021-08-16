#!/bin/bash

# Calculate sequence number after Nextflow Trimming step

# Copyright Emeline Favreau, UCL

# files/path needed as input
thispath="/lustre/scratch/scratch/ucfaeef/comparative-transcriptomics/Ceratina-calcarata/result/star-all-samples/trim_galore/FastQC"

speciespath="/lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/rnaseq-qc/Ceratina-calcarata"

# after trimming
thisanalysis="qc-after-trim"

# dataset name
datasetname="ShellRehan2019"

# make a directory for these calculations
mkdir -p tmp/${datasetname}/${thisanalysis}/
mkdir -p result/${datasetname}/${thisanalysis}/

# copy the zip files
cp ${thispath}/*.zip tmp/${datasetname}/${thisanalysis}/.

# unzip files
cd tmp/${datasetname}/${thisanalysis}

ls > sample-file-names

for filename in `cat sample-file-names`; do
	unzip ${filename}
done

cd ${speciespath}

# obtain list of files e.g. SRR13281979_1_val_1
ls tmp/${datasetname}/${thisanalysis}/*_fastqc.zip \
	| cut -d "/" -f 4 \
	| cut -d "_" -f 1,2,3,4 \
	> tmp/${datasetname}/${thisanalysis}/file-list

# run a loop to obtain number of reads e.g. Total Sequences	15775938
for sample in `cat tmp/${datasetname}/${thisanalysis}/file-list`; do
	grep "Total Sequences" tmp/${datasetname}/${thisanalysis}/${sample}_fastqc/fastqc_data.txt >> tmp/${datasetname}/${thisanalysis}/number-reads-in-file
done

# copy information in a summary text file
paste tmp/${datasetname}/${thisanalysis}/file-list \
	tmp/${datasetname}/${thisanalysis}/number-reads-in-file \
	> tmp/${datasetname}/${thisanalysis}/number-reads-in-files


# obtain list of samples e.g SRR13281979
cat tmp/${datasetname}/${thisanalysis}/file-list \
	| cut -d "_" -f 1 \
	| uniq > tmp/${datasetname}/${thisanalysis}/samples-list

# run a loop to obtain number of reads

for sample in `cat tmp/${datasetname}/${thisanalysis}/samples-list`; do
	grep ${sample} tmp/${datasetname}/${thisanalysis}/number-reads-in-files | awk '{sum+=$4}END{print sum}' >> tmp/${datasetname}/${thisanalysis}/number-reads-in-sample

done

paste tmp/${datasetname}/${thisanalysis}/samples-list  tmp/${datasetname}/${thisanalysis}/number-reads-in-sample > result/${datasetname}/${thisanalysis}/number-reads-in-samples

#sum the read number for the whole experiment: 385141254
awk '{sum+=$2}END{print sum}' result/${datasetname}/${thisanalysis}/number-reads-in-samples