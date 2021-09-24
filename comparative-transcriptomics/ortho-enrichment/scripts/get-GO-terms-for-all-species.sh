#!/bin/bash -l

# Batch script to run a serial array job under SGE: get enriched GO terms

# Request 1h of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=1:00:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=10M

# Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=15G

# Set up the job array.
# wc -l species-list
#$ -t 1-6

# Set the working directory
#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/ortho-enrichment

# set the orthogroup at start of array
species=$(sed -n "${SGE_TASK_ID}p" species-list)

module unload compilers
module unload mpi
module load r/recommended

# example: Rscript makeGOTermsTable.R Ceratina_australensis R BP inputOnScratch/Ceratina_australensis_DESeq2_results.txt inputOnScratch/Ceratina_australensis_upR_DEGs.txt resultOnScratch/Ceratina_australensis_filtered


# loop through the datasets: R reproductives, NR non-reproductives
# loop through the Go Categories: 
# Biological Processes, Metabolic function, Cellular componenet
for dataset in R NR; do
	for Gocategory in BP MF CC; do
		Rscript makeGOTermsTable.R ${species} ${dataset} ${Gocategory} 
	done
done


# result should be a table 
# colnames: "GO.ID" "Term" "Annotated" "Significant" "Expected" "pvalue"      "species"  "dataSubset" "goCategory" 