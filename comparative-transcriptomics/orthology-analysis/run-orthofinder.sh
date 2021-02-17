#!/bin/bash -l

#$ -l h_rt=30:0:0

#$ -l mem=50G

#$ -N Orthofinder_3wasps_apis_droso

#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/orthology-analysis

#$ -e /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/orthology-analysis/tmp

module load python/2.7.9

module load boost/1_54_0/gnu-4.9.2

module load samtools/0.1.19

module load cufflinks/2.2.1

module add mcl/14-137

module add muscle/3.8.31

/lustre/home/ucfaeef/programs/OrthoFinder/orthofinder \
	-t 8 \
	-f /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/orthology-analysis/inputOnScratch/primary_transcripts \
	-A muscle \
	-S diamond \
	-M msa