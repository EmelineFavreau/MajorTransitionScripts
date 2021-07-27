#!/bin/bash -l

#$ -l h_rt=70:0:0

#$ -l mem=20G

#$ -M ucfaeef@ucl.ac.uk -m abe

#$ -N SVM_20rep_Pdom_wasponly

#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -e /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -o /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

# Run svm with wasp only training data

module unload compilers
module unload mpi
module load r/recommended

Rscript leave-one-species-out-svm.R Polistes_dominula SVM_wasps_only Polistes_canadensis Liostenogaster_flavolineata

