#!/bin/bash -l

#$ -l h_rt=70:0:0

#$ -l mem=20G

#$ -M ucfaeef@ucl.ac.uk -m abe

#$ -N SVM_20rep_Lfla

#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -e /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -o /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

# Run your R program

module unload compilers
module unload mpi
module load r/recommended

Rscript leave-one-species-out-svm.R Liostenogaster_flavolineata SVM_all_species Polistes_canadensis Polistes_dominula Ceratina_australensis Ceratina_calcarata Megalopta_genalis

