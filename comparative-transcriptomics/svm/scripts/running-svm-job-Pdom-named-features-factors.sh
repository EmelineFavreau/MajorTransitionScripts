#!/bin/bash -l

#$ -l h_rt=70:0:0

#$ -l mem=20G

#$ -M ucfaeef@ucl.ac.uk -m abe

#$ -N Pdom_classifier_all

#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -e /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -o /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

# Run your R program

module unload compilers
module unload mpi
module load r/recommended

Rscript leave-one-out-svm-factors.R Polistes_dominula Pdom_classifier_all Ceratina_australensis Ceratina_calcarata Megalopta_genalis Polistes_canadensis Liostenogaster_flavolineata

