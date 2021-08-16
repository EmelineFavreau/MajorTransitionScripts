#!/bin/bash -l

#$ -l h_rt=70:0:0

#$ -l mem=20G

#$ -M ucfaeef@ucl.ac.uk -m abe

#$ -N Ccal_bee_regression

#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -e /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -o /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm



module unload compilers
module unload mpi
module load r/recommended

Rscript leave-one-bee-out-svm.R Ceratina_calcarata Ccal_bee_regression Ceratina_australensis Megalopta_genalis
