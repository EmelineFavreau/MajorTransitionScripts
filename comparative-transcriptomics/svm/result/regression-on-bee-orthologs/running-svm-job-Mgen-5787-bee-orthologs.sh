#!/bin/bash -l

#$ -l h_rt=70:0:0

#$ -l mem=20G

#$ -M ucfaeef@ucl.ac.uk -m abe

#$ -N Mgen_bee_regression

#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -e /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -o /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm



module unload compilers
module unload mpi
module load r/recommended

Rscript leave-one-bee-out-svm.R Megalopta_genalis Mgen_bee_regression Ceratina_australensis Ceratina_calcarata
