#!/bin/bash -l

#$ -l h_rt=70:0:0

#$ -l mem=20G

#$ -M ucfaeef@ucl.ac.uk -m abe

#$ -N Pdom_all_regr_species_norm

#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -e /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

#$ -o /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm



module unload compilers
module unload mpi
module load r/recommended

Rscript R-scripts/leave-one-species-out-svm-species-normalised.R Polistes_dominula Pdom_all_regr_norm Polistes_canadensis  Liostenogaster_flavolineata Ceratina_australensis Ceratina_calcarata Megalopta_genalis

