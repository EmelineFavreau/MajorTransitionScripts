#!/bin/bash -l

# Request wallclock time 
#$ -l h_rt=5:0:0

# Request RAM.
#$ -l mem=10G

# Request TMPDIR space (default is 10 GB).
#$ -l tmpfs=10G

# Select the number of threads.
#$ -pe mpi 8

# Set working directory
#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/svm

# blasting 127 predictor genes (in Apis mellifera) against all nr

module load blast+/2.2.30/intel-2015-update2


# Blast
tblastn -query result/Apis_mellifera_127_predictor_genes.faa \
        -db nt \
        -out result/Apis_mellifera_127_predictor_genes \
        -max_target_seqs 1 \
        -outfmt 6 \
        -num_threads 8
