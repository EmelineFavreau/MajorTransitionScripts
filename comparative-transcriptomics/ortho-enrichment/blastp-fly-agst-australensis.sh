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
#$ -wd /lustre/home/ucfaeef/projects/MajorTransitionScripts/comparative-transcriptomics/ortho-enrichment

# blasting proteins of Ceratina australensis against Drosophila proteins

module load blast+/2.2.30/intel-2015-update2


# Blast species against fly
blastp -query inputOnScratch/Ceratina_australensis-longest-isoforms.faa \
	-db inputOnScratch/Drosophila_melanogaster.faa \
	-out resultOnScratch/Ceratina_australensis \
	-max_target_seqs 1 \
	-outfmt 6 \
	-num_threads 8