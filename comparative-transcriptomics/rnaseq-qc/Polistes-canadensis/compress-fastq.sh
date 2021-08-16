#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -t 1-20
#$ -tc 10

INPUT_FILE=$(sed -n "${SGE_TASK_ID}p" list_of_Patalano_files.txt)

# compress them all 
gzip inputOnScratch/Patalano/$INPUT_FILE
