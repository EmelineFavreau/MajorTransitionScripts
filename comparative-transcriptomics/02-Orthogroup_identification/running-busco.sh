#$ -l tmem=2G
#$ -l h_vmem=22G
#$ -l h_rt=02:00:00 
#$ -S /bin/bash
#$ -j y

source /share/apps/source_files/python/python-3.8.5.source

# busco assessement of C. australensis protein set that went in OrthoFinder
busco -m protein -i input/Ceratina_australensis-longest-isoforms.faa -o Caus_protein_test -l hymenoptera_odb10 -f  --config busco-config.ini
