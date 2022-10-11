####################################
#Preparations
####################################

#connect to the SCC with port forwarding
ssh -L 5678:localhost:5678 slurm

#activate environment
conda activate /apps/anvio/7

#load modules 
module load ncbiblast/2.2.26
module load ncbiblastplus/2.11.0

#define working directory
WORK_DIR="/scratch/oceanography/efadeev/MEVIMBA/Alteromonas_phylogeny"

cd $WORK_DIR

mkdir Log
mkdir anvio

####################################
#import genomes to anvio and annotate them
####################################
#generate anvio database from genome and annotate it using anvio
sbatch --array=0-5 ../source/genome2anvio.sh

#move interproscan parsing script into the home directory 
cp /scratch/oceanography/efadeev/anvio-resources/InterProScanParser/iprs2anvio.sh ./anvio/

#run functional annotation of the genomes
sbatch --array=0-5 ../source/genome_annotation.sh

#remove the parsing script and temportary files
rm ./anvio/iprs2anvio.sh ./anvio/*functions-*


#Export genes for proteome analysis
sbatch --array=0-5 ../source/export4R.sh