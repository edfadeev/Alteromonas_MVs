#!/bin/bash
#
#SBATCH --job-name=annotations4R
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=/scratch/oceanography/efadeev/MEVIMBA/Alteromonas_phylogeny/Log/%x-%j.out
#SBATCH --error=/scratch/oceanography/efadeev/MEVIMBA/Alteromonas_phylogeny/Log/%x-%j.err

FILES=(./anvio/*.db)    
FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}
name="$(basename -- ${FILE}| sed 's/\..*//g')"

#export all gene calls
anvi-export-gene-calls -c ${FILE} -o ./anvio/${name}-genes.txt --gene-caller NCBI_PGAP

#export AAs of all the gene calls
anvi-get-sequences-for-gene-calls -c ${FILE} --get-aa-sequences -o ./anvio/${name}-proteins.fasta

#export functions of each gene
DATABASES=("COG20_PATHWAY" "COG20_CATEGORY" "COG20_FUNCTION" \
           "KOfam" "KEGG_Module" "KEGG_Class" \
           "Pfam" "GO" "InterPro" "MetaCyc" "NCBI_PGAP")

for db in ${DATABASES[@]}; do
anvi-export-functions -c ${FILE} --annotation-sources $db -o ./anvio/${name}-$db-functions-raw.txt

#remove spaces
sed -e 's/ /_/g' ./anvio/${name}-$db-functions-raw.txt > ./anvio/${name}-$db-functions-corrected.txt

#sort the table
sed -n '2,$p' ./anvio/${name}-$db-functions-corrected.txt | sort -n -k 1 > ./anvio/${name}-$db-functions.txt

rm ./anvio/${name}-$db-functions-raw.txt ./anvio/${name}-$db-functions-corrected.txt

done
