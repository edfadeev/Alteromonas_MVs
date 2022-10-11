#!/bin/bash
#
#SBATCH --job-name=functional_annotation
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=/scratch/oceanography/efadeev/MEVIMBA/Alteromonas_phylogeny/Log/%x-%j.out
#SBATCH --error=/scratch/oceanography/efadeev/MEVIMBA/Alteromonas_phylogeny/Log/%x-%j.err

FILES=(./anvio/*.db)    
FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}
name="$(basename -- ${FILE}| sed 's/\..*//g')"

#add the original NCBI PGAP annotation
anvi-import-functions -c ./anvio/${name}.db -i ./anvio/${name}-external-functions.txt

#COG annotation
anvi-run-ncbi-cogs -c ./anvio/${name}.db --num-threads 2 --cog-data-dir /scratch/oceanography/efadeev/anvio-resources/anvio-COG

#export AAs of all the gene calls
anvi-get-sequences-for-gene-calls -c ${FILE} --get-aa-sequences -o ./anvio/${name}-proteins.fasta 

#add annotations from interproscan
/scratch/oceanography/efadeev/anvio-resources/my_interproscan/interproscan-5.52-86.0/interproscan.sh -i ./anvio/${name}-proteins.fasta -f tsv --cpu 2 --goterms --pathways --application Pfam -o ./anvio/${name}-interpro-output.tsv -dp 

#parse output from Interproscan
cd ./anvio/

iprs2anvio.sh -i ${name}-interpro-output.tsv \
-o ${name}-interpro-output -g -p -r

cd ..

#import to anvio
anvi-import-functions -c ${FILE} -i ./anvio/${name}-interpro-output_iprs2anvio.tsv

#add annotations for KEGG modules
anvi-run-kegg-kofams -c ${FILE} \
--kegg-data-dir /scratch/oceanography/efadeev/anvio-resources/KOfam -T 2 \
--just-do-it

#produce tables with KO hits and modules
anvi-estimate-metabolism -c ${FILE} \
-O ./anvio/${name}-Kofam \
--kegg-data-dir /scratch/oceanography/efadeev/anvio-resources/KOfam \
--kegg-output-modes kofam_hits,modules
