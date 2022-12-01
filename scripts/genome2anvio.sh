#!/bin/bash
#
#SBATCH --job-name=genomes2anvio
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --mail-user=eduard.fadeev@univie.ac.at
#SBATCH --output=/scratch/oceanography/efadeev/MEVIMBA/Alteromonas_phylogeny/Log/%x-%j.out
#SBATCH --error=/scratch/oceanography/efadeev/MEVIMBA/Alteromonas_phylogeny/Log/%x-%j.err

FILES=(./Genomes/*.gb)    
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
name="$(basename -- $FILE| sed 's/\..*//g')"

#parse genebank file
anvi-script-process-genbank --input-genbank ${FILE} -O ./anvio/${name}

#correct headers of contigs
#anvi-script-reformat-fasta ./anvio/${name}-contigs.fa -o ./anvio/${name}-contigs-fixed.fa -l 0 --simplify-names --seq-type NT

#generate anvio database
anvi-gen-contigs-database -f ./anvio/${name}-contigs.fa -n "${name}" --output-db-path ./anvio/${name}.db --external-gene-calls ./anvio/${name}-external-gene-calls.txt --ignore-internal-stop-codons

#run HMM annotation of the genes for clustering
anvi-run-hmms -c ./anvio/${name}.db --num-threads 1
