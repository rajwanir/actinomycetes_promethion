#!/bin/bash
barcode=$1
strain=$2
ml minimap2
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate racon
reads="data/promethion/plate2/reads/fastq_*/$barcode/*.fastq.gz"
output_folder="data/assemblies/racon/$strain"
mkdir $output_folder
gunzip -c $reads > $output_folder/$strain.fastq
canu_assembly="data/assemblies/canu/$strain/$strain.contigs.fasta"
cp $canu_assembly $output_folder/${barcode}.0.fasta # copying as 0
###################################
#run three rounds of racon
for time in {1..3}
do
#align_reads
minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont -a $output_folder/${barcode}.$((time-1)).fasta $output_folder/$strain.fastq > $output_folder/alignment_${time}.sam
#run racon
racon -m 8 -x -6 --gap -8 --window-length 500 --threads $SLURM_CPUS_PER_TASK \
$output_folder/$strain.fastq $output_folder/alignment_${time}.sam ${output_folder}/${barcode}.$((time-1)).fasta > $output_folder/${barcode}.$((time)).fasta
done 
rm $output_folder/alignment_*.sam
rm $output_folder/$strain.fastq
conda deactivate