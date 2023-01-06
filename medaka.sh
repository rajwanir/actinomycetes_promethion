#!/bin/bash
barcode=$1
strain=$2
racon_assm="data/assemblies/racon/$strain/$barcode.3.fasta"
reads="data/promethion/plate2/reads/fastq_*/$barcode/*.fastq.gz"
output_folder="data/assemblies/medaka/$strain"
module load medaka/1.2.0
medaka_consensus -d $racon_assm -o $output_folder -b 150 -m r941_prom_high_g4011 -t 8 -i $reads 
mv "data/assemblies/medaka/$strain/consensus.fasta" "data/assemblies/medaka/$strain/$strain.medaka.fasta"
rm "data/assemblies/medaka/$strain/calls_to_draft.bam";rm "data/assemblies/medaka/$strain/consensus_probs.hdf"