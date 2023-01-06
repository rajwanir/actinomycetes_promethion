#! /bin/bash
strain=$1
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate antismash6
antismash "data/assemblies/homopolish/$strain/${strain}_homopolished.fasta" \
--cpus ${SLURM_CPUS_PER_TASK} \
--logfile "data/antismash/$strain/antismash.log" \
--output-dir "data/antismash/$strain" \
--cb-general \
--cb-subclusters \
--cb-knownclusters \
--cb-knownclusters \
--genefinding-tool prodigal-m \
--asf --pfam2go --smcog-trees \
--rre