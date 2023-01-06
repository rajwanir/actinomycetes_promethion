#! /bin/bash
strain=$1
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate Prokka
prokka --cpus 8 \
"data/assemblies/homopolish/$strain/${strain}_homopolished.fasta" \
--outdir "data/prokka/$strain" 