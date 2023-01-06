#!/bin/bash
strain=$1
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate homopolish
python3 src/homopolish/homopolish.py polish -a data/assemblies/medaka/$strain/$strain.medaka.fasta -s src/homopolish/bacteria.msh -m R9.4.pkl -o data/assemblies/homopolish/$strain/ -t $SLURM_CPUS_PER_TASK
conda deactivate