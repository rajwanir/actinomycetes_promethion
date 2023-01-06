#! /bin/bash
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate checkm
checkm lineage_wf data/checkm/input/ data/checkm/ -t 32 -x .fasta