#! /bin/bash
ml OpenBabel
strain=$1
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate biocat
awk '{ print $13 " " $5}' data/QTOF/db_screens/npatlas/$strain/significant_unique_matches.tsv | tail -n+2 | obabel -s'[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]' -ismi -osmi | awk '{ print $2"\t"$1}' > temp/biocat/$strain.smi
biocat -file_smiles temp/biocat/$strain.smi -antismash data/antismash/$strain/*.json -out data/QTOF/db_screens/biocat/$strain -cpu 4 