#! /bin/bash
ml OpenBabel
strain=$1
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate nerpa
printf "Name\tSMILES\n" > temp/nerpa/$strain.smi
awk '{ print $13 " " $5}' data/QTOF/db_screens/npatlas/$strain/significant_unique_matches.tsv | tail -n+2 | obabel -s'[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]' -ismi -osmi | awk '{ print $2"\t"$1}' >> temp/nerpa/$strain.smi
nerpa.py -a data/antismash/$strain --smiles-tsv temp/nerpa/$strain.smi --process-hybrids --output_dir data/QTOF/db_screens/nerpa/$strain --col-id Name --force-existing-outdir