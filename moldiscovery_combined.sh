#! /bin/bash
strain=$1
db='prism'
if [[ $db = 'npatlas' ]]
then
db_pre="--db-path ../soil_metagenomics/npdtools/$db --preprocessed_targets ../soil_metagenomics/data/QTOF/moldiscovery_${db}/db_preproc/library.info.prob.mz.target.bin \
--preprocessed_decoys ../soil_metagenomics/data/QTOF/moldiscovery_${db}/db_preproc/library.info.prob.mz.decoy.bin --preprocess"
else
db_pre="--db-path data/QTOF/databases/$strain"
fi
input_data="data/QTOF/scans/*$strain*_MS2.mzXML"
program='../soil_metagenomics/npdtools2/bin/moldiscovery'
general_options="--mode HH --fdr --ppm --fdr-limit 1"
$program.py \
--debug $general_options $input_data -o data/QTOF/db_screens/$db/$strain/ --threads 1 $db_pre
##sbatch --cpus-per-task=16 --ntasks=10 --partition=multinode --mem=200g src/npdtools.sh
##--threads $SLURM_CPUS_PER_TASK --preprocess