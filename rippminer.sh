#! /bin/bash
strain=$1
ripp_class='lantibiotic'
#input_data=$(find data/QTOF/scans/ -wholename "*$strain*_MS2.mzXML")
input_data2="data/QTOF/scans/msdiscovery/*.mzXML"
program='../soil_metagenomics/npdtools/bin/metaminer' ##dereplicator+ or varquest or metaminer or moldiscovery
grep -A1 --no-group-separator $strain data/misc/ripps.fasta > temp/$strain.fasta
metaminer_options="--class $ripp_class --sequence temp/$strain.fasta --max-charge 3"
general_options="--mode HH --fdr"
$program.py \
--debug $general_options $input_data -o data/ripps/metaminer_confirmation/$ripp_class/$strain --threads $SLURM_CPUS_PER_TASK  $metaminer_options --reuse