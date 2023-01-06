#! /bin/bash
ripp_class='lap'
input_data=$(find data/QTOF/scans/ -wholename "*_MS2.mzXML")
program='../soil_metagenomics/npdtools/bin/metaminer' ##dereplicator+ or varquest or metaminer or moldiscovery
sequence="data/misc/ripps_all.fasta" 
metaminer_options="--class $ripp_class --sequence $sequence --max-charge 3"
general_options="--mode HH"
$program.py \
$general_options $input_data -o data/ripps/metaminer_confirmation/ripps_allLAP_nodebug/ --threads $SLURM_CPUS_PER_TASK $metaminer_options --reuse