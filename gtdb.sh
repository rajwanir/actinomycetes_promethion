#! /bin/bash
ml R/4.0
R -e 'gtdb_tbl=dplyr::tibble(FASTA_file=Sys.glob("data/assemblies/homopolish/*/*_homopolished.fasta"),genome_ID=stringr::str_extract(pattern = "\\w+-\\w+",FASTA_file));readr::write_tsv(gtdb_tbl,"data/gtdb/batchfile.tsv",col_names=F)'
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate gtdbtk
gtdbtk classify_wf --batchfile data/gtdb/batchfile.tsv --out_dir data/gtdb/ --cpus $SLURM_CPUS_PER_TASK --scratch_dir data/gtdb/temp/