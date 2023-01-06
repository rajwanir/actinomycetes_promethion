#!/bin/bash
folder=$1
set -e
module load guppy || exit 1
fastq_folder=$1
guppy_barcoder --input_path $folder/fastq --barcode_kits EXP-NBD196 \
--save_path $folder/demultiplexed \
--compress_fastq \
--records_per_fastq 0 \
--trim_barcodes \
--detect_mid_strand_barcodes \
-x cuda:all
##sbatch --partition=gpu --cpus-per-task=14 --mem=16g --gres=lscratch:200,gpu:v100x:2 src/guppy_demultiplex.sh data/flongleQC/plate2_flongle_06092021/