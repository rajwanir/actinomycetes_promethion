#!/bin/bash
folder=$1
set -e
module load guppy || exit 1
guppy_basecaller --input_path $folder/fast5 --flowcell FLO-FLG001 --kit SQK-LSK109 \
--save_path $folder/fastq \
--min_qscore 0 \
-x cuda:all \
--compress_fastq \
--records_per_fastq 0 \
--recursive
##example submission:: sbatch --partition=gpu --cpus-per-task=14 --mem=16g --gres=lscratch:200,gpu:v100x:2 src/guppy_basecall.sh data/flongleQC/plate2_flongle_06092021/