#! /bin/bash
gbk_file=$1
file_path=$(dirname $gbk_file)
gff=$(basename --suffix '.gbk' $gbk_file)
source /gpfs/gsfs11/users/rajwanir2/conda/etc/profile.d/conda.sh
conda activate emboss
seqret -auto $gbk_file -feature yes -outseq $file_path/$gff.gff -osformat2 gff -sformat gb
sed '/\*/,/: / s/: /=/g' $file_path/$gff.gff | sed 's/note=\*//g' > $file_path/$gff.gff.1
sed -i '/##FASTA/,/##gff/d' $file_path/$gff.gff.1
##gpa_files=$(find data/ -wholename 'data/antismash/gpa/*/*.region*.gbk')
##for g in $gpa_files;do bash src/gb2gff.sh $g;done