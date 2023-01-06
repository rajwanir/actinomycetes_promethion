library(tidyverse)

barcode_key = read_csv("data/flongleQC/plate2_flongle_06092021/barcode_key.csv")


#trial1 ####
canu_cmds = barcode_key %>% mutate(
  cmd = sprintf("canu -p %s -d data/assemblies/canu/%s -genomeSize=10m -nanopore data/promethion/plate2/reads/fastq_pass/barcode%02d/*.fastq.gz useGrid=1 corMhapSensitivity=high corMinCoverage=0 gridOptions=\"--time=6:00:00\"",
                strain,strain,barcode)
) %>% select(cmd)


#trial2######

canu_trial2 = read_csv("src/swarms/canu_folders_trial2.csv",col_names = "canu_folder") %>% 
  mutate(strain = str_extract(canu_folder,pattern = "\\w+-\\w+"))

barcode_key = barcode_key %>% filter(strain %in% canu_trial2$strain)

#changing the estimated genome size (10-->8), include fail reads, wall time-->36, Stoponlowoverage<-1
canu_cmds = barcode_key %>% mutate(
  cmd = sprintf("canu -p %s -d data/assemblies/canu/%s -genomeSize=8m -nanopore data/promethion/plate2/reads/fastq_*/barcode%02d/*.fastq.gz useGrid=1 corMhapSensitivity=high corMinCoverage=0 minInputCoverage=1 stopOnLowCoverage=1 gridOptions=\"--time=36:00:00\"",
                strain,strain,barcode)
) %>% select(cmd)


write_tsv(canu_cmds,file = "src/swarms/canu_trial2.swarm",quote_escape = "none",col_names = F)

#swarm --gb-per-process 8 --module canu/2.1 --logdir data/assemblies/canu/swarm_logs --file src/swarms/canu_trial2.swarm



barcode_key %>% mutate(
  cmd = sprintf("bash src/racon.sh %s %s",
                sprintf("barcode%02d",barcode),
                strain)
) %>% select(cmd)



### assm polish

canu_assm_success = str_extract(Sys.glob("data/assemblies/canu/*/*.contigs.fasta"),pattern = "\\w+-\\w+")

barcode_key %>% filter(strain %in% canu_assm_success) %>% 
  mutate(
  cmd = sprintf("bash src/racon.sh %s %s;bash src/medaka.sh %s %s;bash src/homopolish.sh %s",
                sprintf("barcode%02d",barcode),strain,
                sprintf("barcode%02d",barcode),strain,
                strain)
) %>% select(cmd)

assmpolish_cmds = barcode_key %>% filter(strain %in% canu_assm_success) %>% 
  mutate(
    cmd = sprintf("bash src/racon.sh %s %s;bash src/medaka.sh %s %s;bash src/homopolish.sh %s",
                  sprintf("barcode%02d",barcode),strain,
                  sprintf("barcode%02d",barcode),strain,
                  strain)
  ) %>% select(cmd)


write_tsv(assmpolish_cmds,file = "src/swarms/assmpolish.swarm",quote_escape = "none",col_names = F)
#swarm --threads-per-process 8 --gb-per-process 32 --logdir data/assemblies/medaka/swarm_logs --file src/swarms/assmpolish.swarm --gres=lscratch:10

homopolish_cmds = barcode_key %>% filter(strain %in% canu_assm_success) %>% 
  distinct(strain) %>% 
  mutate(
    cmd = sprintf("bash src/homopolish.sh %s",
                  strain)
  ) %>% select(cmd)

write_tsv(homopolish_cmds,file = "src/swarms/homopolish.swarm",quote_escape = "none",col_names = F)

#swarm --threads-per-process 8 --gb-per-process 16 --logdir data/assemblies/homopolish/swarm_logs --file src/swarms/homopolish.swarm --gres=lscratch:10


## antismash ####

antismash_cmds = barcode_key %>% filter(strain %in% canu_assm_success) %>% 
  distinct(strain) %>% 
  mutate(
    cmd = sprintf("bash src/antismash.sh %s",
                  strain)
  ) %>% select(cmd)

write_tsv(antismash_cmds,file = "src/swarms/antismash.swarm",quote_escape = "none",col_names = F)

