library(tidyverse)
library(data.table)

setDTthreads(threads = 0)

seqsummary = data.table::fread("data/promethion/plate2/reads/sequencing_summary_PAH04469_3621a629.txt")


estimated_coverage = seqsummary %>% group_by(barcode_arrangement) %>%
  summarize(
    total_bases = sum(sequence_length_template),
    total_bases_mb = total_bases /1000000,
    estimated_coverage = total_bases_mb /10
  )

write.table(estimated_coverage,
          "data/promethion/plate2/reads/estimated_coverage.tsv",
          sep = '\t',quote = F,row.names = F)


barcode_key = read_csv("data/flongleQC/plate2_flongle_06092021/barcode_key.csv") %>% mutate(barcode = sprintf("barcode%02d",barcode))

estimated_coverage = left_join(estimated_coverage,barcode_key,by = c("barcode_arrangement" = "barcode"))

assm = tibble(FASTA_file=Sys.glob("data/assemblies/canu/*/*.contigs.fasta"),genome_ID=stringr::str_extract(pattern = "\\w+-\\w+",FASTA_file))

estimated_coverage = left_join(estimated_coverage,assm,
          by=c("strain"="genome_ID"))


estimated_coverage %>% 
  filter(is.na(FASTA_file) & 
           !strain%in%c("GA5-005","GA9-001","GA5-002","GA3-004")) %>% 
  mutate(canu_folder = sprintf("data/assemblies/canu/%s/",strain)) %>% 
  select(canu_folder)