library(tidyverse)

biocat_results = lapply(Sys.glob("data/QTOF/db_screens/biocat/*/Results*.tsv"),function(x)
                        read_tsv(x,col_types = cols(.default = "c"))) 
names(biocat_results) = Sys.glob("data/QTOF/db_screens/biocat/*/Results*.tsv")
biocat_results = bind_rows(biocat_results,.id="biocat_path")

biocat_results = biocat_results %>% janitor::clean_names() %>% 
  mutate(relative_score = as.numeric(relative_score)) 

write_tsv(biocat_results,"data/QTOF/db_screens/biocat/allbiocat_results.tsv")


biocat_results = read_tsv("data/QTOF/db_screens/biocat/allbiocat_results.tsv")


biocat_results_filt = biocat_results %>% 
  filter(relative_score>0.8) %>% arrange(-relative_score) %>% 
  distinct(substance,biocat_path,.keep_all=T)


npatlas = read_tsv("../databases/NPatlas_01132020_inchikey.tsv")

biocat_results_filt = left_join(biocat_results,npatlas,
          by=c("substance"="npaid")) %>% 
  filter(relative_score>0.8) %>% arrange(-relative_score) %>% 
  distinct(substance,biocat_path,.keep_all=T)


select(biocat_results_filt,
       c(id=substance,smiles=nonsterosmiles)) %>% jsonlite::toJSON() %>% write(file="figures/rban_input_biocat_filtered.json")

svg(filename = "figures/biocat/hist_relscores.svg")
hist(biocat_results$relative_score,
     xlab = "Score", col="brown", main="Histogram of alignment scores (biocat) \nNRPS A-domain specificity to detected product ")
abline(v=0.8,col="grey",lwd=3, lty=2)
dev.off()