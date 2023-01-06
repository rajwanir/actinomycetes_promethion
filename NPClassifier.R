library(tidyverse)
library(parallel)

npatlas_matches = lapply(Sys.glob("data/QTOF/db_screens/npatlas/*/significant_unique_matches.tsv"),
       read_tsv) %>% bind_rows()

npatlas_matches = read_tsv("data/QTOF/db_screens/npatlas/all_sigmatches_uniq.tsv")

npatlas_matches = npatlas_matches %>% mutate(carbonyl = str_detect(SMILES, pattern = "C\\(=O\\)N"))

npclassifier = lapply(npatlas_matches$SMILES,function(s){
  #print(which(str_detect(npatlas_matches$SMILES,pattern = s)))
  gnps_list = tryCatch(jsonlite::fromJSON(sprintf("https://npclassifier.ucsd.edu/classify?smiles=%s",s)),
                       error=function(e) return(tibble(class_results = NA_character_,
                                                       superclass_results = NA_character_,
                                                       pathway_results = NA_character_,
                                                       isglycoside = NA_real_))) 
  if(!is_tibble(gnps_list)) {
    gnps_list = gnps_list %>% .[1:4] %>% lapply(.,'[',1) %>%  as_tibble() %>% .[1,]}
return(gnps_list)})

npclassifier=data.table::rbindlist(npclassifier)
npclassifier %>%  readr::type_convert()

write_tsv(npclassifier,file = "data/QTOF/db_screens/npatlas/npclassifier.tsv")