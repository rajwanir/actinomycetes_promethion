library(tidyverse)
library(parallel)

npatlas = read_tsv("../databases/NPatlas_01132020_inchikey.tsv")
sigmatches = lapply(Sys.glob("data/QTOF/db_screens/npatlas/*/significant_unique_matches.tsv"),read_tsv) %>% bind_rows()
sigmatches = left_join(sigmatches,npatlas,
                       by=c("Name"="npaid"))


structures = read_tsv("data/prism/prism_structures_npatlasann.tsv")


sigmatches = sigmatches %>% mutate(strain = str_extract(SpecFile,"G\\w+-\\d+"))

sigmatches$prism_similar_name = mclapply(1:nrow(sigmatches),function(idx){
  cmd = sprintf("obabel -ifs data/QTOF/databases/%s/smiles.fs -s\"%s\" -osmi -at1 -xt",
                sigmatches$strain[idx],
                sigmatches$nonsterosmiles[idx]
                )
  system(cmd,intern = T)
    },mc.cores = 60) %>% unlist()

sigmatches$prism_similar_smile = mclapply(1:nrow(sigmatches),function(idx){
  cmd = sprintf("obabel -ifs data/QTOF/databases/%s/smiles.fs -s\"%s\" -osmi -at1",
                sigmatches$strain[idx],
                sigmatches$nonsterosmiles[idx]
  )
  system(cmd,intern = T)
},mc.cores = 60) %>% unlist()

sigmatches$prism_similar_score = mclapply(1:nrow(sigmatches),function(idx){
  cmd = sprintf("obabel -:\"%s\" -:\"%s\" -osmi -ofpt -xfMACCS | grep -o '[0-9]\\.[0-9]\\?[0-9]'",
                sigmatches$nonsterosmiles[idx],sigmatches$prism_similar_smile[idx]
  )
  system(cmd,intern = T)
},mc.cores = 60) %>% unlist()

sigmatches$prism_similar_score=as.numeric(sigmatches$prism_similar_score)

write_tsv(sigmatches,"data/QTOF/db_screens/prism/all_sigmatches_uniq.tsv")


#write json for rban
sigmatches_for_rban  = sigmatches %>% 
                               filter(prism_similar_score>0.85) %>% 
                               mutate(json_id_np = paste(Name,prism_similar_name,"NP",sep = '__'),
                                      json_id_pr = paste(Name,prism_similar_name,"PR",sep = '__'))

bind_rows(select(sigmatches_for_rban,
         c(id=json_id_np,smiles=nonsterosmiles)),
select(sigmatches_for_rban,
       c(id=json_id_pr,smiles=prism_similar_smile))
) %>% jsonlite::toJSON() %>% write(file="rban_input_prismsimilar_np.json")

# similarity scores histogram
svg(filename = "figures/prism_np_similarity_hist.svg")
hist(sigmatches$prism_similar_score,
     main="Histogram of Tanimoto coefficients \n betweeen detected known compound \nand most similar predicted structure",
     xlab = "Tanimoto coefficients",
     col = "brown")
dev.off()


p_shared_comounds_freq = sigmatches %>% add_count(Name) %>% distinct(Name,.keep_all=T) %>% ggplot() + geom_histogram(aes(x=n),binwidth = 1,color="black",fill="brown") + ggprism::theme_prism() + labs(x="Number of strains sharing same compound", y="Frequency")
ggsave(p_shared_comounds_freq,filename = "figures/QTOF/shared_compounds_freq.svg")
p_npatlas_compounds_per_strain = sigmatches %>% group_by(strain) %>% summarize(number_of_compounds=n()) %>% arrange(number_of_compounds) %>%  ggplot() + geom_bar(aes(x=fct_reorder(strain,number_of_compounds),y=number_of_compounds),stat="identity",fill="brown",color="black") + ggprism::theme_prism() + labs(x="Strain",y="Number of known compounds") + theme(axis.text.y = element_text(size=8)) + coord_flip()
ggsave(p_npatlas_compounds_per_strain,filename = "figures/QTOF/npatlas_compounds_per_strain.svg")
