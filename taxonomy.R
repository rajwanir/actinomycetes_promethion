library(tidyverse)


#mash dist data/mash/refseq.genomes.k21s1000.msh -l data/mash/assm_list -p 60 > data/mash/mash.out 
mash = read_tsv("data/mash/mash.out",
                col_names = c("ref","query","ident","pval","sketches")) %>% 
  group_by(query) %>% 
  arrange(.,pval) %>% 
  slice(1)

write_tsv(mash,"data/mash/mash.out.best")


mash_best = read_tsv("data/mash/mash.out.best")
mash_best = mash_best %>% mutate(ref_acc=str_extract(ref,pattern = "G(.)*\\.[0-9]"))
write_tsv(mash_best,"data/mash/mash.out.best")



mash_best = mash_best %>% mutate(ani_cmds = sprintf("../soil_metagenomics/fastANI -t 8 -q %s -r data/mash/ref_genomes/%s* -o /dev/stdout | cut -f3",
                                        query,ref_acc))

mash_best$ani = sapply(mash_best$ani_cmds,function(x) as.numeric(system(x,intern = T)),
              USE.NAMES = F)

mash_best$ani[sapply(mash_best$ani, function(x) if_else(length(x)==0,T,F))] <- NA
mash_best$ani=mash_best$ani %>% unlist()

taxonomy = read_delim("data/mash/ref_genomes_taxon.txt",delim = ":",col_names = c("assm_path","organsim","organism_name"))
taxonomy$genus = str_extract(taxonomy$organism_name,pattern = "\\w+")
taxonomy$assm_acc = str_extract(taxonomy$assm_path,pattern = "G(.)*\\.[0-9]")

mash_best = left_join(mash_best,taxonomy %>% select(assm_acc,genus),by=c("ref_acc"="assm_acc"))
