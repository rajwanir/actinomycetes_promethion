
library(tidyverse)
library(Biostrings)

rippids = bind_rows(read_csv("../databases/MiBiG/ripp_val/ripps_id.csv"),read_csv("../databases/MiBiG/ripp_val/ripps_id_correctclassriip.csv"))
class="Lanthipeptide"
rippids = filter(rippids,biosyn_class %in% class)$mibig_id

mibigmass = read_tsv("../databases/MiBiG/ripp_val/mibig_compounds_edited.tsv")

# filter by pubchem id
#mibigmass=filter(mibigmass,mibig_id %in% rippids & !is.na(pubchem_id)) %>% select(mibig_id,pubchem_id,compound)
# mibigmass = bind_cols(mibigmass,
#                       bind_rows(lapply(mibigmass$pubchem_id,function(x) 
#                         read_csv(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%d/property/MolecularFormula,IsomericSMILES,MonoisotopicMass/CSV",x))))
# ) %>% mutate(MonoisotopicMass = round(MonoisotopicMass,2))

# filter by mol formula
mibigmass = filter(mibigmass,mibig_id %in% rippids & !is.na(molecular_formula)) %>% mutate(molecular_formula = str_remove(molecular_formula,pattern = "\\+"))
mibigmass$mol_mass = sapply(mibigmass$molecular_formula,function(x) rcdk::get.formula(x,charge=0)@mass)
mibigmass = mutate(mibigmass,mol_mass = round(mol_mass,1))


           

# get precursors
mibigprots = Biostrings::readAAStringSet("../databases/MiBiG/mibig_all_prot_seqs/mibig_prot_seqs_2.0.fasta") 
mibigprots = mibigprots[str_detect(names(mibigprots),paste(mibigmass$mibig_id,collapse = "|")) &
                          width(mibigprots) < 200]
precursors = tibble(name=names(mibigprots),sequence=as.character(mibigprots))


#match
precursors =precursors %>% mutate(mass_modified = round(mass_modified,1))
mibigmass %>% filter(!mol_mass %in% precursors$mass_modified)  


mibigmass2 = left_join(mibigmass,mutate(precursors,mibig_id = str_extract(name,"BGC\\d+")),
                       by=c("mol_mass"="mass_modified","mibig_id"="mibig_id")) %>% 
  mutate(correct_prediction = if_else(!is.na(core_seq),T,F)) %>% 
  arrange(-correct_prediction) %>% 
  distinct(mibig_id,compound,.keep_all=T)



# pie chart
correct_predict_summary = mibigmass2 %>% count(correct_prediction)
correct_predict_summary = correct_predict_summary %>%  arrange(desc(correct_prediction)) %>%
  mutate(prop = n / sum(correct_predict_summary$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

p_pie = ggplot(correct_predict_summary, aes(x="",y=prop,fill=correct_prediction)) + 
         geom_bar(width=100,stat="identity",color="black") +
         coord_polar("y",start=0) +theme_void() +
  geom_text(aes(y =ypos, label = sprintf("%d (%.1f %%)",n,prop)), color = "white", size=6) +
  scale_fill_manual(values=c("TRUE"="darkgreen","FALSE"="grey")) + 
  ggtitle(class) +
  theme(legend.position = "bottom",plot.title = element_text(hjust=0.5)) +
  guides(fill = guide_legend(title = "Reported mass predicted",title.position = "top",title.hjust = 0.5))
  

ggsave(p_pie,filename = sprintf("figures/ripps/mibig_val/pie_%s.svg",class))


## checking antismash class
rippids = rippids %>% filter(biosyn_class=="RiPP")
rippids = mutate(rippids,
                 mibiglink=sprintf("https://mibig.secondarymetabolites.org/repository/%s/generated/%s.1.region001.gbk",
                                   mibig_id,mibig_id))
query_class="lanthipeptide"
query_class = sprintf("product=\"%s",query_class)
rippids$lanthipeptide = sapply(rippids$mibiglink, function(link) any(grepl(query_class,readLines(link,n=200))))

rippids = mutate(rippids,biosyn_class = case_when(lanthipeptide==T ~ "Lanthipeptide",
                                                  lassopeptide==T ~ "Lassopeptide",
                                                  LAP==T ~ "LAP",
                                                  T ~ "RiPP"))


write_csv(select(rippids,mibig_id,main_product,biosyn_class),
          "../databases/MiBiG/ripp_val/ripps_id_correctclassriip.csv") 