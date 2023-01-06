


agilent_features = bind_rows(parallel::mclapply(Sys.glob("../actinomycetes_expression/data/msdiscoverLibrary/features/technical_reproducibility/*.csv"),function(x) read_csv(x,comment = "#",skip = 2) %>% janitor::clean_names(),mc.cores=1)) %>% mutate(strain="gresieus")
agilent_features = mutate(agilent_features,rep = if_else(str_detect(file,"tech"),"rep2","rep1"))
replicate1 = agilent_features %>% filter(rep=="rep1")
replicate2 = agilent_features %>% filter(rep=="rep2")

get_ppm_window = function(mass_target){
  mass_window = c( mass_target - (mass_target * 15 / 10^6),
                   mass_target + (mass_target * 15 / 10^6))
  return(mass_window)}


repmatch = lapply(replicate1$mass,function(m){
  replicate2 %>% filter(between(mass,get_ppm_window(m)[1],get_ppm_window(m)[2])) %>% head(1)
}) 

agilent_features$mass_trim = floor(agilent_features$mass*10)/10
agilent_features = mutate(agilent_features,cpd_id=paste0(mass_trim,"@",file),
                          rep = if_else(str_detect(file,"tech"),"rep2","rep1"))
replicate_summary = agilent_features %>% group_by(cpd_id,rep) %>% summarize(height=mean(height)) %>% ungroup() %>% mutate(cpd_id=str_remove(cpd_id,"techrep1_"))
replicate_summary = replicate_summary %>% pivot_wider(id_cols ="cpd_id",names_from ="rep",values_from ="height") %>% replace_na(list(rep1=0,rep2=0))


ggplot(replicate_summary %>% filter(rep1!=0&rep2!=0),aes(x=rep1,y=rep2)) + geom_point() +
  ggprism::theme_prism() + 
  labs(x="replicate-1",y="replicate-2") +
  ggtitle("Technical reproducibility (R2=0.903)") +
  theme(panel.grid.minor = element_line(color="grey"),
        title = element_text(hjust=0.5))