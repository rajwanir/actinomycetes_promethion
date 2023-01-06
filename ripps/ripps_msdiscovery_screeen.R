# getting compound elicitor description for compound hit
msdiscovery = read_csv("../actinomycetes_expression/data/msdiscoverLibrary/NP210401-nih.csv") %>% janitor::clean_names()
msdiscovery = msdiscovery %>% mutate(file=paste0(paste(plate,position,sep="_"),
                                                 ".d") )


hits = read_csv("data/ripps/GA3008/GA3008_msdiscovery_hits.csv")
hits = filter(hits,locus_tag=="ctg1_2547_lanth8_miss2") #ctg1_6461_rs16re23_fl0
hits = left_join(hits,msdiscovery,by="file")
hits = hits %>% arrange(-height) %>% distinct(cas_number,.keep_all=T)


hits$plate = str_remove_all(str_extract(hits$file,"-\\d+_"),"-|_")
hits_summary = hits %>% group_by(massprismdb,locus_tag,biosyn_class,core_seq) %>%
  summarise(n_plate=length(unique(plate)),n_inducers=length(unique(file)),
            height=max(height),rt_avg=mean(rt),
            mz_2=unique(mz_2),mz_3=unique(mz_3)) 
write_csv(hits_summary,"data/ripps/GA6-010/msdiscovery_hits_summary.csv")


ms_collected = read_csv("data/ripps/GA6-010/msdiscovery_hits_msmscollected.csv",comment = "#")$Precursor
hits_summary$msms =  sapply(1:nrow(hits_summary),function(idx){
  any(sapply(ms_collected, function(m)
    (m>get_ppm_window(hits_summary$mz_2[idx])[1] & m < get_ppm_window(hits_summary$mz_2[idx])[2]) |  (m>get_ppm_window(hits_summary$mz_3[idx])[1] & m < get_ppm_window(hits_summary$mz_3[idx])[2])
  ))})

ms_attempt = read_csv("data/ripps/GA6-010/msms_attempt.csv")$mz
hits_summary$msms_attempt = sapply(1:nrow(hits_summary),function(idx){
  any(sapply(ms_attempt, function(m)
    (m>get_ppm_window(hits_summary$mz_2[idx])[1] & m < get_ppm_window(hits_summary$mz_2[idx])[2]) |  (m>get_ppm_window(hits_summary$mz_3[idx])[1] & m < get_ppm_window(hits_summary$mz_3[idx])[2])
  ))})

write_csv(hits_summary,"data/ripps/GA6-010/msdiscovery_hits_summary.csv")




cids = lapply(unique(hits$cas_number), function(cas) webchem::get_cid(query = cas,domain = "compound",match="first",from = "xref/RN")) %>%
  bind_rows() %>% filter(!is.na(cid))

hits = left_join(hits,cids,by=c("cas_number"="query")) %>% filter(!is.na(cid))

#getting strcutures for hits
structures = lapply(hits$cid,function(cid){
  Sys.sleep(1)
  img = EBImage::readImage(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/PNG?record_type=2d&image_size=400x400",cid),
                           type = "png",names = cid)
  # return(readLines(sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/IsomericSMILES/txt",cid)))
  return(img)})

EBImage::display(EBImage::combine(structures),all=T,method="raster")

# sapply(1:6,function(y){
#   ymod = ((y-1)*300) + 20
#   sapply(1:6,function(x){
#    text(x = 300*(x-1), y = ymod, label = x*y, adj = c(0,1), col = "black", cex = 1)
#     #sprintf("%d %d",300*(x-1),(y-1)*300+20)
# })
# })

svg("abc.svg",height=8,width = 16)
par(mfrow=c(3,6))
par(mar = rep(0,4))
for(i in 1:length(structures)){
  EBImage::display(structures[[i]],method="raster",margin=0,spacing=0) 
  text(x=20,y=20,label = i, adj = c(0,1), col = "black", cex = 3)
}
dev.off()


barplot(hits$height, names.arg=as.character(1:nrow(hits)),ylab="Peak intensity",xlab="Coumpound #",col="burlywood1",
        main="Target compound peak intensity\n under elicitors",horiz=TRUE) 




hits = mutate(hits,
              locus=str_extract(locus_tag,"ctg1_\\d+|^K[A-Z]+_\\d+"),
              ms_ms_collected=if_else(!is.na(ms_ms_count),"Yes","No") )

hits = hits %>% mutate(target_peptide=if_else(locus_tag=="ctg1_2547_lanth8_miss2","Yes","No"))

filter(hits,str_detect(locus_tag,"ctg1_6461") ) %>% distinct(massprismdb,file,.keep_all=T) %>% 
ggplot(aes(x=rt,y=height,fill=target_peptide,label=core_seq)) +
  geom_point(aes(size=massagilent),shape=21,color="black") +  
  #ggrepel::geom_label_repel(xlim=c(8,12),color="grey") +
  scale_x_continuous(breaks = 2:8,limits=c(2,8)) + scale_y_continuous(trans='log10') + expand_limits(x = 2, y = 0) +
  scale_fill_manual(values=c("Yes"="lightblue","No"="white")) +
  labs(x="Retention time (min)",y="Peak height") +
#  facet_wrap(~paste(locus,biosyn_class),ncol=1,scales = "free_y") +
  theme_classic() + 
  theme(text=element_text(size=14),axis.text = element_text(size=14,color="black"),
        legend.position = "top") +
 coord_cartesian(xlim = c(2, 12),clip = "off")


#k peptide many matches
imp_locus_tag = filter(hits,str_detect(locus_tag,"^K") ) %>% distinct(massprismdb,file,.keep_all=T) %>% count(locus_tag) %>% filter(n>2) %>% select(locus_tag) %>% unlist()
filter(hits,str_detect(locus_tag,"^K") ) %>% filter(locus_tag %in% imp_locus_tag) %>% distinct(massprismdb,file,.keep_all=T) %>% 
  ggplot(aes(x=rt,y=height,fill=n_lanth,label=core_seq)) +
  geom_point(aes(fill=paste(core_seq,n_lanth)),shape=21,,color="black",size=4) +  
  #ggrepel::geom_label_repel(xlim=c(8,12),color="grey") +
  scale_x_continuous(breaks = seq(2,5,by=1),limits=c(2,5)) + scale_y_continuous(trans='log10') + expand_limits(x = 2, y = 0) +
  scale_fill_brewer(palette = "Set2") +
  labs(x="Retention time (min)",y="Peak height") +
  theme_classic() + 
  theme(text=element_text(size=14),axis.text = element_text(size=14,color="black"),
        legend.position = "none")  + facet_wrap(~paste(core_seq,n_lanth),nrow=1,strip.position = "top")

# repetitive mass
# should subtract masses from different strian first
agilent_features = read_csv("data/ripps/GA13-001/allfeatures_filtered.csv")
mass_filters = read_csv("data/ripps/GA13-001/additonal_mass_filters.txt")$mass_filters
agilent_features$mass_trim = floor(agilent_features$mass*100)/100
agilent_features=agilent_features %>% group_by(mass_trim) %>% mutate(rt_diff = max(rt)-min(rt), n=length(unique(file))) %>% ungroup()
highlight_mass=agilent_features$mass_trim[agilent_features$height>3e+04]
ggplot(agilent_features %>% filter(mass_trim>800 & n>1 & !is.na(max_z) &rt<10 & rt_diff<0.4),aes(x=rt,y=height)) +
  geom_point(aes(size=mass_trim,fill=mass_trim),shape=21,color="black",alpha=0.6) +
  ggrepel::geom_text_repel(aes(label=if_else(mass_trim%in%highlight_mass,mass_trim,NA_real_)))+
  #scale_x_continuous(breaks = seq(2,10,by=0.3),limits=c(2,10)) + 
  scale_y_continuous(trans='log10')+
  labs(x="RT (min)", y="Peak height")+
  ggprism::theme_prism() +
  theme(panel.grid=element_line(color = "grey",
                                size = 0.75,
                                linetype = 2))

hits = read_csv("data/ripps/GA3008/betacarotine_lanthi.csv") %>% janitor::clean_names()
ggplot(hits,aes(x=fct_reorder(as.character(conc),conc),y=height,group=conc,fill=conc)) +
  geom_boxplot(width=0.3) + geom_point(position = "jitter",size=3) + 
 # scale_x_continuous(breaks=c(0,9,18,36,72)) +  
  scale_y_continuous(trans='log10') +
  scale_fill_gradient(low="yellow", high="red") + ggprism::theme_prism() + 
  labs(x="Concentration (uM)", y="Peak height") + theme(legend.position = "none", panel.grid.major = element_line(color="grey",linewidth = 0.1))