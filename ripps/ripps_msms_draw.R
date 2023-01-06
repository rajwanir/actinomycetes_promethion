library(tidyverse)
library(ggtext)




hits = read_csv("data/ripps/hits.csv")

#filter hits
hits = filter(hits,str_detect(BGC_class,"lanthipeptide") & !str_detect(name,pattern = "GA7-009"))
hits = hits %>% separate(name,into=c("strain","locus_tag"),sep="@",remove=F)
hits = hits %>% mutate(gff_path = sprintf("data/antismash/%s/%s_homopolished.gff.1",strain,strain))
hits = hits %>% mutate(precursor = str_extract(locus_tag,"ctg\\d+_\\d+"))
hits = hits %>% mutate(modified_peptide = str_replace_all(putative_peptide,pattern = fixed("S(Dehydrated)"),
                                                          replacement = "<span style='color:#00FF00;'>Dha</span>"),
                modified_peptide = str_replace_all(modified_peptide,pattern = fixed("T(Dehydrated)"),
                                                   replacement = "<span style='color:#0000FF;'>Dhb</span>"))
# 


msms_paths = Sys.glob("data/ripps/msms/*.csv")



draw_msms_ripps=function(idx=1){
msms = read_csv(msms_paths[str_detect(msms_paths,hits$name[idx])]) %>% janitor::clean_names()
peptide=hits$putative_peptide[idx]
modified_peptide=hits$modified_peptide[idx]
peptide_name=hits$name[idx]



msms = msms %>% mutate(internal_frag=case_when(!is.na(ion_type)&!is.na(index)~NA_character_,
                                        !is.na(ion_type)&is.na(index)~str_remove(ion_type,pattern = "-H2O|-NH3|-CO")),
                       adduct=str_extract(ion_type,pattern = "H20|NH3|CO"),
                       ion_class=case_when(!is.na(ion_type)&!is.na(index)~str_extract(ion_type,pattern = "a|b|c|x|y|z"),
                                           !is.na(ion_type)&is.na(index)~"i"))

  

msms = msms %>% mutate(internal_frag_start= str_locate(peptide,fixed(internal_frag))[,1],
                internal_frag_end= str_locate(peptide,fixed(internal_frag))[,2],
                ion_type2=case_when(is.na(internal_frag)&!is.na(index)~paste(ion_class,index,adduct,sep = "-") %>% str_remove("-NA"),
                                    !is.na(internal_frag)~paste(internal_frag_start,internal_frag_end,
                                                                adduct,
                                                                sep = "-") %>% str_remove("-NA")))

msms = msms %>% mutate(ion_label=case_when(is.na(internal_frag)~ion_type2,
                                           !is.na(internal_frag)~NA_character_))





p=ggplot(msms,aes(x=m_z,y=intensity)) +
  geom_bar(stat="identity",width=0.7,fill="black") +
#  scale_y_continuous(limits = c(0,max(msms$intensity,na.rm = T)+200))+
  coord_cartesian(clip = "off") +
  ggrepel::geom_text_repel(aes(label=ion_label,color=ion_class),
                           min.segment.length=0,show.legend=F,
                           ylim = c(max(msms$intensity,na.rm =T),NA),
                           segment.linetype="dashed",segment.size=0.6,fontface="plain") +
  theme_test() +
  theme(axis.line.x = element_line(color="black",size=1),
        axis.text.x = element_text(color="black"),
        axis.title.x = element_text(face="italic"),
        axis.line.y=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  labs(x="m/z") +
  scale_color_manual(values = c("x"="red","y"="red","z"="red",
                                "a"="lightblue","b"="lightblue","c"="lightblue")) +
  labs(title = paste(peptide_name,modified_peptide,sep="<br>")) + theme(plot.title = element_markdown(lineheight = 1.1,hjust =0.5))

return(p)}


p_msms_ripps = lapply(1:nrow(hits), function(idx) draw_msms_ripps(idx))



p_msms_ripps = cowplot::plot_grid(plotlist = p_msms_ripps,ncol=1)


p_lanthi = cowplot::plot_grid(p_lanthi_genes,NULL,p_msms_ripps,rel_heights = c(0.2,0.05,0.8),ncol=1,
                              labels = c("A","B"))
ggsave(p_lanthi,filename = "figures/ripps/lanthi_complete.svg",width = 11,height=10)