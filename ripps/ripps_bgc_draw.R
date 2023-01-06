library(tidyverse)
library(gggenes)

hits = read_tsv("data/ripps/hits.tsv")
msms_pass= read_csv("data/ripps/msms_pass.csv")
hits=left_join(hits,msms_pass) %>% filter(msms_pass==T)
colors = read_csv("data/ripps/bgc_draw_colors.csv")

#filter hits
# hits = filter(hits,str_detect(BGC_class,"lanthipeptide") & !str_detect(name,pattern = "GA7-009"))
#hits = hits %>% separate(name,into=c("strain","locus_tag"),sep="@",remove=F)
hits = hits %>% mutate(strain=strainprismdb,
                       gff_path = sprintf("data/antismash/%s/%s_homopolished.gff.1",strain,strain))
hits = hits %>% mutate(precursor = str_extract(locus_tag,"ctg\\d+_\\d+"))

hits=hits %>% distinct(strain,precursor,.keep_all=T)
# 



get_ripp_bgc = function(strain="x",precursor="y",window=10000, return_region = T){
  gff_path = sprintf("data/antismash/%s/%s_homopolished.gff.1",strain,strain)
  gff=rtracklayer::readGFF(gff_path,filter = list("type"="CDS")) 
  precursor_start=filter(gff,locus_tag==precursor)$start
  precursor_seqid=filter(gff,locus_tag==precursor)$seqid
  bgc_region=rtracklayer::readGFF(gff_path,filter = list("type"="sequence_feature"),tags = c("product","protocluster_number","contig_edge")) %>% 
    filter(!is.na(protocluster_number) & !is.na(contig_edge)) %>% 
    filter(str_detect(product,"lanthi")) %>%
    filter(start< precursor_start & end > precursor_start & 
             seqid == as.character(precursor_seqid)) %>% 
    head(1) # strain GA9-001 has two overlapping lanth BGC
  
  if(return_region==T){return(bgc_region)}
  
  bgc_class = select(bgc_region,product) %>% as.character()
  
  
  gff=filter(gff,str_detect(locus_tag ,pattern=str_extract(precursor,pattern = "ctg\\d+")) &
              start> bgc_region$start & 
               end < bgc_region$end &
               type== "CDS")
  gff=mutate(gff,gene=case_when(locus_tag==precursor ~ "precursor",
                                #str_detect(gene_kind,"transport") ~ "transport",
                               # str_detect(gene_kind,"regulatory") ~ "regulatory",
                                str_detect(sec_met_domain,"Peptidase") ~ "peptidase",
                                str_detect(sec_met_domain,"LANC") ~ "lanC (cyclase)",
                                str_detect(sec_met_domain,"LanB") ~ "lanB (dehydratase)",
                                str_detect(sec_met_domain,"APH") ~ "APH-like (phosphotransferase)",
                                str_detect(sec_met_domain,"micKC ") ~ "micKC (kinase-cyclase)",
                                str_detect(sec_met_domain,"HopA1") ~ "HopA1 (dehydratase)",
                                str_detect(gene_functions,"oxidoreductase") ~ "oxidoreductase",
                                str_detect(gene_functions,"methyltransferase") ~ "methyltransferase",
                                str_detect(gene_functions,"cytochrome") ~ "cytochrome p450"
                               # str_detect(gene_functions,"SMCOG1030") ~ "Ser/thr protein kinase"
                               ),
             bgc_class=bgc_class)
return(gff)
}

gffs=lapply(1:nrow(hits), function(idx) get_ripp_bgc(strain=hits$strain[idx],precursor = hits$precursor[idx]))
names(gffs) = paste(hits$strain,hits$precursor,sep = "@")
gffs=dplyr::bind_rows(gffs,.id="peptide_id")
gffs =mutate(gffs,peptide_id=paste(gffs$peptide_id,
                                 gffs$bgc_class,sep = "\n"))

gffs=mutate(gffs,peptide_id=reorder(peptide_id,bgc_class))

gffs$peptide_id <- factor(gffs$peptide_id, levels = unique(gffs[order(gffs$bgc_class), "peptide_id"]))



p_bgcs = ggplot(gffs) +
  geom_gene_arrow(aes(xmin=start,xmax=end,y=peptide_id,fill=gene),size=0.2,color="black") + 
  facet_wrap(~peptide_id,scales="free",ncol=1)  +
  theme_void() +
  theme(axis.text.y = element_text(size=10),
        strip.text = element_blank(),
        panel.grid.major.y = element_line(color="black",size=0.5),
        legend.position = "top") +
  guides(fill=guide_legend(nrow=2,title.position = "top",title.hjust = 0.5,
                           label.theme=element_text(size=10), title.theme = element_text(size=10))) +
  scale_fill_manual(values=deframe(colors),na.value = "grey86") +
  scale_color_manual(values=deframe(colors),na.value = "grey86")



ggsave(p_bgcs,filename = "figures/ripps/lanthi_bgcs.svg", height = 5, width = 20)


