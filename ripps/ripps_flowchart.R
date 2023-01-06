library(tidyverse)
library(cowplot)

## genomics workflow ############
gff = Sys.glob(file.path(getwd(), "data/antismash/*/*.gff.1"))


regions = lapply(gff,function(x) rtracklayer::readGFF(x,tags=c("product","candidate_cluster_number","protocluster_number","contig_edge"),
                                                      filter = list("type"="sequence_feature"))) %>% 
               setNames(.,gff)
regions = bind_rows(regions,.id="gff_path")
regions = mutate(regions,strain=str_extract(gff_path,pattern = "\\w+-\\w+"))
regions = filter(regions,
       is.na(candidate_cluster_number),!is.na(contig_edge),
       product %in%c("lanthipeptide-class-i","lanthipeptide-class-ii",
                   "lanthipeptide-class-iii","lanthipeptide-class-iv",
                   "lanthipeptide-class-v"))
regions = mutate(regions,across(where(is.factor), as.character))


p_regions = ggplot(regions,aes(x=str_remove_all(product,"lanthipeptide-class-"))) + 
  geom_bar(stat="count",fill="brown",color="black") + labs(x="lanthipeptide class") +
  geom_text(stat="count",aes(label=..count..),size=6) +
  ggprism::theme_prism(base_fontface = "plain") +
  ggtitle(sprintf("# lanthipeptide BGCs : %d",nrow(regions)),
          subtitle = sprintf("# Strains : %d", length(unique(regions$strain))) )



gff=lapply(gff,function(x) rtracklayer::readGFF(x,tags=c("core_sequence","locus_tag","translation"),
                                            filter = list("type"="CDS"))) %>% 
  setNames(.,gff) %>% bind_rows(.id="gff_path") %>% 
  mutate(strain=str_extract(gff_path,pattern = "\\w+-\\w+"),
         width = abs(start-end))
precursorRegions=lapply(1:nrow(regions),function(idx)
  gff %>% filter(
    is.na(core_sequence) &
      gff_path == regions$gff_path[idx] &
      seqid == regions$seqid[idx] &
      start > regions$start[idx] &
      end < regions$end[idx] & 
      type == "CDS" &
      width < 600
    
  )
)
rm(gff)

precursors = select(bind_rows(precursorRegions),strain,locus_tag,translation,width,seqid) %>%
  mutate(name = paste(strain,locus_tag,sep = "@"),
         sequence = translation) %>% select(-c(strain,locus_tag,translation))


precursors = precursors %>% mutate(cterminus_seq = str_sub(sequence,start=-20),
                                   n_ser_thr = str_count(cterminus_seq,pattern = "S|T")) %>% 
  filter(n_ser_thr>2 & str_detect(cterminus_seq,pattern = "C")) %>% distinct()


rm(precursorRegions)

#getting bgcs for precursors
precursors = precursors %>% separate(name,into=c("strain","locus_tag"),sep="@",remove = F)
precursorBGCRegions = lapply(1:nrow(precursors), function(idx) get_ripp_bgc(strain=precursors$strain[idx],
                                                                            precursor = precursors$locus_tag[idx],
                                                                            return_region = T)) 
names(precursorBGCRegions) = precursors$name
precursorBGCRegions=bind_rows(precursorBGCRegions,.id="peptide_id") %>% separate(peptide_id,into=c("strain","locus_tag"),sep="@",remove=F) 
precursors=left_join(precursors,precursorBGCRegions)




p_precursorlength = ggplot(precursors,aes(x=width/3,y="")) + geom_boxplot(fill="brown") + geom_jitter() + labs(y="",x="precursor length") + 
  ggprism::theme_prism(base_fontface = "plain") + theme(axis.line.y = element_blank(),axis.ticks.y = element_blank()) +
  ggtitle(sprintf("# of precursors : %d",
                  nrow(precursors)),
          subtitle = sprintf("containing atleast two Ser/Thr \nand a Cys at C-terminus (last 20 residues)\n# BGCs : %d",
                             nrow(distinct(precursors,protocluster_number,seqid,strain)))
          ) + theme(plot.subtitle = element_text(size=12))

p_cterm_nserthr = ggplot(precursors,aes(x=as.character(n_ser_thr))) + geom_bar(fill="brown",stat="count",color="black") +
  labs(y="",x="# of Ser/Thr at C-terminus") + 
  ggprism::theme_prism(base_fontface = "plain") 


p_precursors = plot_grid(p_precursorlength,p_cterm_nserthr,ncol=1)


cores=sapply(5:30,function(i) str_sub(precursors$sequence,start=-i)) %>% as.vector()
cores=tibble(n_ser_thr = str_count(cores,"S|T"),
             coreseq=cores,
             nW=str_detect(cores,"W"),
             naviC = str_detect(cores,"C$"),
             nasp = str_detect(cores,pattern = "D|B"),
             nlactyl = str_detect(cores,pattern = "^S|^T"),
             nobu= str_detect(cores,pattern = "^T"),
             npyr =str_detect(cores,pattern = "^S")) %>% 
  filter(n_ser_thr>0 & str_detect(coreseq,"C"))

p_cores =  Gmisc::boxGrob(sprintf("# of possible cores (residues 5-30 from C-terminus)\n %d",
                                  nrow(cores))
               )


n_modified_peptides = sum(nrow(cores) + sum(cores$n_ser_thr-1))
p_modification =  Gmisc::boxGrob(sprintf("# of possible modified peptides with atleast\none dehydration and upto n-1 missed sites\n %d",
                                         n_modified_peptides)
)

n_w_cores = sum(cores[cores$nW==T,]$nW) + sum(cores[cores$nW==T,]$n_ser_thr-1)
n_asp_cores = sum(cores[cores$nasp==T,]$nasp) + sum(cores[cores$nasp==T,]$n_ser_thr-1)
n_avic_cores = sum(cores[cores$naviC==T,]$naviC) + sum(cores[cores$naviC==T,]$n_ser_thr-1)
n_lactyl_cores = sum(cores[cores$nlactyl==T,]$nlactyl) + sum(cores[cores$nlactyl==T,]$n_ser_thr-1)
n_obu_cores = sum(cores[cores$nobu==T,]$nobu) + sum(cores[cores$nobu==T,]$n_ser_thr-1)
n_pyr_cores = sum(cores[cores$npyr==T,]$npyr) + sum(cores[cores$npyr==T,]$n_ser_thr-1)

# p_tailor = Gmisc::boxGrob(sprintf("Additional tailoring \n
#                                   Acetylation=%d
#                                   TrpCl=%d
#                                   AspOH=%d
#                                   aviCys=%d
#                                   total modified peptides=%d",
#                                   n_modified_peptides,
#                                   n_w_cores,
#                                   n_asp_cores,
#                                   n_avic_cores,
#                                   (n_modified_peptides*2) + n_w_cores + n_asp_cores +  n_avic_cores)
# )



tailoring_counts = tibble(modification = c("N-Acetyl","TrpCl","AspOH","aviCys","Lactyl","Obu","Pyr"),
       count=c(n_modified_peptides,n_w_cores,n_asp_cores,n_avic_cores,
               n_lactyl_cores,n_obu_cores,n_pyr_cores))

p_tailor = ggplot(tailoring_counts,
       aes(x=modification,y=count)) +
  geom_bar(stat="identity",fill="brown",color="black") + labs(x="Tailoring modification") +
  geom_text(stat="identity",aes(label=count),size=5) +
  ggprism::theme_prism(base_fontface = "plain")+ theme(axis.text.x = element_text(angle = 320))


p_totalpeptides =  Gmisc::boxGrob(sprintf("# total modified peptides : %d",
                                  (n_modified_peptides*2) + 
                                    n_w_cores + n_asp_cores +  n_avic_cores +
                                    n_lactyl_cores+n_obu_cores+n_pyr_cores)
)


p_workflow_genomics = plot_grid(p_regions,p_precursors,
          plot_grid(p_cores,p_modification,ncol=1),
          plot_grid(p_totalpeptides,p_tailor,ncol=1),
          nrow = 1)


ggsave(p_workflow_genomics,filename = "figures/ripps/workflow_genomics.svg",width = 16,height=5)



## metabolomics workflow ####

agilent_features = lapply(Sys.glob("data/QTOF/db_screens/agilent_features/*_peptides.csv"),read_csv) %>% bind_rows() %>% 
  janitor::clean_names() %>% 
mutate(strain = str_extract(file,pattern = "G\\w+-\\d+")) %>% 
  filter(strain %in% unique(regions$strain))


media = read_csv("data/QTOF/db_screens/agilent_features/media/ISP1_MEOH_MC_peptides.csv") %>% 
  janitor::clean_names() 

agilent_features$media_pass = sapply(1:nrow(agilent_features),filter_media)

p_totalfeatures_mass = ggplot(agilent_features,aes(mass)) + geom_histogram(color="black",fill="darkgreen") + 
  scale_x_continuous(breaks = seq(0,max(agilent_features$mass),by=500)) +
  ggprism::theme_prism(base_fontface = "plain") +
  theme(axis.line.y = element_blank(),axis.text.y=element_blank(),axis.ticks = element_blank()) +
  labs(y=NULL)+
  ggtitle(sprintf("# mass features: %d",nrow(agilent_features)))

p_totalfeatures_rt=ggplot(agilent_features,aes(rt)) + geom_histogram(color="black",fill="darkgreen") + 
  scale_x_continuous(breaks = seq(0,max(agilent_features$rt),by=1)) +
  ggprism::theme_prism(base_fontface = "plain") +
  theme(axis.line.y = element_blank(),axis.text.y=element_blank(),axis.ticks = element_blank()) +
  labs(y=NULL,x="retention time (min)")

p_totalfeatures_int=ggplot(agilent_features,aes(height)) + geom_histogram(color="black",fill="darkgreen") + 
  scale_x_continuous(trans='log10') +
  ggprism::theme_prism(base_fontface = "plain") +
  theme(axis.line.y = element_blank(),axis.text.y=element_blank(),axis.ticks = element_blank()) +
  labs(y=NULL,x="peak height")

p_totalfeatures = plot_grid(p_totalfeatures_mass,p_totalfeatures_rt,p_totalfeatures_int,ncol = 1)

#media
p_media_mass = ggplot(agilent_features,aes(mass)) + geom_histogram(color="black",aes(fill=media_pass)) + 
  scale_x_continuous(breaks = seq(0,max(agilent_features$mass),by=500)) +
  ggprism::theme_prism(base_fontface = "plain") +
  theme(axis.line.y = element_blank(),axis.text.y=element_blank(),axis.ticks = element_blank(),legend.position = "none") +
  labs(y=NULL)+
  ggtitle(sprintf("# features in media: %d",nrow(agilent_features)-sum(agilent_features$media_pass))) +
  scale_fill_manual(values=c("TRUE"="darkgreen","FALSE"="grey"))

p_media_rt=ggplot(agilent_features,aes(rt)) + geom_histogram(color="black",aes(fill=media_pass)) + 
  scale_x_continuous(breaks = seq(0,max(agilent_features$rt),by=1)) +
  ggprism::theme_prism(base_fontface = "plain") +
  theme(axis.line.y = element_blank(),axis.text.y=element_blank(),axis.ticks = element_blank(),legend.position = "none") +
  labs(y=NULL,x="retention time (min)")+
  scale_fill_manual(values=c("TRUE"="darkgreen","FALSE"="grey"))

p_media_int=ggplot(agilent_features,aes(height)) + geom_histogram(color="black",aes(fill=media_pass)) + 
  scale_x_continuous(trans='log10') +
  ggprism::theme_prism(base_fontface = "plain") +
  theme(axis.line.y = element_blank(),axis.text.y=element_blank(),axis.ticks = element_blank(),legend.position = "none") +
  labs(y=NULL,x="peak height")+
  scale_fill_manual(values=c("TRUE"="darkgreen","FALSE"="grey"))

p_media = plot_grid(p_media_mass,p_media_rt,p_media_int,ncol = 1)


#hits
hits = read_tsv("data/ripps/hits.tsv")
hits = mutate(hits,precursor_id=str_extract(locus_tag,"(ctg\\d+_\\d+)|(allorf_\\d+_\\d+)"),
              gff_path = sprintf("data/antismash/%s/%s_homopolished.gff.1",strainprismdb,strainprismdb)) 


hits_bgc_regions=lapply(1:nrow(hits), function(idx) get_ripp_bgc(strain=hits$strainprismdb[idx],precursor = hits$precursor_id[idx],return_region = T))
names(hits_bgc_regions) = paste(hits$strainprismdb,hits$precursor_id,sep = "@")
hits_bgc_regions=dplyr::bind_rows(hits_bgc_regions,.id="peptide_id") %>% separate(peptide_id,into=c("strainprismdb","precursor_id"),sep="@") %>% distinct()
hits=left_join(hits,hits_bgc_regions)


p_hits_pepvariants = ggplot(hits %>% add_count(precursor_id,strainprismdb,name="n_variants") ) +
  geom_bar(aes(x=n_variants), stat="count",fill="darkgreen",color="black") +
  scale_y_continuous(limits = c(0,5)) +
  facet_wrap(~str_remove(product,"lanthipeptide-"),nrow=1,scales="free_y")+
ggtitle(sprintf("# peptide mass matches at MS1: %d",nrow(hits)),
        subtitle = sprintf("# precursors: %d\n# BGCs: %d\n# strains: %d",
                           distinct(hits,precursor_id,strainprismdb) %>% nrow(),
                           nrow(distinct(hits,seqid,strainprismdb,protocluster_number)),
                           length(unique(hits$strainprismdb)))) +
  ggprism::theme_prism(base_fontface = "plain") + labs(x="peptide variants from same precursor") 


p_hits_mass = ggplot(hits,aes(x=massprismdb,y=ppm)) + geom_point(color="darkgreen") +
  ggprism::theme_prism(base_fontface = "plain") + labs(x="mass") 



p_hits_ions = ggplot(hits,aes(x=ions,y=height)) + geom_point(color="darkgreen") +
  ggprism::theme_prism(base_fontface = "plain") + labs(x="ions",y="peak height") + 
  scale_x_continuous(breaks = 0:max(hits$ions)) + 
  scale_y_continuous(trans='log10') 


p_hits = plot_grid(p_hits_pepvariants,
                   plot_grid(p_hits_mass,p_hits_ions,nrow=1),
                   ncol=1)


# eic and ms/ms pass
p_eicpass =  Gmisc::boxGrob(sprintf("# eic pass : %d",
                                          sum(read_csv("data/ripps/eic_pass.csv")$eic_pass,na.rm = T))
)

p_msmspass =  Gmisc::boxGrob(sprintf("# MS/MS pass : %d",
                                    sum(read_csv("data/ripps/msms_pass.csv")$msms_pass,na.rm=T))
)




p_metablomicsworkflow = plot_grid(p_totalfeatures,p_media,p_hits,
          plot_grid(p_eicpass,p_msmspass,ncol=1),
          rel_widths = c(0.4,0.4,0.8,0.2),nrow=1)


ggsave(p_metablomicsworkflow,filename = "figures/ripps/workflow_metablomics.svg",width = 16,height=7)



p_workflow = plot_grid(p_workflow_genomics,p_metablomicsworkflow,
          ncol=1,
          rel_heights = c(0.8,1))


ggsave(p_workflow,filename = "figures/ripps/workflow.svg",width = 16,height=10)
