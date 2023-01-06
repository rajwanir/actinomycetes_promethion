library(tidyverse)

#find 'data/assemblies/homopolish/' -name *_homopolished.fasta | parallel -j32 /home/rajwanir2/assembly-stats/build/assembly-stats -t {} 
assm_stats = read_tsv("data/assemblies/assm_stats.tsv")


# get coverage #####
coverage = Sys.glob(file.path("data/assemblies/canu/*/*.contigs.layout.tigInfo"))
names(coverage) = coverage
coverage = lapply(coverage,data.table::fread)
coverage = bind_rows(coverage,.id="tig_info_path")
coverage = coverage %>% filter(tigClass != "unassm")
coverage = coverage %>% group_by(tig_info_path) %>% summarize(coverage = mean(coverage)) 
coverage = coverage %>% mutate(strain = str_extract(tig_info_path,pattern = "\\w+-\\w+")) %>% select(-tig_info_path)

assm_stats = left_join(assm_stats,coverage,by="strain")
assm_stats = assm_stats %>% rename(assembly_length = total_length,contigs=number)
assm_stats = assm_stats %>% arrange(coverage)


# get antismash BGCs ####
gff = Sys.glob(file.path(getwd(), "data/antismash/*/*.gff.1"))
gff = lapply(gff,function(x) rtracklayer::readGFF(x)) %>% setNames(.,gff)
gff = bind_rows(gff,.id="gff_path")
gff = mutate(gff,strain=str_extract(gff_path,pattern = "\\w+-\\w+"))
n_bgcs = gff %>% group_by(strain,seqid) %>%
  filter(type...3 == "sequence_feature") %>% 
  filter(!is.na(protocluster_number)) %>%  
  summarize(bgc_per_contig = max(as.numeric(protocluster_number) )) %>% 
  group_by(strain) %>% 
  summarize(n_bgcs = sum(bgc_per_contig))


assm_stats = left_join(assm_stats,n_bgcs,by="strain")
assm_stats = rename(assm_stats,n_bgcs_antismash=n_bgcs)


# get prism BGCs #####
prism_jsons = Sys.glob("Y:/actinomycetes_promethion/data/prism/*.json")

prism_bgc_freq = lapply(prism_jsons,function(x) as.data.frame(table(unlist(jsonlite::fromJSON(x,flatten = T,simplifyDataFrame =T)$prism_results$clusters$type))))
names(prism_bgc_freq) = prism_jsons

prism_bgc_freq = bind_rows(prism_bgc_freq,.id = "prism_file")
prism_bgc_freq$strain = str_remove(str_extract(prism_bgc_freq$prism_file,pattern = "\\w+-\\w+_"),"_")
n_bgcs_prism = prism_bgc_freq %>% group_by(strain) %>% summarize(n_bgcs_prism = sum(Freq))

assm_stats = left_join(assm_stats,n_bgcs_prism,by="strain")




write_tsv(assm_stats,file = "data/assemblies/assm_stats_detailed.tsv")



p_assm_stats = ggplot(assm_stats) + theme_classic() + labs(x=NULL)
p_assm_stats = list(
p_assm_stats + geom_bar(aes(x=strain,y=assembly_length),stat = "identity") + coord_flip() ,
p_assm_stats + geom_bar(aes(x=strain,y=contigs),stat = "identity") + coord_flip(),
p_assm_stats + geom_bar(aes(x=strain,y=coverage),stat = "identity") + coord_flip(),
p_assm_stats + geom_bar(aes(x=strain,y=n_bgcs_antismash),stat = "identity") + coord_flip(),
p_assm_stats + geom_bar(aes(x=strain,y=n_bgcs_prism),stat = "identity") + coord_flip()
)

p_assm_stats = cowplot::plot_grid(plotlist = p_assm_stats,nrow=1) 
ggsave(p_assm_stats,filename = "figures/p_assmstats.svg",width = 13)  


bgcs_type = gff %>% filter(type...3 == "sequence_feature") %>% 
  filter(!is.na(protocluster_number)) %>% 
  group_by(strain,product) %>% distinct(protocluster_number) %>% 
  summarize(n_bgcs_per_type = n())

p_bgcs_type = ggplot(bgcs_type) +
  geom_tile(aes(x=product,y=strain,fill=n_bgcs_per_type),color="black") +
  geom_text(aes(x=product,y=strain,label=n_bgcs_per_type)) +
  scale_fill_distiller(palette="RdPu",direction = 1) +
  labs(y=NULL, x = "BGC type") +
  guides(fill=guide_colourbar(title = "# of BGCs")) +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2),
        axis.title.x = element_text(face = "bold"),
        legend.position = "top",
        plot.margin = margin(t=0, unit = "pt")
  ) 

ggsave(p_bgcs_type,filename = "figures/p_bgcs_type.svg",width = 11,height = 9)  


p_bgcs_type_prism = ggplot(prism_bgc_freq) +
  geom_tile(aes(x=Var1,y=strain,fill=Freq),color="black") +
  geom_text(aes(x=Var1,y=strain,label=Freq)) +
  scale_fill_distiller(palette="RdPu",direction = 1) +
  labs(y=NULL, x = "BGC type") +
  guides(fill=guide_colourbar(title = "# of BGCs")) +
  theme_linedraw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2),
        axis.title.x = element_text(face = "bold"),
        legend.position = "top",
        plot.margin = margin(t=0, unit = "pt")
  )


ggsave(p_bgcs_type_prism,filename = "figures/p_bgcs_type_prism.svg",width = 11,height = 9)  

