library(tidyverse)
library(ggprism)

folder = "plate1_flongle_07212021"
avg_output=0.167

epi2meqc = read_csv(Sys.glob(sprintf("data/flongleQC/%s/epi2me_qc.csv",folder)))
epi2meqc = epi2meqc %>% select(read_id,seqlen,mean_qscore)

#throughput
throughput=read_csv(Sys.glob(sprintf("data/flongleQC/%s/throughput*.csv",folder))) %>% janitor::clean_names()

#guppy barcode summary
barcode_summary = read_tsv(Sys.glob(sprintf("data/flongleQC/%s/demultiplexed/barcoding_summary.txt",folder)
                                            )) %>% select(read_id,barcode_arrangement)
#guppy sequencing summary
sequencing_summary = read_tsv(Sys.glob(sprintf("data/flongleQC/%s/fastq/sequencing_summary*.txt",folder)
)) %>% select(read_id,sequence_length_template)

#wimp
wimp = read_csv(Sys.glob(sprintf("data/flongleQC/%s/wimp/*.csv",folder)
                                 ))



#merge epi2me qc and barcode
barcode_summary = merge(barcode_summary,epi2meqc,
                        by.x = c("read_id"), by.y = ("read_id"))



#merge guppy and wimp
barcode_summary = merge(barcode_summary,wimp,
      by.x = c("read_id"), by.y = ("readid"))

#extract genus
barcode_summary = mutate(barcode_summary,
       name = case_when(name == "unclassified Streptomyces" ~ "Streptomyces unclassified", TRUE ~ name),
       genus = str_extract(name,"\\w+"))


#collapse low frequency genra
genus_freq = barcode_summary %>% group_by(barcode_arrangement,genus) %>% count() %>%
 group_by(barcode_arrangement) %>% summarize(genus_perc = n/sum(n),genus)
barcode_summary = merge(barcode_summary,genus_freq,
      by= c("barcode_arrangement","genus"))
barcode_summary=mutate(barcode_summary,
                       genus = ifelse(genus_perc < 0.1, "Others (<10% individual)", genus))

barcode_summary = barcode_summary %>% mutate(actino = ifelse(str_detect(lineage,":201174:"),"Actinobacteria","non-Actinobacteria"))

# calculate actino frequency 
actino_freq = barcode_summary %>% group_by(barcode_arrangement,actino) %>% count() %>%
  group_by(barcode_arrangement) %>% summarize(actino_perc = n/sum(n),actino) 

barcode_summary = merge(barcode_summary,actino_freq,
                        by= c("barcode_arrangement","actino"))

# add barcode key
barcode_summary =  mutate(barcode_summary, 
                          barcode_n = str_extract(barcode_arrangement,"\\d+"),
                          barcode_numeric = as.numeric(barcode_n))
barcode_key = read_delim(sprintf("data/flongleQC/%s/barcode_key.csv",folder),delim = ",") %>% select(row,col,barcode,strain)
barcode_summary = merge(barcode_summary,barcode_key,
                        by.x="barcode_numeric",by.y = "barcode",
                        all.x = T, all.y = F)

#add nanodrop data and
barcode_summary = merge(barcode_summary,conc$comb,
      by.x="barcode_n",by.y = "barcode",
      all.x = T, all.y = F)


#to select samples
selected = barcode_summary %>% group_by(barcode_arrangement) %>%
 summarize(n_reads = n(),
          actino_perc = unique(actino_perc),
          actino = unique(actino),
          # x260_230 = unique(x260_230),
          # x260_280 = unique(x260_280),
          strain = unique(strain),
          row = unique(row),
          col = unique(col),
          bases = sum(seqlen)/1000000,
          median_length = median(seqlen),
          vol_barcode = avg_output/bases
         ) %>%
  filter(actino == "Actinobacteria" & actino_perc > 0.5 & n_reads >4) %>%
  arrange(-actino_perc,n_reads) 


#plot ####

#to sort barchart fct_reorder(barcode_arrangement,barcode_arrangement, function(x)-length(x)
p_nreads = ggplot(barcode_summary %>% filter(barcode_arrangement!="unclassified"),
       aes(x=str_extract(barcode_arrangement,'\\d+')
                         )) +
  geom_bar(stat="count", color = "black") +
  labs(x="Barcode",y="Number of reads") +
  scale_y_continuous(expand=c(0,0), labels = scales::comma,
                     #breaks = seq(0,80000,by = 5000)
                     ) +
  theme_prism() +
  viridis::scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text(size=6))


nbases =  barcode_summary %>% filter(barcode_arrangement!="unclassified") %>%
  group_by(barcode = str_extract(barcode_arrangement,'\\d+')) %>% 
  summarize(bases=sum(seqlen))
p_nbases = ggplot(nbases,
                  aes(x=barcode,
                      y=bases/1000000
                  )) +
  geom_bar(stat="identity", color = "black") +
  labs(x="Barcode",y="Megabases") +
  scale_y_continuous(expand=c(0,0)) +
  theme_prism() +
  viridis::scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text(size=6))





#fct_reorder(str_extract(barcode_arrangement,'\\d+'),actino_perc )
p_taxonomy = ggplot(barcode_summary %>% filter(barcode_arrangement!="unclassified"),
                  aes(x=str_extract(barcode_arrangement,'\\d+'),
                  fill = actino,
                  color = actino
                  )) +
  geom_bar(stat="count",position = "fill",  size = 1.5) +
  labs(x="Barcode",y="Reads") +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  theme_prism() +
  viridis::scale_fill_viridis(discrete = T,na.value ="grey") +
  theme(axis.text.x = element_text(size=6),
        legend.position = "bottom") +
  scale_color_manual(values = c("Actinobacteria" = "brown", "non-Actinobacteria" = "white")) +
  guides(color = guide_legend(nrow=2))





p_readlen = ggplot(barcode_summary %>% filter(barcode_arrangement!="unclassified" &
                                                seqlen != 0),
       aes(x=str_extract(barcode_arrangement,'\\d+'),
           y=seqlen/1000
       )) +
  geom_boxplot(color="black",fill="brown") +
  labs(x="Barcode",y="Read length (Kb)") +
  scale_y_continuous(expand=c(0,0), limits = c(0,10), breaks = 1:10) +
  theme_prism() +
  viridis::scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text(size=6)) 









# plotly::ggplotly(p_flongle_qc)


p_throughput = ggplot(throughput,aes(x=experiment_time_minutes/60,y=estimated_bases)) +
  geom_line(size=1.5) + theme_prism() +
  scale_y_continuous(label=scales::comma) +
  labs(y="Estimated bases", x= "Time (hrs)")



#save figures ####
ggsave(p_taxonomy,filename = sprintf("figures/flongleQC/%s/taxonomy.svg",folder), width = 11, height = 7)
ggsave(p_readlen,filename = sprintf("figures/flongleQC/%s/readlength.svg",folder), width = 11, height = 7)
ggsave(p_nreads,filename = sprintf("figures/flongleQC/%s/nreads.svg",folder), width = 11, height = 7)
ggsave(p_nbases,filename = sprintf("figures/flongleQC/%s/nbases.svg",folder), width = 11, height = 7)
ggsave(p_throughput,filename = sprintf("figures/flongleQC/%s/throughput.svg",folder), width = 7, height = 5)

write_csv(slectedpostqc,sprintf("figures/flongleQC/%s/selectedpostqc.csv",folder) ) 