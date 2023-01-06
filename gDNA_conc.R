library(tidyverse)
library(janitor)


plate = 3
date = "06092021"
seqdate = "06152021"

barcodeKey = read_csv(sprintf("data/flongleQC/plate%d_flongle_%s/barcode_key.csv",plate,seqdate))


conc = list()

#read pg #####
conc$pg = read_csv(Sys.glob(sprintf("data/gDNA/picogreen_plate%d_%s.csv",plate,date))) %>% 
  pivot_longer(cols = where(is.numeric), names_to = "col", values_to = "pg_conc",
               names_transform = list(col = as.numeric)) %>% 
  rename(row = X1)  %>% 
  mutate(plate=plate,
         yield = pg_conc*12) 


  

#read nd ######
conc$nd = read_tsv(Sys.glob(sprintf("data/gDNA/nanodrop_plate%d_%s.tsv",plate,date) ))%>% 
  clean_names() 
conc$nd = mutate(conc$nd,
                        row = toupper(str_extract(sample_id,pattern = '[A-Z]')),
                            col = as.numeric(str_extract(sample_id,pattern = '\\d+')),
                        plate = plate)

#combine nd and pg #######
conc$comb = merge(conc$nd,conc$pg,
                         by=c("row","col","plate")) %>% 
  select(row,col,pg_conc,x260_280,x260_230,nucleic_acid,plate,yield)


#add barcode key
conc$comb = merge(conc$comb,barcodeKey, by=c("row","col"),
                  all.x = F,all.y=T)



# check ####
# conc$comb %>% filter(yield > 100) 

#approximate volume
conc$pg = conc$pg %>% mutate(DNA_vol_500 = case_when(500/pg_conc < 5 ~ "4",
                                                     500/pg_conc < 7.14 ~ "6",
                                                     100/pg_conc > 12 ~ "0",
                                                     TRUE ~ "12"))

conc$pg = conc$pg %>% mutate(nfw_vol_500 = 12-as.numeric(DNA_vol_500))


colors_dna_vol =  c("0" = "white", "8" = "green" , "6" = "brown", "12" = "blue")

#disualify
conc$pg = conc$pg %>% mutate(disqualify = ifelse(DNA_vol_500 == "0",T,F))


#amount cut
conc$pg = conc$pg %>% mutate(DNA_amnt_100 = ifelse(100/pg_conc < 13,">100ng","<100ng"))
colors_dna_amt = c(">100ng" = "green","<100ng" = "white")


#convert to 96-well format
p_yield = ggplot(data=conc$pg,aes(x=col,y=fct_relevel(row,levels = LETTERS[8:1]), 
                                  fill=disqualify)) +
  geom_tile(color = "black") +
#  geom_text(aes(label = disqualify), size =12)+
 # scale_fill_manual(values =colors_dna_vol) +
  scale_x_continuous(breaks = 1:12, expand = c(0,0),position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  ggtitle(sprintf("Plate %d",plate)) +
  ggprism::theme_prism() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Yield"))  +
  labs(y="row")


#od cut
conc$nd = conc$nd %>% mutate(OD260_280 = ifelse(x260_280 < 1.5,"impure","pure"))
conc$nd = conc$nd %>% mutate(OD260_230 = ifelse(x260_230 < 1.5,"impure","pure"))

colors_nd_purity = c("impure" = "white","pure" = "green")



p_od260280 = ggplot(data=conc$nd,aes(x=col,y=row, fill=OD260_280)) +
  geom_tile(color = "black") +
  geom_text(aes(label = x260_280))+
  scale_fill_manual(values =colors_nd_purity) +
  scale_x_continuous(breaks = 1:12, expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ggtitle(sprintf("Plate %d",plate)) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Absorbance 260/280")) 



p_od260230 = ggplot(data=conc$nd,aes(x=col,y=row, fill=OD260_230)) +
  geom_tile(color = "black") +
  geom_text(aes(label = x260_230))+
  scale_fill_manual(values =colors_nd_purity) +
  scale_x_continuous(breaks = 1:12, expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  ggtitle(sprintf("Plate %d",plate)) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = "Absorbance 260/230")) 

# DNAvol_ep = conc$pg$plate2 %>% group_by(col) %>% pivot_wider(names_from = col,id_cols = c(row,col), values_from=DNA_vol_500)
# write_csv(file = "data/DNA_vol_ep_plate2_06012021.csv", DNAvol_ep)


ggsave(p_yield,filename = sprintf("figures/dnaQC/yield_plate_%d_%s.svg",plate,date))
ggsave(p_od260280,filename = sprintf("figures/dnaQC/od260280_plate_%d_%s.svg",plate,date))
ggsave(p_od260230,filename = sprintf("figures/dnaQC/od260230_plate_%d_%s.svg",plate,date))
