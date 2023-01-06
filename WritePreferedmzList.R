require(tidyverse)

hits = read_tsv("data/ripps/hits_lanthi.tsv") 
hits = hits %>% filter(ppm<10) %>% select(mass=massprismdb,strain=strainprismdb,rt) 
strains = unique(hits$strain)

#read template for preferedmzlist
preferedmzTemplate = read_csv("D:/NIH/QTOF/templates/prefered_mz_template.csv",skip = 1) %>% head(2) %>%
  mutate(On = if_else(On == T, "True","False"))

WritePreferedmzList = function(strain){
  x =strain
 mzList = hits %>% filter(strain %in% x) %>% distinct(mass,.keep_all=T) %>% mutate(On = "True",
                            "Prec. m/z" = ifelse(mass<1200,mass+1.00784,
                                                 ifelse(mass>2000,
                                                        (mass/3)+1.00784,
                                                        (mass/2)+1.00784
                                                        )
                                                 ),
                            "Delta m/z (ppm)" = 20,
                            Z = ifelse(mass<1200,1,ifelse(mass>2000,
                                                          3,2)),
                            "Prec. Type" = "Preferred",
                            "Ret. Time (min)" = rt,
                            "Delta Ret. Time (min)" = 0.3,
                            "Iso. Width" = "Medium (~4 m/z)",
                            "Collision Energy" = NA_real_
                            ) %>% 
  bind_rows(preferedmzTemplate,.) %>% select(-c(mass,strain,rt))
 
 cat("AutoPreferredExcludeMSMSTable\n",
     file=sprintf("data/QTOF/preferedmzList/%s_preferedmzList.csv",strain))

 write_csv(mzList,sprintf("data/QTOF/preferedmzList/%s_preferedmzList.csv",strain),
           na = "",
           append = T,
           col_names =T,
           eol="\r\n" )

}


lapply(strains,WritePreferedmzList)


# worklists setup

worklists = Sys.glob("D:/NIH/QTOF/promethion/worklists/*.csv")
worklists = bind_rows(lapply(worklists,read_csv)) %>% janitor::clean_names()
worklists = filter(worklists,str_detect(data_file,pattern = "MS1"))
worklists = mutate(worklists,medium=str_extract(data_file,pattern = "ISP1|ISP2|R2A|SFM|SC"),
                   strain = str_extract(sample_name,pattern = "G\\w+-\\d+"))

worklist_fltr = lapply(1:nrow(hits),function(idx) 
  filter(worklists, strain == hits$strainprismdb[idx] & medium == hits$medium[idx]) 
) %>% bind_rows() %>% distinct()



worklists = bind_rows(setNames(lapply(Sys.glob("D:/NIH/QTOF/promethion/*.csv"),read_csv),c("LAP","lanthi"))
                               ,.id="ripp_class")

worklists = mutate(worklists,
                   method=sprintf("D:\\MassHunter\\Methods\\rahim\\hits\\%s\\%s.m",ripp_class,strain),
                   data_file=sprintf("D:\\MassHunter\\Data\\rahim\\hits\\%s\\%s_%s_%s.d",ripp_class,medium,strain,sample_position))


worklists = worklists %>% filter(medium %in% c("ISP1","ISP2")) %>%
  mutate(sample_position=case_when(medium=="ISP1" ~ str_replace(sample_position,pattern = "P2|P1",replacement = "P1"),
                                   medium=="ISP2" ~ str_replace(sample_position,pattern = "P1|P2",replacement = "P2")))


write_csv(worklist_fltr,"D:/NIH/QTOF/promethion/worklists_lanthi.csv")