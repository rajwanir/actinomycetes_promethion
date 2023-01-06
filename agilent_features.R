library(tidyverse)

mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}


#formula to get ppm window around a mass
get_ppm_window = function(mass_target){
  mass_window = c( mass_target - (mass_target * 20 / 10^6),
                   mass_target + (mass_target * 20 / 10^6))
  return(mass_window)}

match_mass_strain = function(idx){
  mass=agilent_features$mass[idx]
  prism_strains_samemass = prismdb$strain[prismdb$mass > get_ppm_window(mass)[1] & prismdb$mass < get_ppm_window(mass)[2]]
  return(agilent_features$strain[idx] %in% prism_strains_samemass)
}


get_prism_match = function(idx){
  mass=hits$mass[idx]
  #agilent_strain=agilent_features$strain[idx]
  miniprismdb = filter(prismdb,strain==hits$strain[idx])
  miniprismdb = miniprismdb[miniprismdb$mass > get_ppm_window(mass)[1] & miniprismdb$mass < get_ppm_window(mass)[2],]
  return(miniprismdb)
}



filter_media = function(idx){
  mass=agilent_features$mass[idx]
  rt=agilent_features$rt[idx]
  media_mass = media[media$mass > get_ppm_window(mass)[1] & media$mass < get_ppm_window(mass)[2],]
  media_check = ifelse(nrow(media_mass)>0, # if  mass matches in media,
         !any(abs(media_mass$rt-rt)<0.5), # check retention time return F if match
         T) # returns T if not detected in media
  return(media_check)
}

count_additional_strains = function(idx){
  strain=agilent_features$strain[idx]
  mass=agilent_features$mass[idx]
  rt=agilent_features$rt[idx]
  additional_strains = agilent_features[agilent_features$strain!=strain & (agilent_features$rt-rt)<0.5 & 
                                               agilent_features$mass > get_ppm_window(mass)[1] &
                                               agilent_features$mass < get_ppm_window(mass)[2],"strain"]
  return(nrow(unique(additional_strains)))
}

agilent_features = lapply(Sys.glob("data/QTOF/db_screens/agilent_features/*.csv"),function(x)read_csv(x,comment="#")) %>%
  bind_rows() %>% 
  janitor::clean_names() # %>% filter(is.na(ms_ms_count))
agilent_features$strain = str_extract(agilent_features$file,pattern = "G\\w+-\\d+")


agilent_features = bind_rows(parallel::mclapply(Sys.glob("../actinomycetes_expression/data/msdiscoverLibrary/features/griseus/*.csv"),function(x) read_csv(x,comment = "#",skip = 2) %>% janitor::clean_names(),mc.cores=50)) %>% mutate(strain="griseus")
agilent_features = agilent_features %>% filter(max_z>1|is.na(max_z))


#filter anything in media
media = read_csv("data/QTOF/db_screens/agilent_features/media/allmedia_peptides.csv") %>% 
  janitor::clean_names() # %>% filter(is.na(ms_ms_count))
media = filter(media,str_detect(file,"R2A"))
media=bind_rows(media,
          read_csv("data/QTOF/db_screens/agilent_features/media/thiostrepton_R2A.csv") %>%janitor::clean_names())
#background in moldiscovery
media = bind_rows(parallel::mclapply(Sys.glob("../actinomycetes_expression/data/msdiscoverLibrary/features/GA13-001/*.csv"),function(x) read_csv(x,comment = "#",skip=2) %>% janitor::clean_names(),mc.cores=50)) %>% mutate(strain="GA13-001")
background2 =parallel::mclapply(Sys.glob("data/QTOF/db_screens/agilent_features/*.csv"),function(x)read_csv(x,comment="#"),mc.cores=50) %>%
  bind_rows() %>% 
  janitor::clean_names() %>% mutate(strain=str_extract(file,pattern = "G\\w+-\\d+")) %>% 
  filter(str_detect(file,"R2A") & strain!="GA3-008")
media=bind_rows(media,background2)
media = media %>% filter(max_z>1|is.na(max_z))

agilent_features = agilent_features[mcsapply(1:nrow(agilent_features),filter_media,mc.cores=50),]


npatlas = read_tsv("../databases/NPAtlas_download_10132021.tsv")


prismdb = read_tsv("data/prism/prism_structures.tsv")
prismdb$strain = str_extract(prismdb$prism_path,pattern = "\\w+-\\d+")


ripps = read_tsv("data/misc/antismash_predicted_core_ripps.tsv") %>%
  rename(strain=strain,mass=mass_modified_core)

lassopeptides = read_tsv("data/misc/lassopeptides.tsv") %>%
  separate(name,into=c("strain","locus_tag"),sep="@") %>% 
  dplyr::rename(mass=mass_modified_core)

LAP = read_tsv("data/misc/LAP_mono.tsv") %>% 
  separate(name,into=c("strain","locus_tag"),sep="@") %>% 
  dplyr::rename(mass=mass_modified)

lanthi = read_tsv("data/misc/lanthi_mono.tsv") %>% 
  separate(name,into=c("strain","locus_tag"),sep="@") %>% 
  dplyr::rename(mass=mass_modified)

sacti = read_tsv("data/misc/sactipeptide_mono.tsv") 

thiopeptide = read_csv("data/ripps/thiopeptide_GA3-008.csv") %>% 
  separate(locus_tag,into=c("strain","locus_tag"),sep="@") 


GA3008db = read_csv("data/ripps/GA3008db.csv")
GA13001db = read_csv("data/ripps/GA13-001/GA13-001db.csv")
GA6010db = read_csv("data/ripps/GA6-010/GA6-010db.csv")
griseus = read_csv("data/ripps/griseus/griseusdb1.csv")

prismdb = read_csv("data/ripps/rippsdb_allclass_allstrains.csv")

# check if any agilent feature mass within 20ppm window of prism mass and matches by strain too
# agilent_features[sapply(agilent_features$mass,function(mass) {
#   any(prismdb$strain[prismdb$mass > get_ppm_window(mass)[1] & prismdb$mass < get_ppm_window(mass)[2]] %in%
#     unique(agilent_features$strain[agilent_features$mass==mass]))
# }),]


# count number of strains with same mass (20ppm window) 
# agilent_features$n_strains = sapply(agilent_features$mass,function(mass) {
#   length(unique(agilent_features$strain[agilent_features$mass > get_ppm_window(mass)[1] & agilent_features$mass < get_ppm_window(mass)[2]]))
# })

agilent_features$n_strains = sapply(1:nrow(agilent_features),count_additional_strains)








#check if mass matches to prismdb
hits = agilent_features[mcsapply(1:nrow(agilent_features),match_mass_strain,mc.cores=50),]

# filter by intensity, select columns and combine
hits = hits %>% 
  mutate(match_idx=as.character(1:n())  ) %>%
  select(mass,rt,height,ions,ms_ms_count,file,strain,match_idx)

  hits = left_join(hits,
  bind_rows(parallel::mclapply(1:nrow(hits),get_prism_match,mc.cores = 50),
            .id="match_idx"),
  by="match_idx",suffix=c("agilent","prismdb")
  )
  
  #write m/z
  hits=mutate(hits,mz_1=massprismdb+1.00784,mz_2=((massprismdb/2)+1.00784),mz_3=((massprismdb/3)+1.00784),
              ppm=abs(((massprismdb-massagilent)/massprismdb)*10^6)  )
  
  hits = hits %>%  select_if(function(x) !(all(is.na(x)) | all(x=="")))
  
hits = distinct(hits,strainprismdb,massprismdb,locus_tag,core_seq,.keep_all = T)

hits=filter(hits,n_strains<20)


write_tsv(hits,"data/ripps/hits.tsv")


##

hits = agilent_features[sapply(1:nrow(agilent_features),match_mass_strain),]

# filter by intensity, select columns and combine
hits = hits %>% filter(height>10000) %>% arrange(-volume) %>% 
  mutate(match_idx=as.character(1:n())  ) %>%
  select(mass,rt,height,volume,precursor,ions,vol_percent,ms_ms_count,file,strain,match_idx)

hits = left_join(hits,
                 bind_rows(lapply(1:nrow(hits),get_prism_match),.id="match_idx"),
                 by="match_idx",suffix=c("agilent","prismdb")
)

#write m/z
hits=mutate(hits,mz_1=massprismdb+1.00784,mz_2=((massprismdb/2)+1.00784),mz_3=((massprismdb/3)+1.00784),
            ppm=abs(((massprismdb-massagilent)/massprismdb)*10^6)  )

###
