  library(tidyverse)
  
  gff = Sys.glob(file.path(getwd(), "data/antismash/*/*.gff.1"))
  gff = lapply(gff,function(x) rtracklayer::readGFF(x)) %>% setNames(.,gff)
  gff = bind_rows(gff,.id="gff_path")
  gff = mutate(gff,strain=str_extract(gff_path,pattern = "\\w+-\\w+"))
  

#
ripps = c("LAP","lanthipeptide","lassopeptide","linaridin","thiopeptide",
          "lipolanthine","microviridin","proteusin","sactipeptide",
          "glycocin","cyanobactin")


antismash6_ripps = c("cyclic-lactone-autoinducer","epipeptide","ranthipeptide","RRE-containing","spliceotide", # new additions
                     "RiPP-like","lanthipeptide-class-i","lanthipeptide-class-ii","lanthipeptide-class-iii","lanthipeptide-class-iv","lanthipeptide-class-v","thioamitides") # updates


ripps = c(ripps,antismash6_ripps)

# NHLP precursor
exclude_for_detection = c("RiPP-like", "RRE-containing")
regions = filter(gff,product %in% ripps &
                   type...3 == "sequence_feature" 
                   ) %>% 
  filter(!product %in% exclude_for_detection &
          is.na(candidate_cluster_number))

gff = mutate(gff, width = abs(start-end))
# core_locations = str_extract_all(regions$core_location,pattern = "\\d+")

precursorRegions = lapply(1:nrow(regions),function(idx)
  gff %>% filter(
    is.na(core_sequence) &
    gff_path == regions$gff_path[idx] &
      seqid == regions$seqid[idx] &
      start > regions$start[idx] &
      end < regions$end[idx] & 
      type...3 == "CDS" &
      width < 600 
  )
)

# to add ripp-like
which(sapply(precursorRegions,function(x) nrow(filter(x,str_detect(sec_met_domain,"lasso|Asn_synthase")|str_detect(gene_functions,"asparagine synthase|lasso|Asn_synthase"))) )>0)


precursors = select(bind_rows(precursorRegions),strain,locus_tag,translation,gene_kind,gene_functions) %>%
  mutate(name = paste(strain,locus_tag,sep = "@"),
         sequence = translation) %>% select(-c(strain,locus_tag,translation))

precursors = precursors %>% filter(is.na(gene_kind)|str_detect(gene_functions,"predicted lant")) 

# lanthipeptides 

lanth_mass=-18.0105

lactyl=2.9992
obu=0.98361
pyr=0.98361 
trpcl=17.96611307 #wrong
acetyl = 42.0106
hydroxyasp=15.99491462
avicys=-46.00578 # incorrectly used 41.0265491 previously

abuorala=2.015650064 # ser/thr after dha/dhb


precursors = precursors %>% mutate(cterminus_seq = str_sub(sequence,start=-20),
                      n_ser_thr = str_count(cterminus_seq,pattern = "S|T")) %>% 
filter(n_ser_thr>2 & str_detect(cterminus_seq,pattern = "C")) %>% distinct()


cores = sapply(5:35,function(i) str_sub(precursors$sequence,start=-i))



precursors = lapply(1:nrow(precursors),function(i){
  
  #for each core  
  lapply(cores[i,],function(core) {
    
    #basic info for each core
    n_lanth_pot=str_count(core,pattern = "S|T")
    if(n_lanth_pot<1) return(tibble(
      n_lanth_pot=0,
      n_lanth=0,
      n_miss=0,
      core_seq=core,
      mass_unmodified_core = Peptides::mw(core_seq,monoisotopic = T),
      mass_modified = mass_unmodified_core + (n_lanth*lanth_mass) + (n_miss*0), 
      name=sprintf("%s_lanth%d_miss%d",
                   precursors$name[i],
                   0,
                   0)
    ))
    
    lanth_perm=gtools::permutations(v=c("lanth","miss"),r=n_lanth_pot,n=2,repeats.allowed=T)
    
    tibble(n_lanth=sapply(1:nrow(lanth_perm), function(p) sum(str_count(lanth_perm[p,],"lanth"))),
           n_miss=sapply(1:nrow(lanth_perm), function(p) sum(str_count(lanth_perm[p,],"miss"))),
    ) %>% distinct() %>% 
      mutate(n_lanth_pot=n_lanth_pot,
             core_seq=core,
             mass_unmodified_core = Peptides::mw(core_seq,monoisotopic = T),
             mass_modified = mass_unmodified_core + (n_lanth*lanth_mass) + (n_miss*0), 
             name=sprintf("%s_lanth%d_miss%d",
                          precursors$name[i],
                          n_lanth,
                          n_miss)
      )
    
  })
  
  
}) 

precursors = bind_rows(precursors)

precursors = bind_rows(precursors,
precursors %>% mutate(mass_modified=mass_modified+trpcl,name=paste0(name,"_trpcl")) %>% filter(str_detect(core_seq,pattern = "W")),
precursors %>% mutate(mass_modified=mass_modified+acetyl,name=paste0(name,"_acetyl")), 
precursors %>% mutate(mass_modified=mass_modified+hydroxyasp,name=paste0(name,"_hydroxyasp")) %>% filter(str_detect(core_seq,pattern = "D|B")),
precursors %>% mutate(mass_modified=mass_modified+avicys,name=paste0(name,"_avicys")) %>% filter(str_detect(core_seq,pattern = "C$")),
precursors %>% mutate(mass_modified=mass_modified+lactyl,name=paste0(name,"_lactyl")) %>% filter(str_detect(core_seq,pattern = "^S|^T")),
precursors %>% mutate(mass_modified=mass_modified+obu,name=paste0(name,"_obu")) %>% filter(str_detect(core_seq,pattern = "^T")),
precursors %>% mutate(mass_modified=mass_modified+pyr,name=paste0(name,"_pyr")) %>% filter(str_detect(core_seq,pattern = "^S")),
precursors %>% mutate(n_cys = str_count(core_seq,"C"),
       disulfide = if_else(n_cys<=2,0,if_else(n_cys<=3,1,2)),
       mass_modified=mass_modified-(2*disulfide*1.00784),
       name = paste0(name,"_ds",disulfide)
) %>% filter(disulfide>0)
)

precursors = bind_rows(precursors,
mutate(filter(precursors,disulfide==2),disulfide=1,mass_modified=mass_modified+(2*disulfide*1.00784))
)

precursors = filter(precursors,n_lanth!=0 & str_detect(core_seq,pattern = "C"))





write_tsv(precursors,"data/misc/lanthi_mono.tsv")


