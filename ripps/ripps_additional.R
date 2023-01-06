jsonify::to_json(precursors) %>% write(.,file = "deepripp/precursors_for_deepripp.json")

deeprip_classifications = jsonify::from_json("deepripp/class_predictions.out_lap.json",simplify = T) 
deeprip_classifications$class_predictions$name = deeprip_classifications$name
deeprip_classifications = deeprip_classifications$class_predictions 
deeprip_classifications$nclassifications = sapply(deeprip_classifications$score,length)
deeprip_classifications = deeprip_classifications %>% filter(nclassifications ==1) %>%
  mutate_if(is.list,as.vector)
deeprip_classifications = deeprip_classifications %>% filter(class=="THIOPEPTIDE")


deeprip_cleavage = jsonify::from_json("deepripp/cleavage_predictions.out_lap.json",simplify = T) 
deeprip_cleavage = deeprip_cleavage$cleavage_prediction
deeprip_cleavage %>% filter(name %in% deeprip_classifications$name)
deeprip_cleavage = deeprip_cleavage %>% filter(name %in% deeprip_classifications$name) %>% unique()
deeprip_cleavage = deeprip_cleavage %>% merge(.,deeprip_classifications,by="name",suffix=c("_cleavage","_class")) %>% unique()


deeprip_cleavage = deeprip_cleavage %>% mutate(n_ser_thr_cys = str_count(sequence,"S|T|C"), ratio = n_ser_thr_cys/length(sequence))


deeprip_cleavage = mutate(deeprip_cleavage,
                          mass_unmodified_core = Peptides::mw(sequence,monoisotopic=T),
                          mass_modifified_core = mass_unmodified_core - (18.0105 * n_ser_thr_cys))





write_tsv(deeprip_cleavage,"deepripp/cleavage_predictions.out.lap.filtered.tsv")


# gets the antismash predicted cores and modified masses
antismash_predict_cores = gff %>% filter(!is.na(core_sequence)) %>% select(strain,locus_tag,peptide,core_sequence,leader_sequence,monoisotopic_mass,predicted_class)
antismash_predict_cores = mutate(antismash_predict_cores,gbk_path = sprintf("data/antismash/%s/%s_homopolished.gbk",strain,strain))


antismash_predict_cores = antismash_predict_cores %>% mutate(mass_unmodified_core = Peptides::mw(core_sequence,monoisotopic=T),
                                                             n_ser_thr = str_count(core_sequence,pattern = "S|T"),
                                                             n_cys = str_count(core_sequence,pattern = "C"),
                                                             mass_modified_core = case_when(peptide == "lassopeptide" ~ mass_unmodified_core - 18.0105,
                                                                                            peptide == "lanthipeptide" ~ mass_unmodified_core - (18.0105 * n_ser_thr) )
)

antismash_predict_cores = mutate(antismash_predict_cores,mass_diff = mass_unmodified_core-monoisotopic_mass)

write_tsv(antismash_predict_cores,file = "data/misc/antismash_predicted_core_ripps.tsv")


#lassopepetides #######


# based on leader
# lasso_leader_motif = "G[A-Z][AVFI][A-Z]{3}T[A-Z]" # PMID: 26079760
# lasso_leader_motif="[A-Z]+P[A-Z][ILVAM]([A-Z]){5,17}T[A-Z]" 
# 
# precursors = precursors %>% filter(str_detect(sequence,pattern = lasso_leader_motif))
# precursors$leaderseq =  str_split(precursors$sequence,pattern = lasso_leader_motif,n=2,simplify = T)[,1]
# precursors$core_sequence = str_split(precursors$sequence,pattern = lasso_leader_motif,n=2,simplify = T)[,2]


precursors = precursors %>% filter(str_detect(core_sequence,pattern="D|E"))
precursors = precursors %>% mutate(core_length = nchar(core_sequence)) %>% 
  filter(core_length > 6) %>% filter(core_length < 50) 



#based on core

get_lasso_rings = function(core_peptide,name){
  start_pos=str_locate_all(core_peptide,"G|I|L|S|C")[[1]][,1]
  end_pos=str_locate_all(core_peptide,"D|E")[[1]][,1]
  
  matches = expand.grid(start_pos,end_pos) %>% dplyr::rename(start_pos=Var1,end_pos=Var2) %>% mutate(ringsize=end_pos-start_pos) %>% 
    filter(between(ringsize,4,9))
  
  matches$ringseq = sapply(1:nrow(matches),function(idx) str_sub(core_peptide,
                                                                 matches$start_pos[idx],
                                                                 matches$end_pos[idx]))
  
  matches$tailseq = sapply(1:nrow(matches),function(idx) str_sub(core_peptide,
                                                                 matches$end_pos[idx]+1))
  matches=mutate(matches,core_peptide = paste0(ringseq,tailseq),
                 name = paste0(name,"_rs",start_pos,"re",end_pos))
  
  return(matches)
}




precursors = precursors %>% mutate(cterminus_seq = str_sub(sequence,start=-30)) %>% 
  filter(str_detect(cterminus_seq,"(G|I|L|S|C)[A-Z]{4,7}(D|E)[A-Z]{6,}")  & 
           (is.na(gene_kind)|str_detect(gene_functions,"predicted lassopept")))
  
  
precursors = lapply(1:nrow(precursors),function(i) get_lasso_rings(precursors$cterminus_seq[i],precursors$name[i])) %>% bind_rows()


cores = sapply(2:7,function(i) str_sub(precursors$core_peptide,end=-i))

precursors=tibble(name = as.vector(sapply(0:6,function(i) paste0(precursors$name,"_fl",i))),
       core_sequence = c(precursors$core_peptide,
                        as.vector(sapply(1:6,function(i)  cores[,i]  )) )
       )



precursors=filter(precursors,str_detect(core_sequence,"(G|I|L|S|C)[A-Z]{4,7}(D|E)[A-Z]{6,}"))

precursors=precursors %>% mutate(mass_unmodified_core = Peptides::mw(core_sequence,monoisotopic=T),
                                 mass_modified_core = mass_unmodified_core - 18.0105)

#add disulfide 
precursors=bind_rows(precursors,
mutate(precursors,
            n_cys = str_count(core_sequence,"C"),
                      difulfide = if_else(n_cys<=2,0,if_else(n_cys<=3,1,2)),
                      mass_modified_core=mass_modified_core-(2*difulfide*1.00784),
                      name = paste0(name,"_ds",difulfide)
) %>% filter(difulfide>0) )


write_tsv(precursors,"data/misc/lassopeptides_ripplike.tsv")


# LAPs

# ser2dha = - (2*2.015650) - 15.994915 # H2 and O
#azole = - (4*2.015650) - 15.994915 # H4 and O
  ser2dha = -18.0105
  azole = -20.02622
  amidine= -18.0105
  acetyl = 42.0106
  dimethyl= 28.0313
  pyridine=-17.026549100 -18.0105 #-NH3 - H20
  azole2azoline=2.01565007
  
  precursors = precursors %>% mutate(cterminus_seq = str_sub(sequence,start=-10),
                                     n_ser = str_count(cterminus_seq,pattern = "S"),
                                     n_cys = str_count(cterminus_seq,pattern = "C"),
                                     n_thr = str_count(cterminus_seq,pattern = "T"),
                                     n_azole = n_ser+n_cys+n_thr) %>% 
    filter(n_azole > 1) %>% distinct()
  
  
  cores = sapply(5:50,function(i) str_sub(precursors$sequence,start=-i))
  
  
  
  
  precursors = lapply(1:nrow(precursors),function(i){
    
    #for each core  
    lapply(cores[i,],function(core) {
      
      #basic info for each core
      n_ser=str_count(core,pattern = "S")
      n_cys = str_count(core,pattern = "C")
      n_thr = str_count(core,pattern = "T")
      n_azole = n_cys+n_thr
      
      if(n_ser<1) return(tibble(
        ser_azole=0,
        ser_dha=0,
        ser_miss=0,
        n_ser=0,
        n_cys=n_cys,
        n_thr=n_thr,
        n_azole=n_azole,
        core_seq=core,
        mass_unmodified_core = Peptides::mw(core_seq,monoisotopic = T),
        mass_modified = mass_unmodified_core + (n_azole*azole) + (ser_azole*azole) + (ser_dha*ser2dha) + (ser_miss*0), 
        name=sprintf("%s_azole%d_serazole%d_serdha%d_sermiss%d",
                     precursors$name[i],
                     n_azole,
                     ser_azole,
                     ser_dha,
                     ser_miss)
      ))
      
      ser_perm=gtools::permutations(v=c("dha","azole","miss"),r=n_ser,n=3,repeats.allowed=T)
      
      tibble(ser_azole=sapply(1:nrow(ser_perm), function(p) sum(str_count(ser_perm[p,],"azole"))),
             ser_dha=sapply(1:nrow(ser_perm), function(p) sum(str_count(ser_perm[p,],"dha"))),
             ser_miss=sapply(1:nrow(ser_perm), function(p) sum(str_count(ser_perm[p,],"miss"))),
      ) %>% distinct() %>% 
        mutate(n_ser=n_ser,
               n_cys=n_cys,
               n_thr=n_thr,
               n_azole=n_azole,
               core_seq=core,
               mass_unmodified_core = Peptides::mw(core_seq,monoisotopic = T),
               mass_modified = mass_unmodified_core + (n_azole*azole) + (ser_azole*azole) + (ser_dha*ser2dha) + (ser_miss*0), 
               name=sprintf("%s_azole%d_serazole%d_serdha%d_sermiss%d",
                            precursors$name[i],
                            n_azole,
                            ser_azole,
                            ser_dha,
                            ser_miss)
        )
      
    })
    
    
  }) 
  
  precursors = bind_rows(precursors)
  
  precursors=bind_rows(precursors,
                       mutate(precursors,mass_modified=mass_modified+acetyl,name=paste(name,"acetyl",sep = "_")),
                       mutate(precursors,mass_modified=mass_modified+amidine,name=paste(name,"amidine",sep = "_")),
                       mutate(precursors,mass_modified=mass_modified+dimethyl,name=paste(name,"dimethyl",sep = "_"))
  )


write_tsv(precursors,"data/misc/LAP_mono.tsv")




# sactipeptides ######

precursors = read_csv("data/ripps/sactipeptides_precursors.csv")


cleave_leader=function(precursor_seq) {
  cores= sapply(1:nchar(precursor_seq),function(s){
    sapply(1:nchar(precursor_seq),function(e)
      str_sub(precursor_seq,start=s,end=e))
  }) %>% as.vector()
  cores = cores[cores!=""] %>% unique()
  return(cores)}


core_tbl = lapply(precursors$sequence,cleave_leader)


core_tbl = tibble(strain = sapply(1:nrow(precursors), function(x) rep(precursors$strain[x],length(core_tbl[[x]])) ) %>% unlist(),
                  locus_tag = sapply(1:nrow(precursors), function(x) rep(precursors$locus_tag[x],length(core_tbl[[x]])) ) %>% unlist(),
                  core_seq = unlist(core_tbl),
                  core_length = nchar(core_seq),
                  unmodified_mass = Peptides::mw(core_seq,monoisotopic = T),
                  n_cys = str_count(core_seq,"C"),
                  mass=unmodified_mass-(n_cys*2*1.00784)) %>% 
  filter(n_cys>0 & core_length>4)

write_tsv(core_tbl,"data/misc/sactipeptide_mono.tsv")



# checking in QTOF through eic ##########
require(parallel)
require(MSnbase)
require(tidyverse)

antismash_ripps = read_tsv("data/misc/antismash_predicted_core_ripps.tsv")

#calc mz and charge
antismash_ripps = antismash_ripps %>% mutate(mz = case_when(mass_modified_core <= 2500 ~ (mass_modified_core/2)+1.00784 ,
                                                            mass_modified_core > 2500 ~ (mass_modified_core/3)+1.00784),
                                             charge = case_when(mass_modified_core <= 2500 ~ 2,
                                                                mass_modified_core > 2500 ~ 3))
#skip sactipeptide - no mass calc.
antismash_ripps = antismash_ripps %>% filter(peptide != "sactipeptide")
#antismash_ripps = antismash_ripps %>% separate(name,into=c("strain","locus_tag"),sep="@")

strains = unique(antismash_ripps$strain)
specPaths = sapply(strains,function(s) Sys.glob(sprintf("data/QTOF/scans/ISP2*%s*_MS2.mzXML",s)[1]))

specFiles = mclapply(specPaths,function(specFile) 
  readMSData(specFile,msLevel.=1,mode="onDisk"),mc.cores = 30)
names(specFiles) = strains


eics = mclapply(1:nrow(antismash_ripps),function(idx){
  strain = antismash_ripps$strain[idx]
  mz_target = antismash_ripps$mz[idx]
  mz = c( mz_target - (mz_target * 20 / 10^6),
          mz_target + (mz_target * 20 / 10^6))
  chromatogram(specFiles[[strain]], mz= mz,missing=1)
},mc.cores=60)


#p_eics = lapply(eics[33:48],function(e) plot(e,main=))
seq = seq(1,64,16)
lapply(seq, function(i) {
  svg(filename = sprintf("figures/QTOF/lassopeptides_ripps_eics_%s_%d.svg",i,i+15))
  plot.new()
  par(mfrow = c(4,4))
  for (idx in i:(i+15)) {
    plot(eics[[idx]],
         main = sprintf(
           "%s m/z %.2f",
           antismash_ripps$strain[idx],
           antismash_ripps$mz[idx]
         ))
  }
  dev.off()
})


write_tsv(antismash_ripps,"data/misc/antismash_predicted_core_ripps_mz_z.tsv")
