"ctg1_183_lanth2_miss2_acetyl"

ms_collected_hits_idx=1
spectra


calc_fragments_adj = function(spectra,ms_collected_hits_idx,ann_prec=T){

pos_allowedlanth = which(cumsum(str_detect(unlist(str_extract_all(ms_collected_hits$core_seq[ms_collected_hits_idx],pattern =  ".")),pattern="S|T"))>ms_collected_hits$n_lanth[ms_collected_hits_idx])[1] 
core_seqlength = str_count(ms_collected_hits$core_seq[ms_collected_hits_idx],"[A-Z]")
pos_allowedlanth_cterm = core_seqlength-pos_allowedlanth


tailoring_mod=str_extract(ms_collected_hits$locus_tag[ms_collected_hits_idx],pattern = "miss\\d+_\\w+") %>% 
  str_remove("miss\\d+_")
cterm_mod=ifelse(tailoring_mod=="avicys"&!is.na(tailoring_mod),41.0265491,0)
nterm_mod=ifelse(tailoring_mod=="acetyl"&!is.na(tailoring_mod),42.0106,0)

additional_middle_tailoring=switch(tailoring_mod,"trpcl"=c("W"=17.96611),"hydroxyasp"=c("D"=15.99491462),c())

fragments = calculateFragments(
  sequence = ms_collected_hits$core_seq[ms_collected_hits_idx],
  modifications = c(
    S = -18.0105,
    'T' = -18.0105,
    C = 0,
    additional_middle_tailoring,
    Nterm = nterm_mod,
    Cterm = cterm_mod
  ),
  type = c("a", "b", "c", "x", "y", "z"),
  z = 1:4,
)

# adjust for missing dehydration sites
dehydrate_reverse = 18.0105
fragments=mutate(fragments,mz=case_when(
                              pos>=pos_allowedlanth&str_detect(type,"a|b|c") ~ mz + dehydrate_reverse*(str_count(seq,"S|T")- ms_collected_hits$n_lanth[ms_collected_hits_idx])/z,
                              pos>pos_allowedlanth_cterm&str_detect(type,"x|y|z") ~ mz + dehydrate_reverse*(ms_collected_hits$n_miss[ms_collected_hits_idx])/z,
                              pos<=pos_allowedlanth_cterm&str_detect(type,"x|y|z") ~ mz + (dehydrate_reverse*(str_count(seq,"S|T"))/z),
                              #pos<pos_allowedlanth&str_detect(type,"a|b|c") ~ mz + (dehydrate_reverse*(ms_collected_hits$n_lanth[ms_collected_hits_idx])/z),
                              TRUE~mz))
fragments = arrange(fragments,mz)

#find matches theoretical vs obs
fragments=fragments[MALDIquant::match.closest(mz(spectra),fragments$mz,tolerance = 0.1,nomatch = 0),] %>% 
  unique() %>% dplyr::rename(mz_th=mz)
fragments$mz = mz(spectra)[MALDIquant::match.closest(fragments$mz,mz(spectra),tolerance = 0.1,nomatch = 0)] 
fragments$intensity = intensity(spectra)[MALDIquant::match.closest(fragments$mz,mz(spectra),tolerance = 0.1,nomatch = 0)] 

if(ann_prec==T){
fragments =bind_rows(fragments,
tibble(mz_th=select(ms_collected_hits[ms_collected_hits_idx,],starts_with("mz_") & contains(as.character(precursorCharge(spectra))) )[[1]],
       mz=mz(spectra)[ifelse(MALDIquant::match.closest(mz(spectra),precursorMz(spectra),tolerance=0.1,nomatch = 0),T,F)],
       intensity=intensity(spectra)[ifelse(MALDIquant::match.closest(mz(spectra),precursorMz(spectra),tolerance=0.1,nomatch = 0),T,F)],
       z=precursorCharge(spectra),
       ion=ifelse(precursorCharge(spectra)==1,"M+H",ifelse(precursorCharge(spectra)==2,"M+2H","M+3H")),
       type="M") %>% head(1)
)
}

#calc and filter by ppm
fragments=mutate(fragments,error=((mz_th - mz)/mz_th)*10^6) %>% filter(abs(error)<50) 


return(fragments)}
