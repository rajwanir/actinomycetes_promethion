
library(tidyverse)
library(MSnbase)
library(cowplot)
library(grid)
library(gridExtra)

check_ms2_colllection = function(idx) {
  return(ifelse(
    length(msLevel(filterPrecursorMz(mzXML[[idx]],mz=hits$mz_1[idx],ppm=20))) >0  | 
      length(msLevel(filterPrecursorMz(mzXML[[idx]],mz=hits$mz_2[idx],ppm=20))) >0  |
      length(msLevel(filterPrecursorMz(mzXML[[idx]],mz=hits$mz_3[idx],ppm=20))) >0   
      ,T,F))
}




calc_fragments = function(ms_collected_hits_idx, sepctra_ce_idx,spectra,cterm_mod,nterm_mod) {
  return(
    calculateFragments(
      sequence = ms_collected_hits$core_seq[ms_collected_hits_idx],
      object = spectra[[sepctra_ce_idx]],
      modifications = c(
        S = -18.0105,
        T = -18.0105,
        Nterm = nterm_mod,
        Cterm = cterm_mod
      ),
      type = c("a", "b", "c", "x", "y", "z"),
      z = 1:4,
      tolerance = 1e-04,
      relative = T
    )
  )
}


annotate_corseq=function(coreseq,max_lanth){
  coreseq_vec = unlist(strsplit(coreseq,split=""))
  max_aaidx = max(str_which(coreseq_vec,"S|T")[1:max_lanth])
  coreseq_vec=sapply(1:length(coreseq_vec),function(aa_idx) 
   # switch(aa,"T"="<span style='color:#0000FF;'>T</span>","S"="<span style='color:#00FF00;'>S</span>",aa))
    ifelse(aa_idx<=max_aaidx&coreseq_vec[aa_idx]%in%c("S","T"),
           ifelse(coreseq_vec[aa_idx]=="T",
                  "<span style='color:#FFA500;'><b>T</b></span>",
                  "<span style='color:#00FF00;'><b>S</b></span>"),
           coreseq_vec[aa_idx]))
  return(paste(coreseq_vec,collapse = ""))
}
annotate_corseq=Vectorize(annotate_corseq)






hits = read_tsv("data/ripps/hits.tsv")
eic_pass= read_csv("data/ripps/eic_pass.csv")
msms_pass= read_csv("data/ripps/msms_pass.csv")

hits=left_join(hits,eic_pass) %>% filter(eic_pass==T)

# hits = mutate(hits,medium=str_extract(file,"ISP1|ISP2"),
#               mzXML=case_when(medium=="ISP1" ~ sprintf("data/QTOF/scans/%s.mzXML",str_replace(file,".d","") %>% str_replace("MS1","MS2") ),
#                               medium=="ISP2" ~ sprintf("data/QTOF/scans/%s.mzXML",str_replace(file,".d","") %>% str_replace("MS1","MS2")   ))
# )  


hits$mzXML=lapply(hits$strainprismdb,function(strain) {
  c(Sys.glob(sprintf("data/QTOF/scans/*%s*_MS2.mzXML",strain)),
    Sys.glob(sprintf("data/QTOF/scans/lanthi/*%s*_MS2.mzXML",strain))) }
  )

#read mzXML files
mzXML = lapply(hits$mzXML,function(x) readMSData(x,mode = "onDisk",msLevel. = 2))
#denoise the spectrum
mzXML = lapply(mzXML,function(x)removePeaks(x,t=80))



#check if MS2 spectra is collected
hits$ms_ms_collect = sapply(1:nrow(hits),function(idx) check_ms2_colllection(idx))


# 
# #redirect mzXML path to lanthi folder if MS2 not collected
# hits=mutate(hits,
#        mzXML=case_when(ms_ms_collect==F ~ str_replace(mzXML,pattern = "/scans/",replacement = "/scans/lanthi/") %>% str_replace("_[A-Z][0-9]+",""),
#                        ms_ms_collect==T ~ mzXML))
# 
# #repeat mzXML reading and ms/ms check
# mzXML = lapply(hits$mzXML,function(x) readMSData(x,mode = "onDisk",msLevel. = 2))
# hits$ms_ms_collect = sapply(1:nrow(hits),function(idx) check_ms2_colllection(idx))




ms_collected_hits = hits %>% filter(ms_ms_collect==T) 
ms_collected_mzXML = mzXML[hits$ms_ms_collect==T]
rm(hits)
rm(mzXML)

ms_collected_hits=mutate(ms_collected_hits,corseq_ann=annotate_corseq(core_seq,n_lanth)) 


p_spectralist = lapply(1:nrow(ms_collected_hits),function(idx){
return_spectraplot=T
print(idx)
#retreiving all spectra for given peptide (all z and all ce)
spectra = c(spectra(filterPrecursorMz(ms_collected_mzXML[[idx]],mz=ms_collected_hits$mz_1[idx],ppm=20)),
spectra(filterPrecursorMz(ms_collected_mzXML[[idx]],mz=ms_collected_hits$mz_2[idx],ppm=20)),
spectra(filterPrecursorMz(ms_collected_mzXML[[idx]],mz=ms_collected_hits$mz_3[idx],ppm=20)))

#denoise the spectrum before selecting best
#spectra=lapply(spectra,function(x) removePeaks(x,t=max(intensity(x)) * 0.05))
#spectra=lapply(spectra,function(x) removePeaks(x,t=200))
#spectra=lapply(spectra,function(x) clean(x,all = TRUE))

#check which spectra ce explains the peptide fragment best
best_spectra_idx = sapply(1:length(spectra),function(spectra_idx) {
   #sum(calc_fragments_adj(ms_collected_hits_idx=idx,spectra_idx,spectra=spectra,cterm_mod,nterm_mod)$intensity)/max(intensity(spectra[[spectra_idx]]))
  sum(calc_fragments_adj(spectra[[spectra_idx]],idx)$intensity,ann_prec=F)/max(intensity(spectra[[spectra_idx]]))
  }) %>% which.max()
spectra = spectra[[best_spectra_idx]]

#denoise the spectrum after selecting best
spectra=removePeaks(spectra,t=max(intensity(spectra)) * 0.05)
spectra=removePeaks(spectra,t=200)
spectra=clean(spectra,all = TRUE)


spectra_tbl=calc_fragments_adj(spectra,idx)
#spectra_tbl = calc_fragments(ms_collected_hits_idx=idx,best_spectra_idx,spectra=spectra,cterm_mod,nterm_mod) %>% mutate(error=(error/mz)*10^6)


write_csv(as.data.frame(spectra),sprintf("data/ripps/msms/msnbase/%s@%s_ce%d_z%d.csv",
                                         ms_collected_hits$strainprismdb[idx],
                                         ms_collected_hits$locus_tag[idx],
                                         collisionEnergy(spectra),
                                         precursorCharge(spectra)))

spectra = left_join(as.data.frame(spectra),
                    spectra_tbl,by="mz") %>% select(-intensity) %>%
  mutate(intensity=i,error=(error/mz)*10^6,
         type2 = str_extract(type,"a|b|c|x|y|z|M")) #%>%
  #filter(mz < precursorMz(spectra)+50)








p_spectra = ggplot(spectra,aes(mz,intensity)) +
  geom_bar(stat="identity",width=0.7,fill="black") +
  #  scale_y_continuous(limits = c(0,max(msms$intensity,na.rm = T)+200))+
  coord_cartesian(clip = "off") +
  ggrepel::geom_text_repel(aes(label=ion,color=type2),
                           min.segment.length=0,show.legend=F,
                           ylim = c(max(spectra$intensity,na.rm =T),NA),
                           segment.linetype="dashed",segment.size=0.6,fontface="plain",
                           direction="y") +
  theme_test() +
  theme(axis.line.x = element_line(color="black",size=1),
        axis.text.x = element_text(color="black"),
        axis.title.x = element_text(face="italic"),
        axis.line.y=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  labs(x="m/z") +
  scale_color_manual(values = c("x"="red","y"="red","z"="red",
                                "a"="lightblue","b"="lightblue","c"="lightblue",
                                "M"="black")) +
  labs(title = sprintf("%s@%s",
                       ms_collected_hits$strainprismdb[idx],
                       ms_collected_hits$locus_tag[idx]
                       ),
       subtitle = sprintf("%s",ms_collected_hits$corseq_ann[idx])) +
  theme(plot.title = element_text(hjust =0.5,size = 13),
        plot.subtitle = ggtext::element_markdown(lineheight = 1.1,hjust =0.5,size=14))


if(isTRUE(return_spectraplot)){return(p_spectra)}


spectra_tbl_lastrow=tibble(seq=sprintf("%d additional ion matches",
                                       ifelse(10<nrow(spectra_tbl),nrow(spectra_tbl)-10,0)  ))
spectra_tbl = arrange(spectra_tbl,-intensity) %>% head(10)
spectra_tbl = bind_rows(spectra_tbl,spectra_tbl_lastrow)
spectra_tbl = mutate_if(spectra_tbl,is.numeric, round,5)
spectra_tbl = replace(spectra_tbl,is.na(spectra_tbl),"")
spectra_tbl = dplyr::rename(spectra_tbl,mz_obs=mz,mz_calc=mz_th)
spectra_tbl = spectra_tbl[,c("ion","z","seq","intensity","mz_obs","mz_calc","error")]

p_spectra_tbl = gridExtra::tableGrob(spectra_tbl,theme=gridExtra::ttheme_minimal(),
                                     rows=NULL)

p_spec_n_tbl = cowplot::plot_grid(p_spectra,p_spectra_tbl,
                                  ncol=1,rel_heights = c(0.4,1))

ggsave(p_spec_n_tbl,filename = sprintf("figures/ripps/msms/%s@%s.svg",
                                       ms_collected_hits$strainprismdb[idx],
                                       ms_collected_hits$locus_tag[idx]),
       width = 10,height = 7)


})



p_spectralist_c = plot_grid(plotlist = p_spectralist,ncol=2,labels = LETTERS[1:length(p_spectralist)])
ggsave(p_spectralist_c,filename = "figures/ripps/lanthi_msmspass.svg",width = 25,height=13)