library(tidyverse)
library(MSnbase)
library(parallel)

get_ppm_window = function(mass_target){
  mass_window = c( mass_target - (mass_target * 20 / 10^6),
                   mass_target + (mass_target * 20 / 10^6))
  return(mass_window)}

peptide_family = read_csv("data/ripps/peptide_families/cirulassin.csv")


strains=as.character(peptide_family$strain)
media=c("ISP1","R2A","R2YE","SFM","ISP2")
group_colors = setNames(RColorBrewer::brewer.pal(length(media),"Dark2"),media)

medium_filepaths = sapply(strains,function(strain){
sapply(media,function(medium){
     medium_filepath = Sys.glob(sprintf("data/QTOF/scans/%s/*%s*",medium,strain))
},USE.NAMES = T)},USE.NAMES = T)

chrs = mclapply(medium_filepaths, function(medium_filepath) {
  
if (length(medium_filepath) != 0) {
    strain=str_extract(medium_filepath,"G?[A-Z]\\d+-\\d+")
    targetmz=as.numeric(peptide_family[peptide_family$strain%in%strain,'mz'][1,1])
  
    msdata = filterMz(readMSData(medium_filepath[[1]], msLevel. = 1, mode = "onDisk"),
                      mz = get_ppm_window(targetmz))
    chr = chromatogram(msdata, mz = get_ppm_window(targetmz), missing = 1)[1, 1] # calling chromatogram
    chr@rtime = rtime(chr) / 60 # transforming seconds to min
  } else{
    chr = medium_filepath
  }
  return(chr)
}, mc.cores = length(medium_filepaths))

names(chrs)=as.character(sapply(colnames(medium_filepaths),function(x) paste(x,rownames(medium_filepaths),sep="-")))

svg("figures/ripps/peptide_familes/cirulassin_eics.svg",width = 16)
par(mfcol = c(length(media),length(strains)), mar=c(4,4,3,3)) #bltr
lapply(strains,function(strain){
lapply(media, function(medium) {
  
  if (length(chrs[[paste(strain,medium,sep="-")]])!= 0) {
    chr = chrs[[paste(strain,medium,sep="-")]]
    max_y = if_else(max(intensity(chr))<2000,2000,max(intensity(chr)))
    plot(chr,
         main=sprintf("%s\n%s",strain,medium),col =group_colors[medium],
         ylab="Intensity",xlab="Retention time (min)", 
         xlim =c(2,10),ylim=c(0,max_y))
  }
  else{plot(0,type='n',axes=FALSE,ann=FALSE,main=sprintf("%s\n%s",strain,medium),col =group_colors[medium],xlab=NULL,ylab=NULL)}
  })
})
dev.off()


#to prepare file
# pep_family = read_csv("data/ripps/peptide_families/sapT.csv")
# n_loss=4
# water=18.0105
# hydrogen=1.00784
# pep_family=mutate(pep_family,
#                   mass=Peptides::mw(core_peptide,monoisotopic = T)-(n_loss*water),
#                   z=2,
#                   mz=(mass/z)+hydrogen)
# 
# write_csv(pep_family,"data/ripps/peptide_families/sapT.csv")