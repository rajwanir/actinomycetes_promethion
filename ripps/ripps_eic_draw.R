library(tidyverse)
library(MSnbase)
library(cowplot)

get_ppm_window = function(mass_target){
  mass_window = c( mass_target - (mass_target * 20 / 10^6),
                   mass_target + (mass_target * 20 / 10^6))
  return(mass_window)}



get_eic=function(idx,z=1){
  mass=ifelse(z==1,
         hits$massprismdb[idx],
          ifelse(z==2,hits$massprismdb[idx]/2,
                 hits$massprismdb[idx]/3))
  chr= chromatogram(mzXML[[idx]],mz=get_ppm_window(mass)+1.00784,missing=1,rt=c(60,600))
  chr=as.data.frame(chr[1,1])
  p=ggplot(chr,aes(rtime/60,intensity)) + geom_line() + theme_classic() +
       scale_x_continuous(breaks = 1:10) +
       labs(x="RT(min)",y="intensity",
       title=paste(hits$strainprismdb[idx],
                   hits$locus_tag[idx],
                   sep = "\n"),
       subtitle = sprintf("mass : %.1f, z : %d",
                          hits$massprismdb[idx],z)
       ) +
       theme(axis.text = element_text(color="black"),
             plot.title = element_text(hjust=0.5,size=11),plot.subtitle = element_text(hjust = 0.5))
  return(p)
}


draw_plot_title = function(z=1) {
  ggdraw() + 
    draw_label(sprintf("z%d",z),fontface = 'bold',x = 0,hjust = 0) +
    theme(plot.margin = margin(0,0,0,7))
}


hits = read_tsv("data/ripps/hits.tsv")

hits = mutate(hits,medium=str_extract(file,"ISP1|ISP2"),
                mzXML=case_when(medium=="ISP1" ~ sprintf("data/QTOF/scans/%s.mzXML",str_replace(file,".d","")),
                                medium=="ISP2" ~ sprintf("data/QTOF/scans/%s.mzXML",str_replace(file,".d","") %>% str_replace("MS1","MS2")   ))
)              

mzXML = lapply(hits$mzXML,function(x) readMSData(x,mode = "onDisk",msLevel. = 1))


#get eic plots
z1=lapply(1:nrow(hits),get_eic)
z2=lapply(1:nrow(hits),function(i) get_eic(idx=i,z=2))
z3=lapply(1:nrow(hits),function(i) get_eic(idx=i,z=3))



for(i in seq(1,nrow(hits),by=5)) {

start_idx=i
end_idx=ifelse((i+4)>nrow(hits),nrow(hits),i+4)
  
#combine eic plots by charge
p_z1=plot_grid(plotlist = z1[start_idx:end_idx],ncol=1)
p_z2=plot_grid(plotlist = z2[start_idx:end_idx],ncol=1)
p_z3=plot_grid(plotlist = z3[start_idx:end_idx],ncol=1)

#add titles
#plot_grid(draw_plot_title(z=1),plotlist = z1,ncol=1,rel_heights = c(0.1,0.9))

p_z1_z2_z3 = plot_grid(
  p_z1,
  p_z2,
  p_z3,
labels=c("A","B","C"),
ncol = 3)


ggsave(p_z1_z2_z3,filename = sprintf("figures/ripps/lanthi_eics_%d_%d.svg",start_idx,end_idx),
       width = 20,height = 11)

}