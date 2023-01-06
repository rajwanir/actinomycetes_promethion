library(tidyverse)
library(MSnbase)

sigmatches = lapply(Sys.glob("data/QTOF/db_screens/prism/*/significant_unique_matches.tsv"),
                    function(x) read_tsv(x,col_types = cols(.default = "c")))

select(bind_rows(sigmatches),
       c(id=Name,smiles=SMILES)) %>% jsonlite::toJSON() %>% write(file="figures/prism_matches/rban_input.json")


sigmatches = bind_rows(sigmatches)

#use exact calc mass
prism_structures = read_tsv("data/prism/prism_structures.tsv") %>% 
  mutate(strain=str_extract(prism_path,pattern = "\\w+-\\d+"),
         name = paste(strain,name,sep="_")) %>% 
  select(name,mass,strain)
sigmatches = left_join(sigmatches,prism_structures,by=c("Name"="name"))

eics = lapply(1:nrow(sigmatches),function(idx) {
  specFile = str_replace(sigmatches$SpecFile[idx],pattern = "MS2",replacement = "MS1")
  mz_target = as.numeric(sigmatches$mass[idx])+1.00784
  mz = c( mz_target - (mz_target * 20 / 10^6),
          mz_target + (mz_target * 20 / 10^6))
  chromatogram(readMSData(specFile,msLevel. = 1,mode = "onDisk"),
                                           mz=mz,missing=1)})


svg("figures/prism_matches/eics.svg",width = 15)
par(mfrow=c(3,6))
for(idx in 1:18) plot(eics[[idx]],main=paste(sigmatches$Name[idx],"\nm/z:",
                                             as.numeric(sigmatches$mass[idx])+1.00784,sep = " "),
                      col="black")
dev.off()

rbanimgs = Sys.glob("figures/prism_matches/compounds/*.png") %>% sort()

plots <- lapply(rbanimgs,function(x){
  img <- as.raster(readPNG(x))
  grid.arrange(rasterGrob(img, interpolate = FALSE),top=str_extract(x,"G(.)*s\\d+"),nrow = 1)
})

svg("figures/prism_matches/rban_imgs_combined.svg",width = 15)
grid.arrange(grobs = plots, nrow=3, ncol=6,top=NULL)
dev.off()


## selected MS/MS
spectrum =  read_csv("figures/prism_matches/GA13-004_c14_s3_ms2_40v_3.3min.csv") %>% janitor::clean_names()
spectrum = new("Spectrum2",mz=spectrum$m_z,intensity=spectrum$abund,centroided=T)
svg("figures/prism_matches/GA13-004_c14_s3_ms2_40v_3.3min.svg")
plot(spectrum,"AVKTV",main="non-ribosomal peptide AVKTV at 3.3min\nprecursor=515.3326,collision energy=40",
     col="black")
dev.off()