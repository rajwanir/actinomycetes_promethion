library(MSnbase)
library(tidyverse)

npatlas = read_tsv("../databases/NPatlas_01132020_inchikey.tsv")
sigmatches = lapply(Sys.glob("data/QTOF/db_screens/npatlas/*/significant_unique_matches.tsv"),read_tsv) %>% bind_rows()
sigmatches = left_join(sigmatches,npatlas,
                       by=c("Name"="npaid"))

gnps_matches = sigmatches %>% filter(gnps_ids != "[]") %>%
  distinct(Name, .keep_all = T) %>% mutate(gnps_id_clean =
                                             str_extract(gnps_ids, pattern = "CCMSLIB\\d+"))
gnps_matches = gnps_matches %>% filter(!gnps_id_clean %in% c("CCMSLIB00004689471","CCMSLIB00000839188",
                                                             "CCMSLIB00004688873","CCMSLIB00004688569"))

loadGNPSSpectra = function(gnps_id){
  #retreives spectra from GNPS and returns as msnbase object
  #example: gnps_id = "CCMSLIB00004722223"
  print(gnps_id)
  json_file = sprintf("https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID=%s",
                      gnps_id)
  json_file =  tryCatch(jsonlite::fromJSON(json_file),error = function(e) NULL)
  peak_list = json_file$spectruminfo$peaks_json %>% jsonlite::fromJSON(flatten = T) 
  GNPSspectra = new("Spectrum2",
                    mz = as.numeric(peak_list[,1]),
                    intensity = as.numeric(peak_list[,2]),
                    precursorMz = as.numeric(sort(json_file$annotations$Precursor_MZ,decreasing = T)[1]),
                    precursorCharge = as.integer(json_file$annotations$Charge)
  )
  
  GNPSspectra@centroided = T
  return(GNPSspectra)
}


gnpsSpectra = lapply(gnps_matches$gnps_id_clean,loadGNPSSpectra)

sampleSpectra = lapply(1:nrow(gnps_matches),function(idx)
  spectra(filterPrecursorScan(readMSData(gnps_matches$SpecFile[idx],mode="onDisk"),
                      gnps_matches$Scan[idx]))[[2]]
)

gnps_matches$gnps_precmz = sapply(gnpsSpectra,precursorMz)
gnps_matches = gnps_matches %>% mutate(gnps_mass_diff =  abs(SpectrumMass - gnps_precmz )) 


gnpsSpectra = gnpsSpectra[gnps_matches$gnps_mass_diff < 0.1]
sampleSpectra = sampleSpectra[gnps_matches$gnps_mass_diff < 0.1]
gnps_matches = gnps_matches %>% filter(gnps_mass_diff < 0.1)

gnps_matches$gnps_common_peaks = sapply(1:nrow(gnps_matches),function(idx) compareSpectra(sampleSpectra[[idx]],gnpsSpectra[[idx]],fun = "common"))
gnps_matches$gnps_dotprot = sapply(1:nrow(gnps_matches),function(idx) compareSpectra(sampleSpectra[[idx]],gnpsSpectra[[idx]],fun = "dotproduct"))

gnps_matches = gnps_matches %>% mutate(strain = str_extract(SpecFile,pattern = "G\\w+-\\d+"))


svg("figures/gnps_alignment_49_58.svg",width = 11, height = 8)
par(mfrow=c(4,4))
for (idx in 49:58) plot(sampleSpectra[[idx]],gnpsSpectra[[idx]],
                       legend.cex = 0.0001,
                       main = sprintf("%s in %s\n[%d common]",
                                      gnps_matches$Name[idx],gnps_matches$strain[idx],gnps_matches$gnps_common_peaks[idx]),
                       cex.main=0.9)

dev.off()

#?`plot,Spectrum,Spectrum-method`
svg("figures/gnps_alignments/commonpeaks_hist.svg")
hist(gnps_matches$gnps_common_peaks,breaks = 100,col="brown",xlab = "# of peaks common",
     main=NULL)
dev.off()
svg("figures/gnps_alignments/dotproduct_hist.svg")
hist(gnps_matches$gnps_dotprot,col="brown",xlab = "Dot product",
     main=NULL)
dev.off()