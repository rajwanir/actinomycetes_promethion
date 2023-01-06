library(tidyverse)
library(GenomicRanges)

#prokka annotation
gff = Sys.glob(file.path(getwd(), "data/prokka/*/*.gff"))
gff = lapply(gff,function(x) as_tibble(rtracklayer::readGFF(x)) )%>% setNames(.,gff)
gff = bind_rows(gff,.id="gff_path")
gff = mutate(gff,strain=str_extract(gff_path,pattern = "\\w+-\\w+"))


antismash = Sys.glob(file.path(getwd(), "data/antismash/*/*.gff.1"))
antismash = lapply(antismash,function(x) rtracklayer::readGFF(x)) %>% setNames(.,antismash)
antismash = bind_rows(antismash,.id="gff_path")
antismash = mutate(antismash,strain=str_extract(gff_path,pattern = "\\w+-\\w+"))
bgcs = antismash %>% filter(type...3 == "sequence_feature") %>% 
  filter(!is.na(protocluster_number))


# filter genes that overlap antismash BGC
gff$bgc_overlap = findOverlaps(makeGRangesFromDataFrame(gff %>% mutate(seqid = paste(seqid,strain,sep = "-")),keep.extra.columns=T),
             makeGRangesFromDataFrame(bgcs %>% mutate(seqid = paste(seqid,strain,sep = "-")),keep.extra.columns=T),select="first")
gff = gff %>% mutate(bgc_overlap_l = if_else(!is.na(bgc_overlap),T,F))
gff = gff %>% filter(bgc_overlap_l==T)


inner_join(antismash,
           gff %>% filter(str_detect(product,pattern = "Glutamyl-")),
           by=c("start","strain","end")) %>% select(product.x,product.y,strain,locus_tag.x,start)

# bgcs = filter(gff,type...3 == "sequence_feature" & !is.na(protocluster_number))
# gff %>% filter(gene_kind == "regulatory" &
#                  !str_detect(gene_functions,pattern = "transcriptional regulator|kinase|transcription regulator|Transcription regulator|response regulator|repressor|cold-shock")) %>%
#   select(gene_functions)
# smcogs = c("SMCOG1032","SMCOG1053","")
# gff %>% filter(gene_kind == "regulatory" & 
#                  !str_detect(gene_functions,pattern = "transcriptional regulator|kinase|transcription regulator|Transcription regulator|response regulator|repressor|cold-shock")) %>%
#   select(gene_functions)

### GNPS network analogs ####

antibacterial = read_tsv("../databases/antibacterial_pubchem_cids.tsv")
mibig_entry = mibig_entry %>% filter(str_detect(chem_acts,pattern = "Antibacterial") & !is.na(chem_struct))
mibig_entry$inchikey = unlist(sapply(mibig_entry$chem_struct,get.inchi.key,USE.NAMES = F))
mibig_entry = mibig_entry[sapply(mibig_entry$inchikey,function(x) !is_null(x)),]
mibig_entry$inchikey = unlist(mibig_entry$inchikey)


#for promethion data
gnps_hits = left_join(read_tsv(Sys.glob("data/QTOF/db_screens/gnps_network/METABOLOMICS-SNETS-V2-*-group_by_compound-main.tsv")),
 read_tsv(Sys.glob("data/QTOF/db_screens/gnps_network/clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.clustersummary")) %>% select(UniqueFileSources,SpectrumID),
 by="SpectrumID")

#for flongle data
gnps_hits = left_join(read_tsv("../soil_metagenomics/data/QTOF/gnps_network/METABOLOMICS-SNETS-V2-f1c83645-group_by_compound-main.tsv"),
                      read_tsv("../soil_metagenomics/data/QTOF/gnps_network/clusterinfosummarygroup_attributes_withIDs_withcomponentID/b75da5998fa148cdb62e1bce55ab249a.clustersummary") %>% select(UniqueFileSources,SpectrumID),
                      by="SpectrumID")

# filtering antimicrobial compounds from gnps network
gnps_hits_antibacterial =  gnps_hits %>% filter(`InChIKey-Planar` %in% str_extract(antibacterial$inchikey,pattern = "[A-Z]+" ) & MZErrorPPM >50) %>% 
  select(Compound_Name,npclassifier_class,UniqueFileSources,contains("MZ"),ExactMass,Adduct,SharedPeaks,MQScore,`InChIKey-Planar`) %>% mutate(specmass=SpecMZ-1.00784 )


require(httr)
urls = sprintf("https://www.npatlas.org/api/v1/compounds/massSearch?mass=%s&type=exact_mass&range=20&unit=ppm&skip=0&limit=10",gnps_hits_antibacterial$specmass)
npatatlas_out = lapply(urls,function(url) content(POST(url,content_type("application/json"))) )
gnps_hits_antibacterial$npatlas_massMatch = sapply(npatatlas_out,function(x) paste(sapply(x,function(y) y$original_name),collapse = ",") )


write_tsv(gnps_hits_antibacterial,"data/QTOF/db_screens/gnps_network/imgs/gnps_hits_antibacterial.tsv")



