
domain = "Pkinase"
class = "iv"
  
target_proteins = bind_rows(precursorRegions) %>% 
  filter(gene_kind == "biosynthetic" & str_detect(sec_met_domain,domain) & str_detect(gene_functions,
                                                                                      sprintf("class-%s",class)) ) %>% 
  mutate(name=paste(strain,locus_tag,sep='@')) %>% distinct(name,translation) 
target_proteins = Biostrings::AAStringSet(setNames(target_proteins$translation,target_proteins$name),use.names = T)
target_proteins = muscle::muscle(target_proteins)

Biostrings::writeXStringSet(target_proteins@unmasked,
                            sprintf("data/ripps/primers/alignments/%s-class%s.fasta",domain,class))
ape::write.tree(ape::nj(ape::dist.aa(ape::as.AAbin(target_proteins))), sprintf("data/ripps/primers/alignments/%s-class%s.nwk",domain,class))
 # select(strain,locus_tag,translation,sec_met_domain,gene_functions)



library(Biostrings)
library(tidyverse)


### extracting nucl sequences


gff=read_csv("data/ripps/primers/sequenced_strains/micKC_complete_metadata.csv") %>% mutate(cds_id=paste(strain,locus_tag,sep="@"))
# gff = mutate(gff,assm=str_replace(str_replace(gff_path,"antismash","assemblies/homopolish"),".gff.1",".fasta" ) )
# gff= filter(gff,!str_detect(seqid,"c"))


genomes = lapply(unique(gff$assm),function(x) readDNAStringSet(x, use.names = T))
for (gi in 1:length(genomes)){
  names(genomes[[gi]]) = paste(unique(gff$strain)[gi], names(genomes[[gi]]), sep = "@")}
genomes = unlist(DNAStringSetList(genomes))


CDS = sapply(1:nrow(gff), function(idx) {
  print(idx)
  genome = readDNAStringSet(gff$assm[idx], use.names = T)
  names(genome) = str_remove_all(names(genome)," .*")
  BSgenome::getSeq(
    genome,
    GenomicRanges::makeGRangesFromDataFrame(gff[idx, ])
  )
})
CDS = unlist(DNAStringSetList(CDS))
names(CDS)=gff$cds_id

coords = (str_locate(gff$translation,"HG") * 3) %>% as_data_frame() %>% mutate(seqid = gff$cds_id)
BSgenome::getSeq(x=CDS,names=GenomicRanges::makeGRangesFromDataFrame(coords %>% filter(!is.na(start)) )) %>% 
  consensusString()
sapply(CDS,function(x) matchPattern(pattern="RCCSGCMCTSTCVCCCGGMTTCTTC",subject=x,fixed=F,max.mismatch=1))
sapply(CDS,function(x) matchLRPatterns(Lpattern = "CGTSCTSAAGGARGCCCGKCCG",Rpattern="CCGSCTCTTCCCSGGMGAC",subject=x,max.Lmismatch=2,max.Rmismatch=2,Lfixed=F,Rfixed = F,max.gaplength=3000)) 


# pfam_scan=read_csv("data/ripps/primers/sequenced_strains/LANC_like_pfamscan.csv") %>% janitor::clean_names()
# pfam_scan=left_join(pfam_scan,gff,by=c("seq_id"="cds_id"))
# pfam_scan=mutate(pfam_scan,
#                          start=case_when(strand=="+" ~ start+(alignment_start*3), 
#                                          strand=="-" ~ end-(alignment_end*3)),
#                          end=case_when(strand=="+" ~ start+(alignment_end*3),
#                                        strand=="-" ~ end-(alignment_start*3)
#                                        ))
# gff$domain_nucl=as(CDS,"character")
# gff$domain_nucl_rvcomp=as(setNames(reverseComplement(CDS),gff$seq_id)
#                                 ,"character")
# mutate(gff,domain_nucl = case_when(strand=="+"~domain_nucl,strand=="-"~domain_nucl_rvcomp))
# 
# writeXStringSet(setNames(as(gff$domain_nucl,"DNAStringSet"),gff$seq_id),
#                          "data/ripps/primers/sequenced_strains/LANC_like_domains_nucl.fasta")



gff=read_csv("data/ripps/primers/sequenced_strains/lanB_metadata.csv") %>% mutate(cds_id=paste(strain,locus_tag,sep="@"))


#################

library(DegeneratePrimerTools)
library(Biostrings)
library(tidyverse)
# pfamdna <- retrieve_PFAM_nucleotide_sequences("PF00733",retrievaltype="rest")

# write_tsv(pfamdna,"data/ripps/primers/pfam/PF04738_nucl.tsv")
# pfamdna=sample_n(pfamdna,600)
# pfamdna = filter(pfamdna,!is.na(domainsequence))
# pfamdnaset = Biostrings::DNAStringSet(setNames(pfamdna$domainsequence,pfamdna$EMBL_ID),use.names = T)

pfam = "PF13575"
pfamdna = readDNAStringSet(Sys.glob(sprintf("data/ripps/primers/pfam/%s_nucl.fasta",pfam)))

pfamdegeprime = degeprimer(pfamdna)
pfamdegeprime@msa = Biostrings::readDNAMultipleAlignment(Sys.glob(sprintf("data/ripps/primers/pfam/%s_nucl.fasta.aligned.fasta",pfam)))
pfamdegeprime@phy_tree = ggtree::read.tree(Sys.glob(sprintf("data/ripps/primers/pfam/%s_nucl.tree",pfam))) 
pfamdegeprime=pfamdegeprime %>% design_primers(maxdegeneracies=c(1, 10, 20), number_iterations=2)


