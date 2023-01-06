

# sapply(1733:1742,function(idx) 
#   any(str_detect(names(ripps_seq[str_detect(ripps_seq,matches$FragmentSeq[idx])]),matches$spec_strain[idx])))

monomers = lapply(Sys.glob("data/QTOF/db_screens/nerpa/*/structures.info"),function(x) read_delim(x,delim = " ",col_names = c("npaid","monomers","extra")))

# for getting npclassifier from API #####
# smile = "CC[C@H](C)C[C@H](C)CCCCCCCCC(=O)N[C@@H]1CCCNC(=O)[C@@H]2[C@H]([C@H](CN2C(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H]3C[C@H](CN3C(=O)[C@@H](NC1=O)[C@@H](C)O)O)[C@@H]([C@H](C4=CC=C(C=C4)O)O)O)[C@@H](CC(=O)N)O)C)O"
# 
# npatlas_matches$npclassifier = sapply(npatlas_matches$SMILES,function(s){
#   print(which(npatlas_matches$SMILES == s))
#   tryCatch(jsonlite::fromJSON(sprintf("https://npclassifier.ucsd.edu/classify?smiles=%s",s))[["class_results"]],
#   error=function(e) return("unknown") )
# }
# ,USE.NAMES = FALSE)
# 
# 
# npatlas_matches$npclassifier = sapply(npatlas_matches$npclassifier,function(x) paste(x,collapse = ","))


npatlas_matches = lapply(Sys.glob("data/QTOF/db_screens/npatlas/*/significant_matches.tsv"),
                         read_tsv) %>% bind_rows()
npatlas_ripps = read_tsv("../databases/NPAtlas_npclassifier_09082022.tsv") %>% filter(str_detect(npatlas_class,"RiPPs"))
npatlas_sequence = read_csv("data/ripps/npatlas_ripps_sequences.csv")
npatlas_matches = inner_join(npatlas_matches,npatlas_ripps,by=c("Name"="npaid"))

npatlas_sequence = read_csv("data/ripps/npatlas_ripps_sequences.csv")
npatlas_matches = left_join(npatlas_matches,select(npatlas_sequence,-compound_names))

npatlas = read_tsv("../databases/NPatlas_01132020_inchikey.tsv")
npatlas_matches = left_join(npatlas_matches,npatlas,by=c("Name"="npaid"))

npatlas_matches = npatlas_matches %>% mutate(strain = str_extract(SpecFile,"G?[A-Z]\\d+-\\d+"))
npatlas_matches = npatlas_matches %>% distinct(strain,Name,.keep_all=T)

# checking in assembly #########
check_pepseq_in_assm = function(strain,query_peptide) {
# strain = "GA6-010"
# query_peptide="ITSVSWCTPGCTSEGGGSGCSHCC"
assm = readDNAStringSet(sprintf("data/assemblies/homopolish/%s/%s_homopolished.fasta",strain,strain))

translation = lapply(1:3, function(pos) 
  subseq(c(assm, reverseComplement(assm)), start=pos)) %>% 
  lapply(., translate)

n_matches = sapply(translation,function(x) vcountPattern(subject = x,pattern = query_peptide,max.mismatch = floor(nchar(query_peptide) * 0.1 ))) %>% 
  sum()

return(n_matches)
}

ripps_matches = npatlas_matches %>% filter(!is.na(sequence))
ripps_matches$assm_seq_matches = sapply(1:nrow(ripps_matches), function(i) check_pepseq_in_assm(strain=ripps_matches$strain[i],query_peptide=ripps_matches$sequence[i]))

ripps_matches %>% filter(between(assm_seq_matches,1,10))



# rippminer core check in assembly ######
rippminer = read_csv("../tools/rippminer_standalone/DB/peptide_info.csv")
rippminer = rippminer %>% filter(!is.na(sequence)) %>% 
  mutate(sequence = case_when(!str_detect(sequence,"^M")  ~ sequence,
                              str_detect(sequence,"^M")  ~ substr(sequence,nchar(sequence)-14, nchar(sequence))  ))

assm_paths=Sys.glob("data/assemblies/homopolish/*/*_homopolished.fasta")
assm_strain=str_extract(assm_paths,"G?[A-Z]\\d+-\\d+")
assm = lapply(assm_paths,readDNAStringSet)
for (gi in 1:length(assm)){
  names(assm[[gi]]) = paste(assm_strain[gi], names(assm[[gi]]), sep = "@")}
assm = unlist(DNAStringSetList(assm))


translation = lapply(1:3, function(pos) 
  subseq(c(assm, reverseComplement(assm)), start=pos)) %>% 
  lapply(., translate)

rippminer$assm_seq_matches = mcsapply(rippminer$sequence,function(query_peptide)
  sapply(translation,function(x) vcountPattern(subject = x,pattern = query_peptide,max.mismatch = floor(nchar(query_peptide) * 0.3 ))) %>% 
  sum()
,mc.cores=50)


rippminer$assm_seq_matches = sapply(rippminer$sequence,function(query_peptide)
  sapply(translation,function(x) vcountPattern(subject = x,pattern = query_peptide,max.mismatch = floor(nchar(query_peptide) * 0.1 ))) %>% 
    sum()
)

write_tsv(rippminer,"data/ripps/rippminer_hits.tsv")

# mibig core check in assembly ######
mibig_prot = readAAStringSet("../databases/MiBiG/mibig_all_prot_seqs/mibig_prot_seqs_2.0.fasta")
mibig_prot = mibig_prot[width(mibig_prot)<100 & width(mibig_prot)>10]
mibig_ripps = ripps = read_csv("../databases/MiBiG/ripp_val/ripps_id.csv")
mibig_prot = mibig_prot[str_detect(names(mibig_prot),paste(mibig_ripps$mibig_id,collapse = "|"))]
mibig_prot = Biostrings::subseq(mibig_prot,start=-15) # pattern should be smaller than the subject

mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}


mibig_matches = mcsapply(mibig_prot,function(query_peptide)
  sapply(translation,function(x) vcountPattern(subject = x,pattern = query_peptide,max.mismatch = floor(nchar(query_peptide) * 0.3 ))) %>% 
    sum()
  ,mc.cores=50)

mibig_matches = stack(mibig_matches)

write_tsv(mibig_matches,"data/ripps/ripps_mibig_hits.tsv")


  
match_key = sapply(matches$sequence,function(x) vcountPattern(subject = ripps,pattern = x,max.mismatch = floor(nchar(x) * 0.2 )))
sapply(1:nrow(matches),function(i) ripps[sapply(match_key[,i],as.logical)])

#retriving core_peptides
match_key= sapply(translation,function(x) vcountPattern(subject = x,pattern = query_peptide,max.mismatch = floor(nchar(query_peptide) * 0.3 ))) 
matches_translation = sapply(1:3, function(i) translation[[i]][sapply(match_key[,i],as.logical)]) 
matches_translation_coords = sapply(matches_translation,function(x) vmatchPattern(subject = x,pattern = query_peptide,max.mismatch = floor(nchar(query_peptide) * 0.3 )))
lapply(1:2, function(i) sapply(1:length(matches_translation[[i]]),function(j)  matches_translation[[i]][[j]][matches_translation_coords[[i]][[j]]]))



translate(subseq(assm,end=1305661*3,width=nchar("CLGIGSCNDFAGCGYAIVCFW")*3)) # traslation coords to assm