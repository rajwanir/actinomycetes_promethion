assm_paths=Sys.glob("data/assemblies/homopolish/GA6-010/GA6-010_homopolished.fasta")
assm_strain=str_extract(assm_paths,"G?[A-Z]\\d+-\\d+")
assm = lapply(assm_paths,readDNAStringSet)
for (gi in 1:length(assm)){
  names(assm[[gi]]) = paste(assm_strain[gi], names(assm[[gi]]), sep = "@")}
assm = unlist(DNAStringSetList(assm))
assm=assm[1]
assm=subseq(assm,start=3695916,end=3717750)

seq=assm
orfs <- findORFs(seq, longestORF = T)
gr <- unlist(orfs, use.names = TRUE)
gr <- GRanges(seqnames = names(seq)[as.integer(names(gr))],
              ranges(gr), strand = "+")
names(gr) <- paste0("ORF_", seq.int(length(gr)), "_", seqnames(gr))
orf_seqs <- getSeq(seq, gr)
orf_seqs