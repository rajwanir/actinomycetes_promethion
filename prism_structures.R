library(tidyverse)
library(rcdk)
library(parallel)
#functions
getStructures = function(json_file = x) {
  prismdf1 = jsonlite::fromJSON(json_file,flatten = T,simplifyDataFrame =T)$prism_results$clusters
  
  ##get the main products #########  
  #list of dfs, each df corresponds to cluster
  main_products = prismdf1$predicted_molecule_masses
  
  if (length(main_products) != 0 & !is.null(unlist(main_products)) ) {
    # append the rownames to include cluster and sample idx
    for (idx in seq(1, length(main_products), 1)) {
      if (nrow(main_products[[idx]]) != 0) {
        #smiles is a list
        main_products[[idx]] = tibble(smiles = unlist(main_products[[idx]]$smiles))
        
        
        #adding cluster and smiles number
        main_products[[idx]] = mutate(
          main_products[[idx]],
          BGC_no = as.character(idx),
          smile_no = rownames(main_products[[idx]]),
          name = sprintf("c%s_s%s", BGC_no, smile_no)
        )
        
        
      }
    }
    #single df per sample including all clusters
    main_products = bind_rows(main_products)
  }
  
  
  ## get the intermediate products #####
  
  intermediate_smiles = lapply(prismdf1$biosynthetic_pathways,function(biosynpathway) biosynpathway$pathway %>% lapply(.,function(pathway) pathway$intermediate_smiles)) 
  
  if (length(intermediate_smiles)!= 0 & !is.null(unlist(intermediate_smiles))) {
  names(intermediate_smiles) = as.character(1:length(intermediate_smiles))
  intermediate_smiles = lapply(intermediate_smiles, function(ismiles) {
    #there are empty lists with no structures so check first
    if (length(ismiles) > 0) {
      if(!is.null(ismiles[[1]])) {
        #add names to list for bind_rows
        names(ismiles) = as.character(1:length(ismiles))
        #bind_rows
        ismiles = bind_rows(ismiles)
        ismiles$intermediate_no = 1:nrow(ismiles)
        ismiles = pivot_longer(
          data = ismiles,
          cols = where(is.character),
          names_to = "smile_no",
          values_to = "smiles"
        )
      }
    }
    return(ismiles)
  }) %>%
    bind_rows(., .id = "BGC_no")
  intermediate_smiles = mutate(intermediate_smiles, name = sprintf("c%s_i%s_s%s", BGC_no, intermediate_no, smile_no))
  
  ## combine main and intermediates 
  all_structures = bind_rows(main_products,intermediate_smiles)
  } else(all_structures = main_products)
  
  ### refine ####
  
  if((length(main_products) != 0 & !is.null(unlist(main_products))) | (length(intermediate_smiles) != 0) & !is.null(unlist(intermediate_smiles))){
  #use rcdk to disconnect molecules, take larger of two and then calc exact mass of the disconnected mols
  allSmiles = parse.smiles(all_structures$smiles)
  allSmiles = lapply(allSmiles,get.largest.component)
  all_structures$mass = as.numeric(sapply(allSmiles,get.exact.mass))
  all_structures$smiles = as.character(sapply(allSmiles,get.smiles))
  all_structures$n_aminoacids = 0
  return(all_structures)} else{return(NULL)}
  
  
  }




prism_jsons = Sys.glob("data/prism/*.json")




cl <- makeCluster(30)
clusterEvalQ(cl, library("tidyverse")) 
clusterEvalQ(cl, library("rcdk")) 
clusterExport(cl, "prism_jsons")
clusterExport(cl, "getStructures")

structures = parLapply(cl,prism_jsons,function(x) {
                  print(x)
                  getStructures(json_file = x)})

names(structures) = prism_jsons
structures = bind_rows(structures,.id="prism_path")

write_tsv(structures,"data/prism/prism_structures.tsv")

structures$inchikey = sapply(structures$smiles,get.inchi.key)



# split db

structures = read_tsv("data/prism/prism_structures.tsv")

structures = structures %>% filter(mass > 200)
structures = structures %>% mutate(strain = str_extract(prism_path,pattern = "\\w+-\\d+"),
                                   name = paste(strain,name,sep = "_"),
                                   cmds = sprintf("obabel --gen3D -:\"%s\" -omol -x3 > data/QTOF/databases/%s/mols/%s.mol",
                                                  smiles,strain,name),
                                   molfile = paste0("mols/",name,".mol"),
                                   metadata = "NA")

lapply(unique(structures$strain),function(given_strain){

dir.create(sprintf("data/QTOF/databases/%s/mols",given_strain),recursive=T)
  
strain_structures=structures %>% filter(strain==given_strain)

write.table(strain_structures[c('molfile','name','mass',
                             'n_aminoacids','metadata')] ,
            file = sprintf("data/QTOF/databases/%s/library.info",given_strain),
            col.names = F,
            row.names = F,
            quote = F)

write.table(strain_structures['smiles'],
            file = sprintf("data/QTOF/databases/%s/library.smiles",given_strain),
            col.names = F,
            row.names = F,
            quote = F)

})



write.table(structures$cmds,
            file = "src/swarms/mol_files_obabel.swarm",
            quote = F,row.names = F,col.names = F)

# construct a finger print
lapply(unique(structures$strain),function(given_strain){
  db_folder = sprintf("data/QTOF/databases/%s",given_strain)
  cmd = sprintf("cut -d' ' -f2 %s/library.info | paste %s/library.smiles - | obabel -ismi -osmi > %s/smiles.smi",
                db_folder,db_folder,db_folder)
  system(cmd)
  return(NULL)
})

lapply(unique(structures$strain),function(given_strain){
  db_folder = sprintf("data/QTOF/databases/%s",given_strain)
  cmd = sprintf("obabel -ismi %s/smiles.smi -ofs -xfMACCS",
                db_folder,db_folder,db_folder)
  system(cmd)
  return(NULL)
})