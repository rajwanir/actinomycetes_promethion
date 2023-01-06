




sequence="AAAASCNCVCGFCCSCSPSA"

dehydration_mass = 18.0105 # either cyclodehydration or not
az_reduction_mass = 2.01565007 #-H2
N_ter_rem=16.018724 # -NH2
pyridine_mass = (18.0105*3) + N_ter_rem # assuming one is dehydrated already
thioether_mass=2.015650 #-H2
proton=1.00784 

cores = sapply(5:50,function(i) str_sub(sequence,start=-i))

core_tbl=lapply(cores,function(core){
n_sites = str_count(core,"S|T|C")
dehydrated_mass=sapply(0:n_sites, function(i) Peptides::mw(core,monoisotopic = T) - (i*dehydration_mass))
az_reduced_mass=sapply(0:(length(dehydrated_mass)-1), function(i) dehydrated_mass - (i*az_reduction_mass) )

tibble(
  core_seq=core,
  dehydration=rep(0:n_sites,n_sites+1),
  az_reduction=as.vector(sapply(0:n_sites,function(s) rep(s,n_sites+1))),
  mass=as.vector(az_reduced_mass)
)

})

core_tbl = bind_rows(core_tbl)

core_tbl=mutate(core_tbl,
                pyridine=case_when(dehydration>0 ~ 1,
                                   TRUE ~ 0),
                mass=mass-(pyridine_mass*pyridine))

core_tbl=bind_rows(core_tbl,
              mutate(core_tbl,thioether=case_when(pyridine>0 ~ 1,
                                            TRUE ~ 0),
                mass=mass-(thioether_mass*thioether))
)

core_tbl=mutate(core_tbl,mass=mass+proton,proton=1)
