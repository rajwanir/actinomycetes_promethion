library(tidyverse)

nerpa_paths = Sys.glob("data/QTOF/db_screens/nerpa/*/report.csv")
nerpa = lapply(nerpa_paths,function(x) read_csv(x,col_types = cols(.default = "c")))