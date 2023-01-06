



barcode_key = read_csv("data/flongleQC/plate2_flongle_06092021/barcode_key.csv")

plate_positions = as.vector(outer(LETTERS[1:8],1:9,paste0)) %>% sort()
  qtof_worklists = barcode_key %>% mutate(
    culture_position = paste0(row, col),
    qtof_position = case_when(col<=9&!row %in% c("G","H") ~ paste0("P1",culture_position)),
    method = "SO_10mingradient_pos_v3.m",
    sample_name = paste(strain, culture_position, sep =
                          "_"),
    data_file = sprintf(
      "D:\\MassHunter\\Data\\Rahim\\promethion\\plate2\\ISP1_MEOH_%s.d",
      sample_name
    )
  ) 
  
  qtof_worklists[is.na(qtof_worklists$qtof_position),]$qtof_position = paste0("P2",qtof_worklists[is.na(qtof_worklists$qtof_position),]$row,1:3)
  
  qtof_worklists = qtof_worklists %>% select(sample_name, qtof_position, method, data_file)
  
  
  write_csv(qtof_worklists, file = "data/QTOF/worklists/plate2_isp1_10212021.csv")
