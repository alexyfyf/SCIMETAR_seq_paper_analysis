pcr2 <- read.csv('output/LINE1_qPCR_all.csv')

test <- pcr2 %>% select(Plate_Well.x, `Amplicon_Barcode_5._Read1`, `Amplicon_Barcode_5._Read2`) %>% 
  rename(name = 'Plate_Well.x',
         forward = 'Amplicon_Barcode_5._Read1',
         reverse = 'Amplicon_Barcode_5._Read2') 

test %>%
  mutate(
    forward = str_remove_all(forward, "^H|VS$|HW$|V$|H$"),
    reverse = str_remove_all(reverse, "^H|VS$|HW$|V$|H$")
  ) %>% 
  # mutate(
  #   forward_len = str_length(forward),
  #   reverse_len = str_length(reverse)
  # ) %>%
  # select(forward_len, reverse_len) %>%
  # summary %>%
  write.table('barcode_renamed.tsv', quote = F, row.names = F, sep = '\t')
  
