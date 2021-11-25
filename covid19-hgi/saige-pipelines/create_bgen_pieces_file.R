rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)

setwd('~/Projects/covid19-hgi/release_15042021/')

pieces <- fread('bgen_pieces.txt', header = F) %>% 
  mutate(chr_piece = gsub(".bgen", "", gsub(".*_","",V1))) %>% 
  separate(chr_piece, c("chr", "piece"), convert = TRUE) %>%
  mutate(chr = as.integer(gsub("X", "23", chr))) %>% 
  arrange(chr, piece) %>% 
  pull(V1)

N_per_line <- 16
multi_pieces <- NULL

i <- 1
while (i <= length(pieces)) {
  multi_pieces <- c(multi_pieces, paste(pieces[i:(i+N_per_line-1)], collapse = "\t"))
  
  i <- i+N_per_line

}

# remove NA and extra tabs added at the end of the last line (needs to manually count how many)
# TODO: find some better way for this

multi_pieces <- gsub("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA$","",multi_pieces)

fwrite(list(multi_pieces), paste0("bgen_pieces_", N_per_line, ".txt"), col.names = FALSE, na = "")
system("gsutil cp bgen_pieces_16.txt gs://dsge-covid19-data/15042021/conf/")

fwrite(list(multi_pieces[1]), paste0("bgen_pieces_", N_per_line, "_test.txt"), col.names = FALSE, na = "")
system("gsutil cp bgen_pieces_16_test.txt gs://dsge-covid19-data/15042021/conf/")
