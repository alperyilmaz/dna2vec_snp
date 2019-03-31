#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lsa))
suppressPackageStartupMessages(library(optparse))

# parsing file name

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# import word2vector data

write("** importing dna2vec data..", stderr())

w2v <- read_delim(opt$file, 
                  delim = " ", 
                  skip = 1, 
                  col_names = FALSE)

# keep vector nested, now there are only two columns

write("** nesting vector data..", stderr())

w2v_nested <- w2v %>% 
  rename(word=X1) %>% 
  group_by(word) %>% 
  nest() %>% 
  mutate(vec = map(data,~ as.numeric(.x[1,]))) %>% 
  select(-data)

# function to replace nth nucleotide with given nucleotide

replace_nt <- function(word,coord,alt){
  str_sub(word,coord,coord) <- alt
  word
}

# tidy way of calculating neighbors (hamming distance 1)
# needs improvement, kind of slow

write("** generating neighbors list..", stderr())

eightmers_neighbors <- w2v %>% 
  filter(str_length(X1)==8) %>% 
  select(word=X1) %>% 
  mutate(coordinate=map(word,~ seq(nchar(.x)))) %>% 
  unnest(coordinate) %>% 
  mutate(original=str_sub(word,coordinate,coordinate)) %>% 
  mutate(variant=list(c("A","C","G","T"))) %>% 
  unnest(variant) %>% 
  filter(original != variant) %>% 
  mutate(var_word=pmap_chr(list(word,coordinate,variant), ~ replace_nt(..1,..2,..3)) )


# join vector data, then calculate cos_sim for every pair

write("** calculating cosine similarity", stderr())

result <- eightmers_neighbors %>% 
  #unnest(neighbor) %>% 
  #filter(word != neighbor) %>% 
  left_join(w2v_nested) %>% 
  left_join(w2v_nested, by=c("var_word"="word")) %>% 
  mutate(cos_sim=map2_dbl(vec.x, vec.y, ~cosine(.x, .y))) %>%
  select(word,coordinate,original,variant,var_word,cos_sim)

# write to stdout

write("** writing results to output..", stderr())

cat(format_csv(result))
