#!/usr/bin/env Rscript

library(stringi)

#### Calculate degenerated oligonucleotide frequencies ####

files <- commandArgs(TRUE) # pass bash input to R

for(file in files){
  # Read raw oligonucleotide frequency matrix (e.g. output of `ONF_calc.R`)
  if(grepl("\\.rds$", file, ignore.case = TRUE)){
    data <- readRDS(file)
  } else {
    data <- data.matrix(read.delim(file, row.names = 1))
    saveRDS(data, sub("\\.\\w*$", ".rds", file)) # save ONF matrix as RDS file if it isn't already
  }
  
  # Find column indices of reverse complement for each kmer
  kmers <- colnames(data)
  dgen.key <- stri_replace_all(kmers, fixed = c("a","c","g","t"), c("T","G","C","A"), vectorize_all = FALSE)
  dgen.key <- stri_reverse(tolower(dgen.key))
  dgen.key <- match(dgen.key, kmers)
  cols2sum <- !is.na(dgen.key) & dgen.key > 1:length(kmers)
  cols2del <- !is.na(dgen.key) & dgen.key < 1:length(kmers)
  
  data[,cols2sum] <-  data[,cols2sum] + data[,dgen.key[cols2sum]] # sum kmers + reverse complements
  data <- data[,!cols2del] # remove reverse complements from matrix
  # do nothing to kmers that are their own reverse complement
  data <- data/rowSums(data) # normalize degenerated kmer counts
  
  # Save degenerated ONF matrix as RDS file
  saveRDS(data, sub("\\.\\w*$", "_dgen.rds", file))
}