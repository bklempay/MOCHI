#!/usr/bin/env Rscript

library(stringi)

#### Calculate degenerated oligonucleotide frequencies ####

files <- commandArgs(TRUE) # pass command line arguments to R

for (f in files) {
  if (!(file.exists(f) && grepl("\\.(rds|tsv)$", f, ignore.case = TRUE))) {
    stop("illegal file '", f, "'")
  }
  
  # Read raw oligonucleotide frequency matrix (e.g. output of `ONF_calc.R`)
  if (grepl("\\.rds$", f, ignore.case = TRUE)) onf.matrix <- readRDS(f)
  if (grepl("\\.tsv$", f, ignore.case = TRUE)) {
    onf.matrix <- data.matrix(read.delim(f, row.names = 1))
    # assumes row names in first column
  }
  
  # Find column indices of reverse complement for each kmer
  kmers <- tolower(colnames(onf.matrix))
  dgen.key <- stri_replace_all(kmers,
                               fixed = c("a", "c", "g", "t"),
                               replacement = c("T", "G", "C", "A"),
                               vectorize_all = FALSE)
  dgen.key <- stri_reverse(tolower(dgen.key))
  dgen.key <- match(dgen.key, kmers)
  cols2sum <- dgen.key > 1:length(kmers)
  cols2del <- dgen.key < 1:length(kmers)
  
  # Combine reverse complement kmers
  onf.dgen <- onf.matrix
  # sum kmers + reverse complements
  onf.dgen[, cols2sum] <- onf.dgen[, cols2sum] + onf.dgen[, dgen.key[cols2sum]]
  # remove reverse complements from matrix
  onf.dgen <- onf.dgen[, !cols2del]
  # do nothing to kmers that are their own reverse complement
  
  # Save degenerated ONF matrix as RDS file
  saveRDS(onf.dgen, sub("\\.\\w*$", "_dgen.rds", f))
}
