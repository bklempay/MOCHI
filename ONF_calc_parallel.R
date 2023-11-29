#!/usr/bin/env Rscript

library(ape)
library(stringi)
library(parallel)

#### Default options ####

prefix <- "mochi"     # output file prefix
format <- "rds"       # output file format ["rds"|"tsv"]
k <- 4                # kmer length
window_size <- 5000   # sliding window dimensions (to bypass sliding windows,
window_step <- 1000   # make the window and step sizes arbitrarily large)
rm_plasmids <- FALSE  # exclude plasmids due to uncharacteristic ONFs?
cat_contigs <- FALSE  # concatenate contigs into a single, long sequence?
subsample <- FALSE    # subsample sliding windows? [integer]
rand_seed <- 4444     # set random seed for reproducibility if subsampling
nthreads <- detectCores() # default all available cores


#### Command line arguments ####

input <- commandArgs(TRUE) # pass command line arguments to R
files <- character(0) # empty vector to append sequence files

# Parse command line input
for (arg in input) {
  if (!grepl("=", arg)) {
    files <- c(files, arg)
  } else {
    # split argument into variable name (arg[1]) and value (arg[2])
    # if multiple equal signs, everything after the second will be disregarded
    arg <- unlist(strsplit(arg, "="))
    # flag invalid variable names
    if (!(arg[1] %in% ls())) {
      stop("invalid command line argument (", arg[1], ")")
    }
    # assign value to variable name
    assign(arg[1], arg[2])
    # convert values to numeric if possible
    if (grepl("^\\d+$", arg[2])) assign(arg[1], as.numeric(arg[2]))
  }
}

stopifnot(file.exists(files), 
          grepl("\\.(fasta|fas|fa|fna|ffn)", files, ignore.case = TRUE), 
          format %in% c("rds", "tsv"))


#### Calculate oligonucleotide frequencies ####

# Generate possible kmers
kmers <- expand.grid(rep(list(c("a", "c", "g", "t")), k),
                     stringsAsFactors = FALSE)
kmers <- apply(kmers, 1, paste, collapse = "")

# Create parallel socket cluster
clust <- makeCluster(nthreads)
clusterExport(clust, c("window_size", "window_step", "rm_plasmids",
                       "cat_contigs", "subsample", "k", "kmers"))

# Set random seed (recommended if subsampling)
if (rand_seed) clusterSetRNGStream(clust, rand_seed)

cat("Generating sliding windows...\n")

# Run all sequence files in parallel
seqs <- parLapply(clust, files, function(f) {
  # read sequence file (fasta format)
  fasta <- ape::read.dna(f,
                         format = "fasta",
                         as.character = TRUE,
                         as.matrix = FALSE)
  # exclude contigs with headers containing the word "plasmid" (optional)
  if (rm_plasmids) {
    fasta <- fasta[!grepl("plasmid", names(fasta), ignore.case = TRUE)]
  }
  # collapse DNA bases into long character strings
  fasta.cat <- sapply(fasta, paste, collapse = "")
  # concatenate contigs into a single, long string (optional)
  # separated by k undetermined (n) bases
  if (cat_contigs) fasta.cat <- paste(fasta.cat, collapse = strrep("n", k))
  # rename contigs using the format path:</full/local/path>.contig:<contig name>
  names(fasta.cat) <- paste0("path:<", f, ">.contig:<", names(fasta.cat), ">")
  
  # Generate sliding windows
  lapply(fasta.cat, function(contig) {
    # determine the number of sliding windows in each contig
    # if the final window is shorter by more than half of step size, discard it
    n.window <- max(1, round((nchar(contig) - window_size) / window_step) + 1)
    # determine the start and end positions for each window
    window.pos <- list(start = (1:n.window - 1) * window_step + 1,
                       end = (1:n.window - 1) * window_step + window_size)
    # trim each contig from start position to end position
    windows <- mapply(substr,
                      x = contig,
                      start = window.pos$start,
                      stop = window.pos$end)
    # rename windows using the format start:<int>.end:<int>
    names(windows) <- paste0("start:<", as.integer(window.pos$start),
                             ">.end:<", as.integer(window.pos$end), ">")
    # discard sliding windows with > 10% undermined (n) bases
    n.percent <- stri_count_fixed(windows, "n") / sapply(windows, nchar)
    windows <- windows[n.percent <= 0.1]
    # subsample due to memory/time limitations for downstream steps (optional)
    if (subsample && subsample < length(windows)) {
      windows <- windows[sort(sample(1:length(windows), subsample))]
    }
    return(windows)
  })
})

cat("Calculating oligonucleotide frequencies...\n")

# Calculate oligonucleotide frequencies for each sliding window
onf.matrix <- t(parSapplyLB(clust, unlist(seqs), stringi::stri_count,
                            fixed = kmers,
                            overlap = TRUE))
onf.matrix <- onf.matrix/rowSums(onf.matrix) # normalize kmer counts
colnames(onf.matrix) <- kmers

stopCluster(clust)

# Save ONF matrix
switch(format,
       rds = saveRDS(onf.matrix, paste0(prefix, "_ONF_matrix.rds")),
       tsv = write.table(onf.matrix, paste0(prefix, "_ONF_matrix.tsv"),
                         sep = "\t",
                         quote = FALSE))


#### Calculate degenerated oligonucleotide frequencies ####

# Find column indices of reverse complement for each kmer
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

# Save degenerated ONF matrix
switch(format,
       rds = saveRDS(onf.dgen, paste0(prefix, "_ONF_dgen.rds")),
       tsv = write.table(onf.dgen, paste0(prefix, "_ONF_dgen.tsv"),
                         sep = "\t",
                         quote = FALSE))