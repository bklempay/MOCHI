#!/usr/bin/env Rscript

library(ape)
library(stringi)
library(parallel)

#### Default options ####

outdir <- "."         # output file directory
prefix <- "mochi"     # output file prefix
format <- "rds"       # output file format ["rds"|"tsv"]
k <- 4                # kmer length
window_size <- 3000   # sliding window dimensions (to bypass sliding windows,
window_step <- 1000   # make the window and step sizes arbitrarily large)
rm_plasmids <- FALSE  # exclude plasmids due to uncharacteristic ONFs?
trim_bp <- FALSE      # trim head and tail of each contig? [integer]
cat_contigs <- FALSE  # concatenate contigs into a single, long sequence?
min_length <- TRUE    # minimum contig length? (default 0.9 * window_size)
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
      stop("invalid command line argument \"", arg[1], "\"")
    }
    # assign value to variable name
    assign(arg[1], arg[2])
  }
}

# Convert variables to numeric/integer/logical as appropriate
for (var in ls()) assign(var, type.convert(get(var), as.is = TRUE))

# Check for illegal user input
stopifnot(
  file.exists(files),
  grepl("\\.(fasta|fas|fa|fna|ffn)(\\.gz)?$", files, ignore.case = TRUE),
  format %in% c("rds", "tsv"),
  is.integer(c(k, window_size, window_step, nthreads)),
  is.logical(c(rm_plasmids, cat_contigs)),
  class(c(min_length, subsample, rand_seed)) %in% c("integer", "logical")
)

# Create output file directory (if it doesn't already exist)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
# append output file directory to prefix
prefix <- paste0(outdir, "/", prefix)

# By default, minimum contig length equals 0.9 * window_size
if(min_length == TRUE) min_length <- 0.9 * window_size


#### Calculate oligonucleotide frequencies ####

# Generate possible kmers
kmers <- expand.grid(rep(list(c("a", "t", "c", "g")), k))
kmers <- apply(kmers, 1, paste, collapse = "")

# Create parallel socket cluster
clust <- makeCluster(nthreads)
clusterExport(clust, c("k", "window_size", "window_step", "rm_plasmids",
                       "cat_contigs", "min_length", "subsample", "kmers"))

cat("Generating sliding windows...\n")

# Run all sequence files in parallel
seqs <- parLapplyLB(clust, files, function(f) {
  # read sequence file (fasta format)
  fasta <- ape::read.dna(f,
                         format = "fasta",
                         as.character = TRUE,
                         as.matrix = FALSE)
  # exclude contigs with headers containing the word "plasmid" (optional)
  if (rm_plasmids) {
    fasta <- fasta[!grepl("plasmid", names(fasta), ignore.case = TRUE)]
  }
  # trim head and tail of each contig
  if (trim_bp) {
    fasta <- sapply(fasta, function(contig) {
      contig <- contig[-(1:trim_bp)]
      contig <- contig[-((length(contig) - trim_bp + 1):length(contig))]
      return(na.omit(contig))
    })
  }
  # collapse DNA bases into long character strings
  fasta.cat <- sapply(fasta, paste, collapse = "")
  # concatenate contigs into a single, long string (optional)
  # separated by k undetermined (n) bases
  if (cat_contigs) fasta.cat <- paste(fasta.cat, collapse = strrep("n", k))
  # discard contigs shorter than minimum length (optional)
  if (min_length) fasta.cat <- fasta.cat[nchar(fasta.cat) >= min_length]
  # rename contigs using the format path:'/full/local/path'_contig:'contig name'
  names(fasta.cat) <- paste0("path:'", f, "'_contig:'", names(fasta.cat), "'")
  
  # Generate sliding windows
  lapply(fasta.cat, function(contig) {
    # determine the number of sliding windows in each contig
    # if the final window is shorter by more than half of step size, discard it
    n.window <- max(1, round((nchar(contig) - window_size) / window_step) + 1)
    # determine the start and end positions for each window
    window.pos <- list(start = (1:n.window - 1) * window_step + 1,
                       end = (1:n.window - 1) * window_step + window_size)
    window.pos$end[n.window] <- min(window.pos$end[n.window], nchar(contig))
    # trim each contig from start position to end position
    windows <- mapply(substr,
                      x = contig,
                      start = window.pos$start,
                      stop = window.pos$end)
    # rename windows using the format start:[int]_end:[int]
    names(windows) <- paste0("start:", as.integer(window.pos$start),
                             "_end:", as.integer(window.pos$end))
    # discard sliding windows with > 10% undermined (n) bases
    n.percent <- stringi::stri_count_fixed(windows, "n") / nchar(windows)
    windows <- windows[n.percent <= 0.1]

    return(windows)
  })
})

# Subsample due to memory/time limitations for downstream steps (optional)
if (subsample) {
  # set random seed (recommended if subsampling)
  if (rand_seed) clusterSetRNGStream(clust, rand_seed)
  seqs <- parLapply(clust, seqs, function(contigs) {
    lapply(contigs, function(windows) {
      if (subsample < length(windows)) {
        windows <- windows[sort(sample(1:length(windows), subsample))]
      }
      return(windows)
    })
  })
}

cat("Calculating oligonucleotide frequencies...\n")

# Calculate oligonucleotide frequencies for each sliding window
onf.calc <- t(parSapplyLB(clust, unlist(seqs), stringi::stri_count,
                          fixed = kmers,
                          overlap = TRUE))
onf.calc <- onf.calc / rowSums(onf.calc) # normalize kmer counts
colnames(onf.calc) <- kmers

stopCluster(clust)

# Update sliding window names to the format:
# path:'/full/local/path'_contig:'contig name'_start:[int]_end:[int]
rownames(onf.calc) <- gsub("'.start", "'_start", rownames(onf.calc))

# Save ONF matrix
switch(format,
       rds = saveRDS(onf.calc, paste0(prefix, "_ONF_calc.rds")),
       tsv = write.table(onf.calc, paste0(prefix, "_ONF_calc.tsv"),
                         sep = "\t",
                         quote = FALSE))


#### Calculate degenerated oligonucleotide frequencies ####

# Find column indices of reverse complement for each kmer
dgen.key <- stri_replace_all(kmers,
                             fixed = c("a", "t", "c", "g"),
                             replacement = c("T", "A", "G", "C"),
                             vectorize_all = FALSE)
dgen.key <- stri_reverse(tolower(dgen.key))
dgen.key <- match(dgen.key, kmers)
cols2sum <- dgen.key > 1:length(kmers)
cols2del <- dgen.key < 1:length(kmers)

# Combine reverse complement kmers
onf.dgen <- onf.calc
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