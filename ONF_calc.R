#!/usr/bin/env Rscript

library(ape)
library(stringi)

#### Default options ####

outdir <- "."         # output file directory
prefix <- "mochi"     # output file prefix
format <- "rds"       # output file format ["rds"|"tsv"]
k <- 4                # kmer length
window_size <- 5000   # sliding window dimensions (to bypass sliding windows,
window_step <- 1000   # make the window and step sizes arbitrarily large)
rm_plasmids <- FALSE  # exclude plasmids due to uncharacteristic ONFs?
cat_contigs <- FALSE  # concatenate contigs into a single, long sequence?
min_length <- TRUE    # minimum contig length? (default 0.9 * window_size)
subsample <- FALSE    # subsample sliding windows? [integer]
rand_seed <- 4444     # set random seed for reproducibility if subsampling


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
  }
}

# Convert variables to numeric/integer/logical as appropriate
for (var in ls()) assign(var, type.convert(get(var), as.is = TRUE))

# Check for illegal user input
stopifnot(
  file.exists(files),
  grepl("\\.(fasta|fas|fa|fna|ffn)(\\.gz)?$", files, ignore.case = TRUE),
  format %in% c("rds", "tsv"),
  is.integer(c(k, window_size, window_step)),
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

# Create output file
if (file.exists(paste0(prefix, "_ONF_calc.tsv"))) {
  stop(prefix, "_ONF_calc.tsv already exists")
}
write.table(matrix(c("", kmers), nrow = 1), paste0(prefix, "_ONF_calc.tsv"),
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Set random seed (recommended if subsampling)
if (rand_seed) set.seed(rand_seed)

cat("Calculating oligonucleotide frequencies...\n")

# Run for loop on all sequence files...
# this could and should be done in parallel for large datasets
for (i in 1:length(files)) {
  f <- files[i]
  cat("(", i, "/", length(files), ") Working on ", f, "\n", sep = "")
  # read sequence file (fasta format)
  fasta <- read.dna(f, format = "fasta", as.character = TRUE, as.matrix = FALSE)
  # exclude contigs with headers containing the word "plasmid" (optional)
  if (rm_plasmids) {
    fasta <- fasta[!grepl("plasmid", names(fasta), ignore.case = TRUE)]
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
  seqs <- lapply(fasta.cat, function(contig) {
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
    # rename windows using the format start:[int]_end:[int]
    names(windows) <- paste0("start:", as.integer(window.pos$start),
                             "_end:", as.integer(window.pos$end))
    # discard sliding windows with > 10% undermined (n) bases
    n.percent <- stri_count_fixed(windows, "n") / sapply(windows, nchar)
    windows <- windows[n.percent <= 0.1]
    # subsample due to memory/time limitations for downstream steps (optional)
    if (subsample && subsample < length(windows)) {
      windows <- windows[sort(sample(1:length(windows), subsample))]
    }
    return(windows)
  })
  
  # Calculate oligonucleotide frequencies for each sliding window
  onf <- t(sapply(unlist(seqs), stri_count, fixed = kmers, overlap = TRUE))
  onf <- onf / rowSums(onf) # normalize kmer counts
  
  # Update sliding window names to the format:
  # path:'/full/local/path'_contig:'contig name'_start:[int]_end:[int]
  rownames(onf) <- gsub("'.start", "'_start", rownames(onf))
  
  # Output kmer counts to 'ONF_calc.tsv'
  write.table(onf, paste0(prefix, "_ONF_calc.tsv"),
              append = TRUE,
              sep = "\t", 
              quote = FALSE,
              row.names = TRUE,
              col.names = FALSE)
}


#### Calculate degenerated oligonucleotide frequencies ####

# Read full ONF matrix
onf.calc <- data.matrix(
  read.delim(paste0(prefix, "_ONF_calc.tsv"), row.names = 1)
)
if (format == "rds") saveRDS(onf.calc, paste0(prefix, "_ONF_calc.rds"))

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