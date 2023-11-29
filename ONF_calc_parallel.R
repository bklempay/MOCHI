#!/usr/bin/env Rscript

library(ape)
library(stringi)
library(parallel)

#### Default options ####

prefix <- "mochi" # output file prefix
format <- "tsv" # output file format ["tsv"|"rds"]
k <- 4 # choose kmer length (default 4)
nthreads <- detectCores() # default all available cores
window_size <- 5000 # sliding window size (default 5 kb)
window_step <- 1000 # sliding window step size (default 1 kb)
# if you don't wish to use sliding windows, make the window and step sizes arbitrarily large
rm_plasmids <- FALSE # exclude plasmids due to uncharacteristic ONFs?
cat_contigs <- FALSE # concatenate contigs into a single, long sequence?
subsample <- FALSE # subsample sliding windows? [integer]
rand_seed <- 4444 # set random seed for reproducibility if subsampling


#### Command line arguments ####

input <- commandArgs(TRUE) # pass command line arguments to R
files <- character(0) # empty vector to append sequence files

# Parse command line input
for(arg in input){
  if(grepl("=", arg)){
    arg <- unlist(strsplit(arg, "=")) # split argument into variable name (arg[1]) and value (arg[2])
    if(!(arg[1] %in% ls())) stop(paste0("invalid command line argument (",arg[1],")")) # flag invalid arguments
    assign(arg[1], arg[2]) # if there are multiple equal signs, everything after the second will be disregarded
    if(grepl("^\\d+$", arg[2])) assign(arg[1], as.numeric(arg[2])) # convert to numeric if possible
  } else files <- c(files, arg)
}
stopifnot(file.exists(files),
          grepl("\\.(fasta|fas|fa|fna|ffn)", files, ignore.case = TRUE),
          format %in% c("tsv","rds"))


#### Calculate oligonucleotide frequencies ####

# Generate possible kmers
kmers <- expand.grid(rep(list(c("a","c","g","t")), k), stringsAsFactors = FALSE)
kmers <- apply(kmers, MARGIN = 1, FUN = paste, collapse = "")

# Create parallel socket cluster
clust <- makeCluster(nthreads)
clusterExport(clust, c("window_size","window_step","rm_plasmids","cat_contigs","subsample","k","kmers"))
if(rand_seed) clusterSetRNGStream(clust, rand_seed) # set random seed (recommended if subsampling)

cat("Generating sliding windows...\n")

# Run all sequence files in parallel
seqs <- parLapply(clust, files, function(f){
  # read sequence file (fasta format)
  fasta <- ape::read.dna(f, format = "fasta", as.character = TRUE, as.matrix = FALSE)
  # exclude contigs with sequence headers containing the word "plasmid" (optional)
  if(rm_plasmids) fasta <- fasta[!grepl("plasmid", names(fasta), ignore.case = TRUE)]
  fasta.cat <- sapply(fasta, FUN = paste, collapse = "")
  # concatenate contigs into a single, long string separated by k undetermined (n) bases (optional)
  if(cat_contigs) fasta.cat <- paste(fasta.cat, collapse = strrep("n",k))
  # rename contigs using the format path:`/full/local/path`.contig:`contig name`
  names(fasta.cat) <- paste0("path:`",f,"`.contig:`",names(fasta.cat),"`")
  
  # generate sliding windows
  lapply(fasta.cat, function(contig){
    # determine the number of sliding windows in each contig
    n.window <- max(1, round((nchar(contig)-window_size)/window_step) + 1)
    # if the final window is shorter by more than half of step size, discard it entirely
    windows <- sapply(1:n.window, FUN = function(i){
      substr(contig, (i-1)*window_step+1, (i-1)*window_step+window_size)
    })
    # rename windows using the format start:int.end:int
    names(windows) <- paste0("start:[",as.integer((1:n.window-1)*window_step+1),
                             "].end:[",as.integer((1:n.window-1)*window_step+window_size),"]")
    # discard sliding windows with > 10% undermined (n) bases
    windows <- windows[stringi::stri_count_fixed(windows, "n") < sapply(windows, nchar)*0.1]
    # in some cases, it may be necessary to subsample due to memory/time limitations for downstream steps
    if(subsample && subsample < length(windows)) windows <- windows[sort(sample(1:length(windows), subsample))]
    return(windows)
  })
})

cat("Calculating oligonucleotide frequencies...\n")

# Calculate oligonucleotide frequencies for each sliding window
onf.matrix <- t(parSapplyLB(clust, unlist(seqs), FUN = stringi::stri_count, fixed = kmers, overlap = TRUE))
onf.matrix <- onf.matrix/rowSums(onf.matrix) # normalize kmer counts
colnames(onf.matrix) <- kmers

stopCluster(clust)

# Save ONF matrix
if(format == "tsv") write.table(onf.matrix, paste0(prefix,"_ONF_matrix.tsv"), sep = "\t", quote = FALSE)
if(format == "rds") saveRDS(onf.matrix, paste0(prefix,"_ONF_matrix.rds"))


#### Calculate degenerated oligonucleotide frequencies ####

# Find column indices of reverse complement for each kmer
dgen.key <- stri_replace_all(kmers, fixed = c("a","c","g","t"), c("T","G","C","A"), vectorize_all = FALSE)
dgen.key <- stri_reverse(tolower(dgen.key))
dgen.key <- match(dgen.key, kmers)
cols2sum <- dgen.key > 1:length(kmers)
cols2del <- dgen.key < 1:length(kmers)

onf.dgen <- onf.matrix
onf.dgen[,cols2sum] <- onf.dgen[,cols2sum] + onf.dgen[,dgen.key[cols2sum]] # sum kmers + reverse complements
onf.dgen <- onf.dgen[,!cols2del] # remove reverse complements from matrix
# do nothing to kmers that are their own reverse complement

# Save degenerated ONF matrix
if(format == "tsv") write.table(onf.dgen, paste0(prefix,"_ONF_dgen.tsv"), sep = "\t", quote = FALSE)
if(format == "rds") saveRDS(onf.dgen, paste0(prefix,"_ONF_dgen.rds"))