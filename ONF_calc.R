options(repos = list(CRAN = "http://cran.rstudio.com/"))
if(!require(ape, quietly = TRUE)) install.packages("ape"); library(ape)
if(!require(stringi, quietly = TRUE)) install.packages("stringi"); library(stringi)


#### User parameters ####

input <- commandArgs(TRUE) # pass bash input to R
prefix <- "refseq" # output file prefix
seq_meta <- "refseq_genomes_meta.tsv" # sequence metadata

window.size <- 5000 # sliding window size (default 5 kb)
window.step <- 1000 # sliding window step size (default 1 kb)
subsample <- FALSE # default FALSE

k <- 4 # choose kmer length (default 4)
kmers <- expand.grid(rep(list(c("a","c","g","t")), k), stringsAsFactors = FALSE)
kmers <- apply(kmers, MARGIN = 1, FUN = paste, collapse = "")

# Careful not to accidentally overwrite anything!
if(paste0(prefix,"_ONF_matrix.rds") %in% list.files()) stop(paste0(prefix,"_ONF_matrix.tsv already exists"))
if(paste0(prefix,"_ONF_dgen.rds") %in% list.files()) stop(paste0(prefix,"_ONF_dgen.rds already exists"))


#### Get sequence metadata ####

# Read tab-separated sequence metadata file (assume headers in line 1). You can create a metadata file
# in the correct format using `get_meta.R`. At a minimum, must contain a column named `local_path` with
# the full path to each sequence file. Sequence files must be in gzipped fasta format!
meta <- read.delim(seq_meta, quote = "", comment.char = "")
files <- meta$local_path


#### Calculate oligonucleotide frequencies ####

# Create output file
write.table(matrix(c("",kmers), nrow = 1), paste0(prefix,"_ONF_matrix.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Run for loop on all sequence files... this could and should be done in parallel for large data sets
for(i in 1:length(files)){
  f <- files[i]
  print(paste0("(",i,"/",length(files),") Working on ",f))
  fasta <- read.dna(gzfile(f), format = "fasta", as.character = TRUE, as.matrix = FALSE)
  fasta.sub <- fasta[!(grepl("plasmid", names(fasta)))] # exclude plasmids due to uncharacteristic ONFs
  fasta.cat <- paste(sapply(fasta.sub, FUN = paste, collapse = ""), collapse = strrep("n",k))
  # concatenate contig sequences into a single, long string separated by k undetermined (n) bases
  
  # determine the number of sliding windows in concatenated sequence
  n.window <- max(1, round((nchar(fasta.cat)-window.size)/window.step) + 1)
  # if the final window is shorter by more than half of step size, discard it entirely
  
  # generate sliding windows
  seqs <- sapply(1:n.window, FUN = function(i){
    substr(fasta.cat, (i-1)*window.step+1, (i-1)*window.step+window.size)
  })
  # discard sliding windows with > 10% undermined (n) bases
  seqs <- seqs[stri_count_fixed(seqs, "n") < window.size*0.1]
  if(length(seqs) == 0) next
  # in some cases, it will be necessary to subsample due to memory/time limitations for downstream steps
  if(is.numeric(subsample) && subsample < length(seqs)) seqs <- sample(seqs, subsample)
  
  # calculate oligonucleotide frequencies for each sliding window
  onf.matrix <- t(sapply(seqs, FUN = function(seq){
    onf <- sapply(kmers, FUN = stri_count_fixed, str = seq, overlap = TRUE) # count kmers
    onf/sum(onf) # normalize kmer counts
  }, USE.NAMES = FALSE))

  # output kmer counts to `ONF_matrix.tsv`
  write.table(onf.matrix, file = paste0(prefix,"_ONF_matrix.tsv"), append = TRUE, sep = "\t",
              quote = FALSE, row.names = paste0(f,"_",1:length(seqs)), col.names = FALSE)
}

# Read full ONF matrix and save as RDS file
onf.full <- data.matrix(read.delim(paste0(prefix,"_ONF_matrix.tsv"), row.names = 1))
saveRDS(onf.full, paste0(prefix,"_ONF_matrix.rds"))


#### Calculate degenerated oligonucleotide frequencies ####

# Find column indices of reverse complement for each kmer
dgen.key <- stri_replace_all(kmers, fixed = c("a","c","g","t"), c("T","G","C","A"), vectorize_all = FALSE)
dgen.key <- stri_reverse(tolower(dgen.key))
dgen.key <- match(dgen.key, kmers)
cols2sum <- dgen.key > 1:length(kmers)
cols2del <- dgen.key < 1:length(kmers)

onf.dgen <- onf.full
onf.dgen[,cols2sum] <-  onf.dgen[,cols2sum] + onf.dgen[,dgen.key[cols2sum]] # sum kmers + reverse complements
onf.dgen <- onf.dgen[,!cols2del] # remove reverse complements from matrix
# do nothing to kmers that are their own reverse complement

# Save degenerated ONF matrix as RDS file
saveRDS(onf.dgen, paste0(prefix,"_ONF_dgen.rds"))