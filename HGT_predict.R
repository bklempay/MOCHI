#!/usr/bin/env Rscript

library(ape)
library(stringi)

#### Default options ####

outdir <- "."         # output file directory
prefix <- "mochi"     # output file prefix
mode <- "sensitive"   # conservative, sensitive, or ultrasensitive
reference <- NULL     # compare input to a reference distribution [ONF matrix]
hgt_pval <- 0.001     # p-value cutoff for predicting HGT
grp_bins <- TRUE      # genomes/MAGs -> TRUE; unbinned contigs -> FALSE
one_fasta <- TRUE     # output one single FASTA, or one per input (meta)genome?


#### Command line arguments ####

input <- commandArgs(TRUE) # pass command line arguments to R
file <- character(0) # empty vector to append ONF matrix file

# Parse command line input
for (arg in input) {
  if (!grepl("=", arg)) {
    file <- c(file, arg)
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
  file.exists(file),
  file.exists(as.character(reference)),
  grepl("\\.rds$", c(file, reference), ignore.case = TRUE),
  0 < hgt_pval && hgt_pval < 1,
  is.logical(c(grp_bins, one_fasta))
)

# Create output file directory (if it doesn't already exist)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
# append output file directory to prefix
prefix <- paste0(outdir, "/", prefix)


#### Get ONF matrices ####

# For now, HGT_predict.R can only handle one ONF matrix at a time
if (length(file) > 1) {
  cat("HGT_predict.R can only handle one ONF matrix at a time.\n",
      "Proceeding with ", file[1], "...\n", sep = "")
}
file <- file[1]

# Read ONF matrix (preferably degenerated)
onf.matrix <- readRDS(file)
# Get kmer/sliding window parameters
k <- attributes(onf.matrix)$k
window_size <- attributes(onf.matrix)$window_size
window_step <- attributes(onf.matrix)$window_step

# Extract sequence file paths and contig names from row names
paths <- stri_match(rownames(onf.matrix), regex = "path:'(.*)'_contig")[, 2]
contigs <- stri_match(rownames(onf.matrix), regex = "contig:'(.*)'_start")[, 2]

# Group sliding windows by genome/MAG or by contig depending on analysis mode:
# each group is expected to share a null (non-HGT) ONF distribution
onf.list <- split(data.frame(onf.matrix),
                  f = if (grp_bins) paths else list(paths, contigs),
                  drop = TRUE,
                  sep = ">")

# Use query ONF matrix as training data for the null ONF distribution
# else (if provided) get reference ONF matrix (assumed to be a single genome)
onf.data <- if (is.null(reference)) onf.list else list(readRDS(reference))

#### Model the ONF distribution of non-HGT sliding windows ####

# High-dimensional ONF mode estimation
onf.mode <- lapply(onf.data, function(onf) {
  # kmer counts can be modeled as a Poisson distribution, though the tails will
  # be slightly fatter due to non-independence of overlapping sliding windows.
  # Determine the log-likelihood of each sliding window as if the true null ONF
  # distribution were the composite of independent, Poisson-distributed kmers
  # (not the case, but it's a good way to rank distance from the ONF mode).
  kmer.loglik <- apply(onf, 2, function(kmer) {
    kmer <- round(kmer * (window_size - k + 1))
    # begin by truncating outliers (probable HGT events)
    kmer.sub <- kmer[kmer >= qpois(0.001, mean(kmer))] # use hgt_pval here?
    kmer.sub <- kmer.sub[kmer.sub <= qpois(0.999, mean(kmer))]
    # fit a new Poisson distribution to the truncated data
    kmer.lambda <- mean(kmer.sub) # no need to do anything fancy!
    # return the log-likelihood for each kmer
    return(dpois(kmer, kmer.lambda, log = TRUE))
  })
  onf.loglik <- rowSums(kmer.loglik)
  onf.rankdist <- rank(-onf.loglik) / nrow(onf)
  
  # Compute a weighted average of sliding window ONFs, where weights decrease
  # logistically with rank distance (inflection point @ 80th percentile)
  onf.weight <- 1 / (1 + exp(80 * (onf.rankdist - 0.8)))
  # experimentally, I think that 80 is a good scale coefficient but this might
  # need to be tuned in the future
  
  return(colSums(onf * onf.weight) / sum(onf.weight))
})

# PCA ordination in order to reduce correlation between dimensions in ONF space
pca.list <- mapply(prcomp, x = onf.data, center = onf.mode, SIMPLIFY = FALSE)
# Overlapping kmers, e.g. ATCA and TCAA, are strongly correlated with each other
# (assume linear correlations). Manhattan distance in ONF space therefore double
# counts the influence of such kmers. PCA effectively weights each dimension in
# order to reduce double counting, while preserving distances between points.

# Calculate scaled-PCA Manhattan distances of sliding windows to the ONF mode
# i.e. Manhattan distances from the origin in scaled PCA ordination space
pca.dist <- lapply(pca.list, function(pca) {
  pca.scale <- sweep(pca$x, 2, pca$sdev, "/")
  return(rowSums(abs(pca.scale)))
})

# Thanks to the Central Limit Theorem, the scaled-PCA Manhattan distances of
# non-HGT sliding windows should be normally distributed. HGT events will have
# a greater distance from the ONF mode, so the combined distributions of non-HGT
# and HGT sliding windows should be slightly right-skewed (depending on the rate
# of HGT). The most parsimonious way to parametrize the distribution of non-HGT
# sliding windows is to fit a normal distribution using only the left half of
# the combined distribution. This method assumes that all sliding windows in the
# left half of the combined distribution are non-HGT events.

# Estimate normal distribution parameters for non-HGT sliding windows
par.norm <- lapply(pca.dist, function(dist) {
  # 2 * median - mean is a simple and effective heuristic for estimating the
  # center (mode) of the non-HGT normal distribution, but this can (and should)
  # be verified visually by looking at mochi_figures.pdf
  dist.mode <- 2 * median(dist) - mean(dist)
  # estimate SD using Median Absolute Deviation to the left of the mode
  mad.left <- mad(dist[dist <= dist.mode], center = dist.mode)
  return(list(mean = dist.mode, sd = mad.left))
})

# Gaussian mixture modeling: this is an unimplemented idea for parametrizing the
# normal distribution of non-HGT sliding windows from a previous version
# gmm.list <- lapply(pca.dist, mclust::Mclust, G = 1:2, verbose = FALSE)
# par.norm <- lapply(gmm.list, function(gmm) {
#   return(list(mean = gmm$parameters$mean[1],
#               sd = sqrt(gmm$parameters$variance$sigmasq[1])))
# })


#### Figures ####

pdf(paste0(prefix, "_figures.pdf"), width = 6, height = 6)

for (i in 1:length(pca.dist)) {
  hist(pca.dist[[i]], breaks = 50, freq = FALSE,
       main = basename(names(pca.dist)[i]),
       xlab = "Sliding window distance from ONF mode")
  curve(dnorm(x, par.norm[[i]]$mean, par.norm[[i]]$sd), add = TRUE)
  abline(v = qnorm(1 - hgt_pval, par.norm[[i]]$mean, par.norm[[i]]$sd),
         lty = "dashed")
}

dev.off()


#### Predict candidate HGT events ####

hgt.pred <- function(pca, par, newdata = NULL) {
  # Convert new sliding windows to PCA ordination space (if required)
  pca.data <- if (!is.null(newdata)) predict(pca, newdata) else pca$x
  pca.dist <- rowSums(abs(sweep(pca.data, 2, pca$sdev, "/")))
  
  # Parse sliding window names
  windows <- data.frame(
    stri_match(str = rownames(pca.data),
               regex = "^path:'(.*)'_contig:'(.*)'_start:(\\d+)_end:(\\d+)$")
  )
  names(windows) <- c("name", "path", "contig", "start", "end")
  # convert start and end positions to integers
  windows <- type.convert(windows, as.is = TRUE)
  
  # Calculate p-value for each sliding window
  windows$pval <- pnorm(pca.dist, par$mean, par$sd, lower.tail = FALSE)
  
  # Break sliding windows into chunks the length of window_step
  windows <- split(windows, ~ contig)
  chunks <- lapply(windows, function(window) {
    starts <- seq(min(window$start), max(window$end), window_step)
    ends <- starts + window_step - 1
    ends[length(ends)] <- max(window$end)
    return(data.frame(path = unique(window$path),
                      contig = unique(window$contig),
                      start = starts,
                      end = ends))
  })
  
  # Predict candidate HGT events (p < hgt_pval)
  hgt.preds <- mapply(
    function(chunk, window) {
      # get indices for sliding windows which overlap with each chunk
      which.windows <- mapply(
        function(start, end) {
          return(which(window$start < end & window$end > start))
        },
        start = chunk$start,
        end = chunk$end
      )
      
      # Calculate p-value for each chunk
      chunk$pval <- sapply(which.windows, function(ind) {
        switch(mode,
               # ultraconservative mode: maximum of sliding window p-values
               ultraconservative = max(window[ind, "pval"]),
               # conservative mode: arithmetic mean of sliding window distances
               conservative = pnorm(mean(pca.dist[ind]), par$mean, par$sd,
                                    lower.tail = FALSE),
               # sensitive mode: geometric mean of sliding window p-values
               sensitive = exp(mean(log(window[ind, "pval"]))),
               # ultrasensitive mode: minimum of sliding window p-values
               ultrasensitive = min(window[ind, "pval"]))
      })
      return(chunk[chunk$pval < hgt_pval, ])
    },
    chunk = chunks,
    window = windows,
    SIMPLIFY = FALSE
  )
  
  # Merge contiguous chunks
  hgt.heads <- lapply(hgt.preds, function(preds) {
    return(which(c(preds$start, 0) > c(0, preds$end + 1)))
  })
  hgt.tails <- lapply(hgt.preds, function(preds) {
    return(which(c(Inf, preds$end) < c(preds$start - 1, Inf)) - 1)
  })
  hgt.merge <- mapply(
    function(preds, heads, tails) {
      return(data.frame(
        path = preds$path[heads],
        contig = preds$contig[heads],
        start = preds$start[heads],
        end = preds$end[tails]
      ))
    },
    preds = hgt.preds,
    heads = hgt.heads,
    tails = hgt.tails,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )
  
  # Return merged HGT predictions
  return(do.call(rbind, hgt.merge))
}

hgt.list <- 
  if (is.null(reference)) {
    mapply(hgt.pred, pca = pca.list, par = par.norm, SIMPLIFY = FALSE)
  } else {
    lapply(onf.list, hgt.pred, pca = pca.list[[1]], par = par.norm[[1]])
  }

# Discard empty elements of hgt.list (i.e. no HGT predicted)
hgt.list <- hgt.list[sapply(hgt.list, nrow) > 0]

if (length(hgt.list) == 0) {
  cat("Zero HGT candidates were detected.\n")
  quit(save = "no")
}

# Output HGT predictions matrix
write.table(do.call(rbind, hgt.list), paste0(prefix, "_HGT_pred.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


#### Extract candidate HGT sequences and output FASTA file(s) ####

# Extract candidate HGT sequences
hgt.seqs <- lapply(hgt.list, function(preds) {
  # read sequence file (fasta format)
  fasta <- ape::read.FASTA(unique(preds$path))
  # subset candidate HGT sequences
  fasta.sub <- mapply(
    function(contig, start, end) {
      return(fasta[[contig]][start:end])
    },
    contig = preds$contig,
    start = preds$start,
    end = preds$end
  )
  # rename sequence headers using the format:
  # path:'/full/local/path/'_contig:'contig name'_start:[int]_end:[int]
  names(fasta.sub) <- paste0("path:'", preds$path, "'_contig:'", preds$contig,
                             "'_start:", preds$start, "_end:", preds$end)
  return(fasta.sub)
})

# Output FASTA file(s)
if (one_fasta) {
  # either one single FASTA (default)
  hgt.seqs <- unlist(hgt.seqs, recursive = FALSE)
  class(hgt.seqs) <- "DNAbin"
  ape::write.FASTA(hgt.seqs, paste0(prefix, "_HGT_pred.fasta"))
} else {
  # or one per input (meta)genome
  invisible(sapply(unique(paths), function(path) {
    is.query <- grepl(path, names(hgt.seqs))
    hgt.seqs <- unlist(hgt.seqs[is.query], recursive = FALSE)
    class(hgt.seqs) <- "DNAbin"
    # truncate directory name and file extension from (meta)genome path
    name <- gsub("^.*/|\\.(fasta|fas|fa|fna|ffn)(\\.gz)?$", "", path)
    ape::write.FASTA(hgt.seqs, paste0(outdir, "/", name, "_HGT_pred.fasta"))
  }))
}