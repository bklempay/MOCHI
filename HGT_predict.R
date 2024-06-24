#!/usr/bin/env Rscript

library(ape)
library(stringi)
library(robustbase)
library(vegan)

#### Default options ####

outdir <- "."         # output file directory
prefix <- "mochi"     # output file prefix
mode <- "sensitive"   # conservative, sensitive, or ultrasensitive
reference <- NULL     # compare input to a reference genome [ONF matrix]
nmds_dims <- 4        # number of NMDS dimensions
mcd_alpha <- 0.5      # MCD subset size parameter (between 0.5 and 1)
hgt_pval <- 0.001     # p-value cutoff for predicting HGT
grp_bins <- TRUE      # genomes/MAGs -> TRUE; unbinned contigs -> FALSE
one_fasta <- TRUE     # output one single FASTA, or one per input (meta)genome?
rand_seed <- 4444     # set random seed for reproducibility of NMDS and MCD
nthreads <- if (require(parallel, quietly = TRUE)) detectCores() else 1


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

# Extract sequence file paths and contig names from row names
paths <- stri_match(rownames(onf.matrix), regex = "path:'(.*)'_contig")[, 2]
contigs <- stri_match(rownames(onf.matrix), regex = "contig:'(.*)'_start")[, 2]

# Group sliding windows by genome/MAG or by contig depending on analysis mode:
# each group is expected to share a null (non-HGT) ONF distribution
group.by <-
  if (grp_bins) {
    as.factor(paths)
  } else {
    interaction(paths, contigs, drop = TRUE, sep = ">")
  }
group.by <- factor(group.by, levels = unique(group.by))
onf.list <- split(data.frame(onf.matrix), group.by)

# Read reference ONF matrix if provided (assumed to be a single genome)
if (!is.null(reference)) onf.ref <- readRDS(reference)

# Use query ONF matrix as training data for the null ONF distribution
# else (if provided) append query to reference ONF matrix
onf.data <- 
  if (is.null(reference)) {
    onf.list
  } else {
    stopifnot(!any(rownames(onf.ref) %in% rownames(query)))
    # you really shouldn't use a genome as a reference for itself, but if you
    # absolutely had to, you could just alter the reference ONF matrix row names
    lapply(onf.list, function(query) rbind(onf.ref, query))
  }


#### Dimension reduction with NMDS ####

cat("Reducing ONF dimensions... (this might take a while)\n")

if (rand_seed) set.seed(rand_seed)

# Non-metric multidimensional scaling
nmds.list <- lapply(onf.data, metaMDS,
                    distance = "euclidean",
                    k = nmds_dims,
                    autotransform = FALSE,
                    trace = FALSE,
                    parallel = nthreads)

# Use full NMDS coordinates to estimate the null ONF distribution
# else (if provided) use only the reference data points
nmds.null <- 
  if (is.null(reference)) {
    lapply(nmds.list, function(nmds) return(nmds$points))
  } else {
    lapply(nmds.list, function(nmds) return(nmds$points[1:nrow(onf.ref), ]))
  }


#### Parameterize null ONF distribution in NMDS space ####

cat("Parameterizing ONF distribution...\n")

# Estimate center and spread using Minimum Covariance Determinant
mcd.list <- lapply(nmds.null, function(nmds) {
  covMcd(nmds, alpha = mcd_alpha)
})

# Calculate squared Mahalanobis distances of query data points
mah.list <- 
  if (is.null(reference)) {
    lapply(mcd.list, function(mcd) return(mcd$mah))
  } else {
    mapply(
      function(nmds, mcd) {
        query <- nmds$points[-(1:nrow(onf.ref)), ]
        return(mahalanobis(query, mcd$center, mcd$cov))
      },
      nmds = nmds.list,
      mcd = mcd.list,
      SIMPLIFY = FALSE
    )
  }


#### Predict candidate HGT events ####

cat("Predicting candidate HGT events...\n")

hgt.list <- lapply(mah.list, function(mah) {
  # Parse sliding window names
  windows <- data.frame(
    stri_match(str = names(mah),
               regex = "^path:'(.*)'_contig:'(.*)'_start:(\\d+)_end:(\\d+)$")
  )
  names(windows) <- c("name", "path", "contig", "start", "end")
  # convert contigs to factor and start/end positions to integers
  windows <- type.convert(windows, as.is = FALSE)
  windows$contig <- factor(windows$contig, levels = unique(windows$contig))
  # there should only be one level in windows$path
  # if this is not the case, something went very wrong somewhere along the way
  
  windows$mah <- mah
  # Calculate p-value for each sliding window
  windows$pval <- pchisq(mah, nmds_dims, lower.tail = FALSE)
  
  windows <- split(windows, ~ contig)
  # Break sliding windows into chunks the length of window_step
  chunks <- lapply(windows, function(window) {
    window_step <- unique(window$start[-1] - window$start[-nrow(window)])
    # same as above, there should really only be one unique step size
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
          return(which(window$start <= start & window$end >= end))
          # consider including partially overlapping windows in the future
          # for now, I recommend making window_size a multiple of window_step
        },
        start = chunk$start,
        end = chunk$end
      )
      
      # Calculate p-value for each chunk
      chunk$pval <- sapply(which.windows, function(ind) {
        switch(mode,
               # ultraconservative mode: maximum of sliding window p-values
               ultraconservative = max(window$pval[ind]),
               # conservative mode: arithmetic mean of mahalanobis distances
               conservative = pchisq(mean(sqrt(window$mah[ind]))^2, nmds_dims,
                                     lower.tail = FALSE),
               # sensitive mode: geometric mean of sliding window p-values
               sensitive = exp(mean(log(window$pval[ind]))),
               # ultrasensitive mode: minimum of sliding window p-values
               ultrasensitive = min(window$pval[ind]))
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
})

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
  fasta <- ape::read.FASTA(levels(preds$path))
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


#### Figures ####

# Define mochi-themed color palette
palette.pval <- function(x) {
  ind <- round(x * 99 + 1)
  ramp <- colorRampPalette(c("#F8F8E0", "#672422")) # 小豆色 red bean color
  return(ramp(100)[ind])
}

pdf(paste0(prefix, "_figures.pdf"), width = 6, height = 6)

for (i in 1:length(nmds.list)) {
  # Plot the first 2 NDMS dimensions
  plot(
    if (is.null(reference)) {
      nmds.list[[i]]$points
    } else {
      nmds.list[[i]]$points[-(1:nrow(onf.ref)), ] # only the query data points
    },
    col = "#000000A0",
    bg = palette.pval(pchisq(mah.list[[i]], nmds_dims, lower.tail = FALSE)),
    pch = 21,
    main = gsub("^.*/|\\.(fasta|fas|fa|fna|ffn)(\\.gz)?$", "",
                basename(names(nmds.list)[i]))
  )
  mtext(paste0("NMDS: k = ",  nmds_dims,", stress = ",
               round(nmds.list[[i]]$stress, 3)))
  
  # Add contours for Mahalanobis distances at which p = 0.05, 0.01, and 0.001
  grid <- expand.grid(seq(par("usr")[1], par("usr")[2], length.out = 100),
                      seq(par("usr")[3], par("usr")[4], length.out = 100))
  for (j in which(1:nmds_dims > 2)) grid[, j] <- mcd.list[[i]]$center[j]
  grid$mah <- mahalanobis(grid, mcd.list[[i]]$center, mcd.list[[i]]$cov)
  grid$pval <- pchisq(grid$mah, df = nmds_dims, lower.tail = FALSE)
  
  contour(unique(grid$Var1), unique(grid$Var2),
          z = matrix(grid$pval, 100, 100),
          levels = c(0.05, 0.01, 0.001),
          add = TRUE,
          lty = "dashed")
}

dev.off()