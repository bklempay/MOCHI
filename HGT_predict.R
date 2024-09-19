#!/usr/bin/env Rscript

library(ape)
library(stringi)
library(robustbase)
library(vegan)

#### Default options ####

outdir <- "."           # output file directory
prefix <- "mochi"       # output file prefix
method <- "nmds"        # ordination method ["nmds"|"pcoa"|"pca"]
ndims <- 9              # number of ordination dimensions
distance <- "euclidean" # dissimilarity index (NMDS and PCoA only)
mcd_alpha <- 0.5        # MCD subset size parameter (between 0.5 and 1)
hgt_pval <- 0.001       # p-value cutoff for predicting HGT
grp_bins <- TRUE        # genomes/MAGs -> TRUE; unbinned contigs -> FALSE
one_fasta <- TRUE       # output one single FASTA, or one per input genome/MAG
rand_seed <- 4444       # set random seed for reproducibility of NMDS and MCD
reference <- NULL


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
group.by <- if (grp_bins) paths else paste0(paths, ">" , contigs)
group.by <- factor(group.by, levels = unique(group.by))
onf.list <- split(data.frame(onf.matrix), group.by)


#### Reduce ONF dimensionality ####

cat("Reducing ONF dimensions... (this might take a while)\n")

if (rand_seed) set.seed(rand_seed)

# Perform chosen ordination method
onf.reduce <- lapply(
  # use reference ONF matrix if provided
  X = if (is.null(reference)) onf.list else list(readRDS(reference)),
  FUN = function(onf) switch(
    tolower(method),
    nmds = metaMDS(onf,
                   distance = distance,
                   k = ndims,
                   autotransform = FALSE,
                   trace = FALSE),
    pcoa = cmdscale(vegdist(onf, method = distance),
                    k = ndims,
                    list. = TRUE),
    pca = prcomp(onf, rank. = ndims)
  )
)

# Extract sliding window coordinates in lower-dimensional ordination space
onf.coord <-
  if (method == "pca") {
    lapply(onf.reduce, function(pca) return(pca$x))
  } else {
    lapply(onf.reduce, function(mds) return(mds$points))
  }


#### Parameterize null ONF distribution in ordination space ####

cat("Parameterizing ONF distribution...\n")

# Estimate center and spread using Minimum Covariance Determinant
mcd.list <- lapply(onf.coord, function(coord) {
  covMcd(coord, alpha = mcd_alpha)
})

# Calculate squared Mahalanobis distances
mah.list <-
  if (is.null(reference)) {
    lapply(mcd.list, function(mcd) return(mcd$mah))
  } else {
    # if a reference was used for the ordination, first project query sliding
    # windows into the reference ordination space...
    lapply(onf.list, function(query) {
      query.pca <- predict(onf.reduce[[1]], query)
      # then calculate Mahalanobis distances to the reference distribution
      return(mahalanobis(query.pca, mcd.list[[1]]$center, mcd.list[[1]]$cov))
    })
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
  
  # Calculate p-value for each sliding window
  windows$mah <- mah
  windows$pval <- pchisq(mah, ndims, lower.tail = FALSE)
  
  # Break sliding windows into chunks the length of window_step
  windows <- split(windows, ~ contig)
  chunks <- lapply(windows, function(window) {
    window_step <-
      if (nrow(window) > 1) {
        unique(window$start[-1] - window$start[-nrow(window)])
        # same as above, there should really only be one unique step size
      } else window$end
    starts <- seq(min(window$start), max(window$end), window_step)
    ends <- starts + window_step - 1
    ends[length(ends)] <- max(window$end)
    return(data.frame(path = unique(window$path),
                      contig = unique(window$contig),
                      start = starts,
                      end = ends))
  })
  
  # Predict candidate HGT events
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
        end = chunk$end,
        SIMPLIFY = FALSE
      )
      
      # Give each chunk the p-value of its central sliding window
      chunk$pval <- sapply(which.windows, function(ind) {
        return(window$pval[round(median(ind))])
      })
      # return the chunks for which p < hgt_pval
      return(chunk[chunk$pval < hgt_pval, ])
    },
    chunk = chunks,
    window = windows,
    SIMPLIFY = FALSE
  )
  
  # Discard empty elements of hgt.preds (i.e. no HGT predicted)
  hgt.preds <- hgt.preds[sapply(hgt.preds, nrow) > 0]
  
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
        path = unique(preds$path),
        contig = unique(preds$contig),
        start = preds$start[heads],
        end = preds$end[tails],
        pval = exp(mean(log(preds$pval))) # geometric mean
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

# Discard NULL elements of hgt.list (i.e. no HGT predicted)
hgt.list <- hgt.list[!sapply(hgt.list, is.null)]

if (length(hgt.list) == 0) {
  cat("Zero HGT candidates were detected.\n")
  quit(save = "no")
}

# Output HGT predictions table
write.table(do.call(rbind, hgt.list), paste0(prefix, "_HGT_preds.tsv"),
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
      return(fasta[[as.character(contig)]][start:end])
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
names(hgt.seqs) <- NULL # keep only sub-list names
hgt.seqs <- unlist(hgt.seqs, recursive = FALSE)
class(hgt.seqs) <- "DNAbin"

if (one_fasta) {
  # either one single FASTA (default)
  ape::write.FASTA(hgt.seqs, paste0(prefix, "_HGT_preds.fasta"))
} else {
  # or one per input (meta)genome
  invisible(sapply(unique(paths), function(path) {
    is.query <- grepl(path, names(hgt.seqs))
    # truncate directory name and file extension from (meta)genome path
    name <- gsub("^.*/|\\.(fasta|fas|fa|fna|ffn)(\\.gz)?$", "", path)
    ape::write.FASTA(hgt.seqs[is.query],
                     paste0(outdir, "/", name, "_HGT_preds.fasta"))
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

for (i in 1:length(onf.coord)) {
  # Plot the first 2 ordination dimensions
  plot(
    onf.coord[[i]][, 1:ndims],
    col = "#000000A0",
    bg = palette.pval(pchisq(mah.list[[i]], ndims, lower.tail = FALSE)),
    pch = 21,
    main = gsub("^.*/|\\.(fasta|fas|fa|fna|ffn)(\\.gz)?$", "",
                names(onf.coord)[i])
  )
  mtext(switch(tolower(method),
               nmds = bquote(NMDS: ~ italic(N) ~ "=" ~ .(ndims) * ", stress =" ~
                               .(round(onf.reduce[[i]]$stress, 3))),
               pcoa = bquote(PCoA: ~ italic(N) ~ "=" ~ .(ndims)),
               pca = bquote(PCA: ~ italic(N) ~ "=" ~ .(ndims))
  ))
  
  # Add contours for Mahalanobis distances at which p = 0.05, 0.01, and 0.001
  grid <- expand.grid(seq(par("usr")[1], par("usr")[2], length.out = 100),
                      seq(par("usr")[3], par("usr")[4], length.out = 100))
  for (j in which(1:ndims > 2)) grid[, j] <- mcd.list[[i]]$center[j]
  grid$mah <- mahalanobis(grid, mcd.list[[i]]$center, mcd.list[[i]]$cov)
  grid$pval <- pchisq(grid$mah, df = ndims, lower.tail = FALSE)
  
  contour(unique(grid$Var1), unique(grid$Var2),
          z = matrix(grid$pval, 100, 100),
          levels = c(0.05, 0.01, 0.001),
          add = TRUE,
          lty = "dashed")
}

dev.off()
