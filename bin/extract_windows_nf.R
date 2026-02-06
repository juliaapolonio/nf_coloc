#!/usr/bin/env Rscript

library(dplyr)
library(GenomicRanges)
library(vroom)

# Sumstats must have the following columns: chromosome, base_pair_location, p_value

# Add parameters: chr_col, bp_col, pval_col
# todo: parametrizar window_size, overlap_size e pval
args <- commandArgs(trailingOnly = TRUE)

# Inputs
sumstats_path <- args[1]

sumstats <- vroom(sumstats_path)
index_variants <- subset(sumstats, select = c("chr", "bp", "p"))
colnames(index_variants) = c("chr", "pos", "pval")
index_variants <- index_variants %>% filter(pval < 5e-8)
index_variants <- na.omit(index_variants)


# Defining function to extract windows
define_coloc_regions <- function(index_variants) {
  # Step 1: sort by p-value
  index_variants <- index_variants %>%
    arrange(pval)
  
  # Step 2: create 500kb windows
  gr <- GRanges(seqnames = index_variants$chr,
                ranges = IRanges(start = index_variants$pos - 250000,
                                 end = index_variants$pos + 250000),
                pval = index_variants$pval)
  
  # Step 3: remove less significant overlapping windows
  retained <- GRanges()
  for (i in seq_along(gr)) {
    overlaps <- findOverlaps(gr[i], retained)
    if (length(overlaps) == 0) {
      retained <- c(retained, gr[i])
    }
  }
  
  # Step 4: resolve region overlaps (>200kb = merge, ≤200kb = split)
  # Sort retained by start
  retained <- sort(retained)
  
  resolved <- GRanges()
  i <- 1
  while (i <= length(retained)) {
    current <- retained[i]
    j <- i + 1
    while (j <= length(retained)) {
      next_region <- retained[j]
      ov <- pintersect(current, next_region, drop.nohit.ranges = TRUE)
      if (length(ov) > 0) {
        ov_width <- width(ov)
        if (ov_width > 200000) {
          # Merge and move on
          merged <- reduce(c(current, next_region))
          current <- merged[1]
          j <- j + 1
        } else {
          # Split the overlapping part in half
          midpoint <- start(ov) + floor(width(ov) / 2)
          end(current) <- midpoint
          start(next_region) <- midpoint + 1
          retained[j] <- next_region
          break  # Only split with immediate neighbor
        }
      } else {
        break
      }
    }
    resolved <- c(resolved, current)
    i <- j
  }
  
  return(as.data.frame(resolved)[, c("seqnames", "start", "end")])
}


# Call function with input
coloc_regions <- define_coloc_regions(index_variants)

vroom_write(coloc_regions, "windows.tsv")