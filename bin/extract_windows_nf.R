#!/usr/bin/env Rscript

library(dplyr)
library(GenomicRanges)
library(vroom)

args <- commandArgs(trailingOnly = TRUE)
clumped_path <- args[1]

clumped_data <- read.table(clumped_path, header = T) 

# Selecionando e renomeando as colunas necessárias
index_variants <- clumped_data %>%
  select(chr = CHR, pos = BP, pval = P) %>%
  filter(!is.na(chr) & !is.na(pos))

# Mantivemos a função idêntica à sua lógica original de mitigação de overlaps
define_coloc_regions <- function(index_variants) {
  # Step 1: sort by p-value
  index_variants <- index_variants %>%
    arrange(pval)
  
  # Step 2: create 500kb windows (±250kb ao redor da variante independente)
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
          # Merge e avança
          merged <- reduce(c(current, next_region))
          current <- merged[1]
          j <- j + 1
        } else {
          # Split o overlap no meio
          midpoint <- start(ov) + floor(width(ov) / 2)
          end(current) <- midpoint
          start(next_region) <- midpoint + 1
          retained[j] <- next_region
          break  # Splita apenas com o vizinho imediato
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

# Executa a função
coloc_regions <- define_coloc_regions(index_variants)

# Renomeia para manter o padrão de saída limpo
colnames(coloc_regions) <- c("chr", "start", "end")

# Salva o resultado
vroom_write(coloc_regions, "windows.tsv", delim = "\t")
