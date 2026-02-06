#!/usr/bin/env Rscript

library(coloc)
library(dplyr)
library(vroom)
library(stringr)
library(tidyr)
library(locuszoomr)
library(cowplot)
library(EnsDb.Hsapiens.v75)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

#qtl_path <- "/home/jamorim/scripts/nf_coloc/work/e0/29ccb1ea1422dee7a976bed00f416f/toy_eqtl_chr17.tsv"
#gwas_path <- "/home/jamorim/scripts/nf_coloc/work/e0/29ccb1ea1422dee7a976bed00f416f/toy_gwas_chr17.tsv"
#regions_path <- "/home/jamorim/scripts/nf_coloc/work/e0/29ccb1ea1422dee7a976bed00f416f/windows.tsv"

qtl_path <- args[1]
gwas_path <- args[2]
regions_path <- args[3]

qtl <- vroom(qtl_path, col_select = c("snp", "chr", "bp", "freq", "b", "se", "p", "N", "symbol"))
gwas <- vroom(gwas_path, col_select = c("snp", "chr", "bp", "freq", "b", "se", "p", "N"))
regions <- vroom(regions_path)

###### STEP 1 - Function to extract regions, split by gene and run coloc #######

# Adotar convenção de header e mudar no código
run_coloc_per_region <- function(regions, gwas_df, qtl_df) {
  results <- list()
  
  for (i in seq_len(nrow(regions))) {
    region_chr <- as.character(regions$seqnames[i])
    region_start <- regions$start[i]
    region_end <- regions$end[i]
    
    message("Processing region: ", region_chr, ":", region_start, "-", region_end)
    
    # Subset GWAS and QTL to region
    gwas_sub <- gwas_df %>%
      dplyr::filter(chr == region_chr & bp >= region_start & bp <= region_end)
    
    qtl_sub <- qtl_df %>%
      dplyr::filter(chr == region_chr & bp >= region_start & bp <= region_end)
    
    # Inner join by position to match variants
    merged <- inner_join(gwas_sub, qtl_sub, by = "snp", suffix = c(".gwas", ".qtl"))
    
    # Skip region if no overlap
    if (nrow(merged) == 0) next
    
    split_by_gene <- split(merged, merged$symbol)
    
    for (gene in names(split_by_gene)) {

      df <- split_by_gene[[gene]] %>%
        na.omit() %>%
        distinct(snp, .keep_all = TRUE)
      
      # Skip if too few SNPs
      if (nrow(df) < 5) {
        message("Skipped gene ", gene, ": too few overlapping variants")
        next
      }
      
      # Format for coloc.abf
      dataset1 <- list(
        beta = df$b.gwas,
        varbeta = df$se.gwas^2,
        snp = df$snp,
        MAF = df$freq.gwas,
        N = unique(df$N.gwas),
        type = "quant"
      )
      
      dataset2 <- list(
        beta = df$b.qtl,
        varbeta = df$se.qtl^2,
        snp = df$snp,
        MAF = df$freq.qtl,
        N = unique(df$N.qtl),
        type = "quant"
      )
      
      coloc_out <- tryCatch({
        coloc.abf(dataset1, dataset2)
      }, error = function(e) {
        message("coloc failed for gene ", gene, ": ", e$message)
        return(NULL)
      })

      if (is.null(coloc_out)) next
      
      key <- paste(
        region_chr,
        region_start,
        region_end,
        gene,
        sep = ":"
      )

      results[[key]] <- coloc_out
    }
  }

  results
}

# Run function
coloc_results <- run_coloc_per_region(regions, gwas, qtl)

###### END OF STEP 1 ######


####### STEP 2 - Function to extract results and save to a df #########
extract_coloc_summary <- function(coloc_results) {

  summaries <- lapply(names(coloc_results), function(name) {
    result <- coloc_results[[name]]
    if (is.null(result)) return(NULL)
    
    # Split region name into components
    parts <- strsplit(name, "[:]")[[1]]
    if (length(parts) < 4) {
      warning("Skipping malformed name: ", name)
      return(NULL)
    }
    
    chr <- parts[1]
    start <- parts[2]
    end <- parts[3]
    gene <- paste(parts[4:length(parts)], collapse = "-")
    
    # Extract H3 and H4 using correct names
    summary_vals <- result$summary
    if (is.null(summary_vals) || !all(c("PP.H3.abf", "PP.H4.abf") %in% names(summary_vals))) {
      warning("Missing H3 or H4 for: ", name)
      return(NULL)
    }
    
    PP.H3 <- summary_vals["PP.H3.abf"]
    PP.H4 <- summary_vals["PP.H4.abf"]
    
    lead_variant <- NA
    if (!is.null(result$results) && "SNP.PP.H4" %in% names(result$results)) {
      top_snp <- result$results[which.max(result$results$SNP.PP.H4), ]
      lead_variant <- top_snp$snp
    }
    
    return(data.frame(
      gene = gene,
      chr = chr,
      start = as.numeric(start),
      end = as.numeric(end),
      PP.H3 = as.numeric(PP.H3),
      PP.H4 = as.numeric(PP.H4),
      variant_id = lead_variant,
      stringsAsFactors = FALSE
    ))
  })
  
  coloc_summary_df <- do.call(rbind, summaries)
  return(coloc_summary_df)
}

coloc_summary_df <- extract_coloc_summary(coloc_results)

####### STEP 3 - Regional plot for significant results #########
plot_regional_coloc <- function(
  coloc_result, gwas_df, qtl_df,
  region, start, end, gene,
  h4_threshold = 0.8
) {

  # SNP causal
  causal_snp <- coloc_result$results$snp[
    which.max(coloc_result$results$SNP.PP.H4)
  ]

  # Subset GWAS / QTL da região
  gwas_sub <- gwas_df %>%
    dplyr::filter(chr == region, bp >= start, bp <= end) %>%
    dplyr::select(chr, bp, snp, pgwas = p)

  qtl_sub <- qtl_df %>%
    dplyr::filter(chr == region, bp >= start, bp <= end, symbol == gene) %>%
    dplyr::select(snp, p_eqtl = p)

  plot_df <- inner_join(
    gwas_sub,
    qtl_sub,
    by = "snp",
    suffix = c(".gwas", ".qtl")
  )

  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df <- as.data.frame(plot_df)

  # Região centrada no SNP causal
  center <- plot_df$bp[plot_df$snp == causal_snp][1]
  xrange <- c(start, end)

  if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) return(NULL)

  loc_gwas <- locus(
    data = plot_df,
    seqname = plot_df$chr[1],
    xrange = xrange,
    ens_db = "EnsDb.Hsapiens.v75",
    chrom = "chr",
    pos = "bp",
    p = "pgwas"
  ) |> link_LD(token = "5c4d1f5eeb21")

  loc_qtl <- locus(
    data = plot_df,
    seqname = plot_df$chr[1],
    xrange = xrange,
    ens_db = "EnsDb.Hsapiens.v75",
    chrom = "chr",
    pos = "bp",
    p = "p_eqtl"
  ) |> link_LD(token = "5c4d1f5eeb21")

  g  <- gg_genetracks(loc_gwas, highlight = gene)
  pg <- gg_scatter(loc_gwas, labels = causal_snp, nudge_x = 0.1,
  nudge_y = 0.1) + labs(title = "GWAS")
  pq <- gg_scatter(loc_qtl) + labs(title = paste0("QTL - ", gene))

  plot <- cowplot::plot_grid(pq, pg, g, ncol = 1, rel_heights = c(2, 2, 1), align = "v")

  outfile <- paste0(gene, "_", region, "_", start, "_", end, "_regional.png")

  ggsave(outfile, plot, width = 10, height = 8)
  message("Plot saved on ", outfile)

}

h4_threshold <- 0.8


coloc_to_plot <- coloc_summary_df %>%
  dplyr::filter(
    !is.na(PP.H4),
    PP.H4 >= h4_threshold
  )

for (i in seq_len(nrow(coloc_to_plot))) {

  row <- coloc_to_plot[i, ]

  res_key <- paste(
    row$chr,
    row$start,
    row$end,
    row$gene,
    sep = ":"
  )

  if (!res_key %in% names(coloc_results)) {
    message("Chave não encontrada: ", res_key)
    next
  }

  plot_regional_coloc(
    coloc_result = coloc_results[[res_key]],
    gwas_df = gwas,
    qtl_df  = qtl,
    region  = row$chr,
    start   = row$start,
    end     = row$end,
    gene    = row$gene
  )
}


write.csv(coloc_summary_df, paste0("coloc_summary.csv"), row.names = FALSE)

# Outputs desse módulo: coloc_summary.csv e *_regional.png
