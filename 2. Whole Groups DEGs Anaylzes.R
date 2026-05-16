# =====================================================================
# Integrated Bulk RNA-seq Analysis Pipeline
# Project: HN00273522 — BMDM × Macrophage × TNFα × Intestinal Organoid
# Design : ~ group, 8 levels, n=3 each (24 samples), 6-Tier 14-contrast
# Species: Mus musculus (mouse-only, consistent throughout)
# Modules: DESeq2 + GSEA .gct + EnhancedVolcano + clusterProfiler GO/KEGG
#        + pathview + gprofiler2 + limma + edgeR + WGCNA
# =====================================================================

# ===== 0. Packages ===================================================
# ---- Install if missing ----
# CRAN
cran_pkgs <- c("tidyverse","pheatmap","ggplot2","ggrepel","corrplot",
               "gprofiler2","pROC","ggforce","matrixStats","Hmisc",
               "foreach","doParallel","survival","dynamicTreeCut",
               "fastcluster","Rcpp","knitr","ggfortify","RColorBrewer",
               "patchwork","ashr","BiocManager")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
# Bioconductor
bioc_pkgs <- c("DESeq2","annotate","pathview","KEGGREST","org.Mm.eg.db",
               "gage","gageData","limma","apeglm","EnhancedVolcano",
               "edgeR","clusterProfiler","topGO","impute","preprocessCore",
               "enrichplot","AnnotationDbi","EnsDb.Mmusculus.v79","WGCNA")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

# ---- Load ----
suppressPackageStartupMessages({
  library(DESeq2);   library(ashr);   library(apeglm)
  library(ggplot2);  library(ggrepel); library(ggfortify)
  library(pheatmap); library(RColorBrewer); library(patchwork)
  library(corrplot); library(EnhancedVolcano)
  library(AnnotationDbi); library(org.Mm.eg.db); library(EnsDb.Mmusculus.v79)
  library(clusterProfiler); library(enrichplot); library(DOSE)
  library(gprofiler2)
  library(gage); library(gageData); library(pathview)
  library(limma); library(edgeR)
  library(WGCNA); library(matrixStats); library(dynamicTreeCut); library(fastcluster)
  library(dplyr); library(tidyr); library(tibble); library(knitr)
})

# WGCNA multi-thread
allowWGCNAThreads()

# ===== Configuration =================================================
# DEG cutoffs
padj_cutoff   <- 0.05
lfc_cutoff    <- 1.0          # log2(2)
fdr           <- 0.05

# Pre-filter
atleast_dds_count_for_PCA      <- 1
min_dds_count                   <- 10
min_dds_count_min_rowSums       <- 3

# Top genes for visualization
n_top_heatmap   <- 1000
n_top_oldPCA    <- 1000

# Default contrasts for "single-contrast" cross-validation analyses
# (limma, edgeR, gprofiler2, pathview)
default_contrasts <- c("T4_OrgTNFaBMDM5K_vs_OrgTNFa",
                       "T6_OrgTNFaBMDM5K_vs_MacTNFa")

# Random seed for WGCNA
wgcna_seed <- 6689

# Output root
out_root <- "Step1_Integrated_Pipeline"
dir.create(out_root, showWarnings = FALSE)
dir.create(file.path(out_root, "QC"),                showWarnings = FALSE)
dir.create(file.path(out_root, "DEG_results"),       showWarnings = FALSE)
dir.create(file.path(out_root, "Volcano"),           showWarnings = FALSE)
dir.create(file.path(out_root, "GO_KEGG"),           showWarnings = FALSE)
dir.create(file.path(out_root, "Pathview"),          showWarnings = FALSE)
dir.create(file.path(out_root, "gProfiler2"),        showWarnings = FALSE)
dir.create(file.path(out_root, "limma"),             showWarnings = FALSE)
dir.create(file.path(out_root, "edgeR"),             showWarnings = FALSE)
dir.create(file.path(out_root, "WGCNA"),             showWarnings = FALSE)
dir.create(file.path(out_root, "GSEA_export"),       showWarnings = FALSE)
dir.create(file.path(out_root, "tables"),            showWarnings = FALSE)

# 8-group color palette
group_colors_8 <- c(
  Organoid              = "#A6CEE3",
  Organoid_TNFa         = "#1F78B4",
  Organoid_BMDM_1K      = "#B2DF8A",
  Organoid_TNFa_BMDM_1K = "#33A02C",
  Organoid_BMDM_5K      = "#FDBF6F",
  Organoid_TNFa_BMDM_5K = "#FF7F00",
  Macrophage            = "#FB9A99",
  Macrophage_TNFa       = "#E31A1C"
)

# =====================================================================
# ===== 1. Load gene-level count matrix ===============================
# =====================================================================
cat("===== Loading gene-level count matrix =====\n")
data_raw <- read.csv("gene_count_matrix.csv", header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)

# First column = gene_id / gene_symbol → use as row names
gene_id_col <- colnames(data_raw)[1]
cat("Gene ID column: ", gene_id_col, "\n", sep = "")

frame <- data.frame(data_raw[, -1], row.names = data_raw[, 1],
                    check.names = FALSE)
cat("Raw matrix: ", nrow(frame), " genes × ", ncol(frame),
    " samples\n", sep = "")

# =====================================================================
# ===== 2. Sample mapping (24 samples → 8 full-name groups) ===========
# =====================================================================
sample_map <- c(
  Con1R = "Organoid_1",   Con2R = "Organoid_2",   Con3R = "Organoid_3",
  TNFa2 = "Organoid_TNFa_1", TNFa3 = "Organoid_TNFa_2", TNFa4 = "Organoid_TNFa_3",
  `1K1R`   = "Organoid_BMDM_1K_1",
  `1K2R`   = "Organoid_BMDM_1K_2",
  `1K3R`   = "Organoid_BMDM_1K_3",
  `5K1R`   = "Organoid_BMDM_5K_1",
  `5K2R`   = "Organoid_BMDM_5K_2",
  `5K3R`   = "Organoid_BMDM_5K_3",
  TNFa1K1  = "Organoid_TNFa_BMDM_1K_1",
  TNFa1K2  = "Organoid_TNFa_BMDM_1K_2",
  TNFa1K3  = "Organoid_TNFa_BMDM_1K_3",
  TNFa5K1  = "Organoid_TNFa_BMDM_5K_1",
  TNFa5K2  = "Organoid_TNFa_BMDM_5K_2",
  TNFa5K3  = "Organoid_TNFa_BMDM_5K_3",
  BMDM1    = "Macrophage_1",   BMDM2 = "Macrophage_2",   BMDM3 = "Macrophage_3",
  TNFaBMDM1 = "Macrophage_TNFa_1",
  TNFaBMDM2 = "Macrophage_TNFa_2",
  TNFaBMDM3 = "Macrophage_TNFa_3"
)
stopifnot(all(names(sample_map) %in% colnames(frame)))

counts <- frame[, names(sample_map), drop = FALSE]
colnames(counts) <- sample_map[colnames(counts)]
cat("After sample mapping: ", ncol(counts), " samples\n", sep = "")

# colData (8 groups, n=3 each)
group_levels <- c(
  "Organoid", "Organoid_TNFa",
  "Organoid_BMDM_1K", "Organoid_BMDM_5K",
  "Organoid_TNFa_BMDM_1K", "Organoid_TNFa_BMDM_5K",
  "Macrophage", "Macrophage_TNFa"
)
colData <- data.frame(
  sample = colnames(counts),
  group  = factor(sub("_[0-9]+$", "", colnames(counts)),
                  levels = group_levels),
  row.names = colnames(counts),
  stringsAsFactors = FALSE
)
stopifnot(all(!is.na(colData$group)))
stopifnot(identical(rownames(colData), colnames(counts)))

cat("\n===== Group × replicate balance =====\n")
print(table(colData$group))

# =====================================================================
# ===== 3. Build DESeqDataSet + pre-filtering =========================
# =====================================================================
counts_mat <- as.matrix(counts)
mode(counts_mat) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = colData,
                              design    = ~ group)
cat("\nDESeqDataSet created: ", nrow(dds), " genes × ",
    ncol(dds), " samples\n", sep = "")

# Pre-filter (≥10 counts in ≥3 samples)
keep <- rowSums(counts(dds) >= min_dds_count) >= min_dds_count_min_rowSums
cat("Pre-filter: ", nrow(dds), " → ", sum(keep),
    " (", round(100*mean(keep), 2), "% kept)\n", sep = "")
dds <- dds[keep, ]

# =====================================================================
# ===== 4. Density / Boxplot QC (raw vs filtered) =====================
# =====================================================================
cat("\n===== Generating raw vs filtered QC plots =====\n")
epsilon <- 1

# Sample colors per group
sample_color <- group_colors_8[as.character(colData$group)]

# Raw histogram + boxplot
png(file.path(out_root, "QC", "QC_raw_distribution.png"),
    width = 2200, height = 1400, res = 200)
par(mfrow = c(1, 2), mar = c(8, 5, 3, 2))
hist(log2(as.matrix(counts) + epsilon),
     main = "Raw counts (log2)", xlab = "log2(count + 1)",
     col = "grey80", border = "grey50")
boxplot(log2(as.matrix(counts) + epsilon),
        col = sample_color, las = 2, cex.axis = 0.6,
        main = "Raw counts (log2) by sample",
        ylab = "log2(count + 1)")
dev.off()

# Filtered (low-count removed)
filt_counts <- as.matrix(counts[keep, ])
png(file.path(out_root, "QC", "QC_filtered_distribution.png"),
    width = 2200, height = 1400, res = 200)
par(mfrow = c(1, 2), mar = c(8, 5, 3, 2))
hist(log2(filt_counts + epsilon),
     main = "After pre-filter (≥10 in ≥3 samples)",
     xlab = "log2(count + 1)", col = "grey80", border = "grey50")
boxplot(log2(filt_counts + epsilon),
        col = sample_color, las = 2, cex.axis = 0.6,
        main = "Pre-filtered counts (log2) by sample",
        ylab = "log2(count + 1)")
dev.off()

# Detected gene count per sample
png(file.path(out_root, "QC", "QC_detected_genes_per_sample.png"),
    width = 1800, height = 1200, res = 200)
par(mar = c(8, 5, 3, 2))
detected <- colSums(filt_counts > 1)
barplot(detected, las = 2, col = sample_color, cex.names = 0.6,
        ylab = "Detected genes (count > 1)",
        main = "Detected genes per sample (after pre-filter)")
abline(h = median(detected), col = "red", lty = 2)
dev.off()

# =====================================================================
# ===== 5. DESeq2 run + VST + PCA + sample distance + dispersion =====
# =====================================================================
cat("\n===== Running DESeq2 =====\n")
dds <- DESeq(dds)
cat("resultsNames:\n"); print(resultsNames(dds))

# Normalized counts (for boxplot QC + GSEA export)
normalized_counts_mat <- counts(dds, normalized = TRUE)

# Normalized boxplot
png(file.path(out_root, "QC", "QC_normalized_boxplot.png"),
    width = 2200, height = 1400, res = 200)
par(mar = c(8, 5, 3, 2))
boxplot(log2(normalized_counts_mat + epsilon),
        col = sample_color, las = 2, cex.axis = 0.6,
        main = "DESeq2-normalized counts (log2)",
        ylab = "log2(normalized count + 1)")
abline(h = median(log2(normalized_counts_mat + epsilon)),
       col = "blue", lty = 2)
dev.off()

# VST
vsd <- vst(dds, blind = FALSE)

# PCA
pca_df <- DESeq2::plotPCA(vsd, intgroup = "group", returnData = TRUE)
pv     <- round(100 * attr(pca_df, "percentVar"))
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 2.6, max.overlaps = 30, show.legend = FALSE) +
  scale_color_manual(values = group_colors_8) +
  xlab(paste0("PC1: ", pv[1], "% variance")) +
  ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw(base_size = 12) +
  ggtitle("PCA — VST normalized (8 groups, 24 samples)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(out_root, "QC", "QC_PCA.png"),
       p_pca, width = 9, height = 6, dpi = 300)

# Sample distance heatmap
sd_mat <- as.matrix(dist(t(assay(vsd))))
ann_col <- data.frame(group = colData$group, row.names = rownames(colData))
pheatmap(sd_mat,
         clustering_distance_rows = as.dist(sd_mat),
         clustering_distance_cols = as.dist(sd_mat),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         annotation_col = ann_col,
         annotation_colors = list(group = group_colors_8),
         fontsize = 7,
         filename = file.path(out_root, "QC", "QC_sample_distance.png"),
         width = 9, height = 7.5)

# Dispersion plot
png(file.path(out_root, "QC", "QC_dispersion.png"),
    width = 1800, height = 1500, res = 300)
plotDispEsts(dds, main = "Dispersion estimates (gene-level, n=24)")
dev.off()

# =====================================================================
# ===== 6. Contrast definitions (6 tiers, 14 contrasts) ===============
# =====================================================================
contrasts_tier1 <- list(
  T1_OrgBMDM1K_vs_Org = c("group", "Organoid_BMDM_1K", "Organoid"),
  T1_OrgBMDM5K_vs_Org = c("group", "Organoid_BMDM_5K", "Organoid")
)
contrasts_tier2 <- list(
  T2_OrgBMDM1K_vs_Mac = c("group", "Organoid_BMDM_1K", "Macrophage"),
  T2_OrgBMDM5K_vs_Mac = c("group", "Organoid_BMDM_5K", "Macrophage")
)
contrasts_tier3 <- list(
  T3_OrgTNFa_vs_Org           = c("group", "Organoid_TNFa",           "Organoid"),
  T3_OrgTNFaBMDM1K_vs_Org     = c("group", "Organoid_TNFa_BMDM_1K",   "Organoid"),
  T3_OrgTNFaBMDM5K_vs_Org     = c("group", "Organoid_TNFa_BMDM_5K",   "Organoid")
)
contrasts_tier4 <- list(
  T4_OrgTNFaBMDM1K_vs_OrgTNFa = c("group", "Organoid_TNFa_BMDM_1K", "Organoid_TNFa"),
  T4_OrgTNFaBMDM5K_vs_OrgTNFa = c("group", "Organoid_TNFa_BMDM_5K", "Organoid_TNFa")
)
contrasts_tier5 <- list(
  T5_MacTNFa_vs_Mac           = c("group", "Macrophage_TNFa",         "Macrophage"),
  T5_OrgTNFaBMDM1K_vs_Mac     = c("group", "Organoid_TNFa_BMDM_1K",   "Macrophage"),
  T5_OrgTNFaBMDM5K_vs_Mac     = c("group", "Organoid_TNFa_BMDM_5K",   "Macrophage")
)
contrasts_tier6 <- list(
  T6_OrgTNFaBMDM1K_vs_MacTNFa = c("group", "Organoid_TNFa_BMDM_1K", "Macrophage_TNFa"),
  T6_OrgTNFaBMDM5K_vs_MacTNFa = c("group", "Organoid_TNFa_BMDM_5K", "Macrophage_TNFa")
)

all_contrasts <- c(contrasts_tier1, contrasts_tier2, contrasts_tier3,
                   contrasts_tier4, contrasts_tier5, contrasts_tier6)

contrast_tiers <- data.frame(
  contrast = c(names(contrasts_tier1), names(contrasts_tier2),
               names(contrasts_tier3), names(contrasts_tier4),
               names(contrasts_tier5), names(contrasts_tier6)),
  tier = c(rep("Tier1_basal_OrgRef",   length(contrasts_tier1)),
           rep("Tier2_basal_MacRef",   length(contrasts_tier2)),
           rep("Tier3_TNFa_OrgAbs",    length(contrasts_tier3)),
           rep("Tier4_TNFa_OrgTNFa",   length(contrasts_tier4)),
           rep("Tier5_TNFa_MacAbs",    length(contrasts_tier5)),
           rep("Tier6_TNFa_MacTNFa",   length(contrasts_tier6))),
  control_group = c(rep("Organoid",          length(contrasts_tier1)),
                    rep("Macrophage",        length(contrasts_tier2)),
                    rep("Organoid",          length(contrasts_tier3)),
                    rep("Organoid_TNFa",     length(contrasts_tier4)),
                    rep("Macrophage",        length(contrasts_tier5)),
                    rep("Macrophage_TNFa",   length(contrasts_tier6))),
  stringsAsFactors = FALSE
)
write.csv(contrast_tiers, file.path(out_root, "contrast_tiers.csv"),
          row.names = FALSE)

cat("\n===== 6-Tier Contrast structure (14 contrasts) =====\n")
print(contrast_tiers, row.names = FALSE)

# =====================================================================
# ===== 7. Gene annotation (mouse) — handles "ENSMUSG.version|SYMBOL"
#         and "MSTRG.X|SYMBOL" composite IDs
# =====================================================================
cat("\n===== Building mouse gene annotation =====\n")

all_gene_ids_raw <- rownames(dds)
cat("First 5 raw gene IDs:\n")
print(head(all_gene_ids_raw, 5))

# ---- Parse composite IDs: "PREFIX.version|SYMBOL" ----
parse_id <- function(id) {
  parts <- strsplit(id, "\\|", fixed = FALSE)[[1]]
  ensembl <- NA_character_
  symbol  <- NA_character_
  mstrg   <- NA_character_
  
  for (p in parts) {
    if (grepl("^ENSMUSG", p)) {
      ensembl <- sub("\\..*$", "", p)  # strip version
    } else if (grepl("^MSTRG\\.", p)) {
      mstrg <- p
    } else {
      symbol <- p  # everything else assumed to be SYMBOL
    }
  }
  c(raw = id, ensembl = ensembl, symbol = symbol, mstrg = mstrg)
}

parsed_mat <- t(sapply(all_gene_ids_raw, parse_id))
parsed_df  <- as.data.frame(parsed_mat, stringsAsFactors = FALSE)
rownames(parsed_df) <- NULL

cat("\n===== ID composition =====\n")
cat("  Total IDs                : ", nrow(parsed_df), "\n", sep = "")
cat("  With ENSEMBL             : ", sum(!is.na(parsed_df$ensembl)),
    " (", round(100*mean(!is.na(parsed_df$ensembl)), 1), "%)\n", sep = "")
cat("  With SYMBOL              : ", sum(!is.na(parsed_df$symbol)),
    " (", round(100*mean(!is.na(parsed_df$symbol)), 1), "%)\n", sep = "")
cat("  MSTRG (StringTie novel)  : ", sum(!is.na(parsed_df$mstrg)),
    " (", round(100*mean(!is.na(parsed_df$mstrg)), 1), "%)\n", sep = "")
cat("    of which with SYMBOL   : ",
    sum(!is.na(parsed_df$mstrg) & !is.na(parsed_df$symbol)), "\n", sep = "")
cat("    of which without SYMBOL: ",
    sum(!is.na(parsed_df$mstrg) & is.na(parsed_df$symbol)), "\n", sep = "")

# ---- FILTER: remove MSTRG-only entries (no SYMBOL, no ENSEMBL) ----
keep_anno <- !is.na(parsed_df$symbol) | !is.na(parsed_df$ensembl)
cat("\n===== Filtering =====\n")
cat("  Removing MSTRG-only entries (no annotation): ",
    sum(!keep_anno), "\n", sep = "")
cat("  Keeping: ", sum(keep_anno), " genes\n", sep = "")

parsed_df <- parsed_df[keep_anno, ]
dds       <- dds[parsed_df$raw, ]   # subset dds to match
cat("  dds subset to: ", nrow(dds), " genes × ", ncol(dds),
    " samples\n", sep = "")

# ---- Query org.Mm.eg.db ----
# Use SYMBOL as primary key (faster + cleaner), fall back to ENSEMBL
valid_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")
valid_ensembl <- keys(org.Mm.eg.db, keytype = "ENSEMBL")

# Query SYMBOL (primary)
symbols_to_query <- unique(na.omit(parsed_df$symbol))
symbols_to_query <- intersect(symbols_to_query, valid_symbols)
anno_sym <- if (length(symbols_to_query) > 0) {
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = symbols_to_query,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID", "GENENAME")
  )
} else NULL

# Query ENSEMBL (for those without SYMBOL match)
ensembl_to_query <- unique(na.omit(parsed_df$ensembl[is.na(parsed_df$symbol) |
                                                       !parsed_df$symbol %in% symbols_to_query]))
ensembl_to_query <- intersect(ensembl_to_query, valid_ensembl)
anno_ens <- if (length(ensembl_to_query) > 0) {
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = ensembl_to_query,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL", "ENTREZID", "GENENAME")
  )
} else NULL

# ---- Merge annotations back to parsed_df ----
anno <- parsed_df
colnames(anno) <- c("gene_id", "ENSEMBL", "SYMBOL_parsed", "MSTRG")

# Merge SYMBOL-based annotation
if (!is.null(anno_sym)) {
  anno_sym <- anno_sym[!duplicated(anno_sym$SYMBOL), ]
  anno <- merge(anno, anno_sym,
                by.x = "SYMBOL_parsed", by.y = "SYMBOL",
                all.x = TRUE, sort = FALSE)
  # Keep SYMBOL column
  anno$SYMBOL <- anno$SYMBOL_parsed
} else {
  anno$ENTREZID  <- NA_character_
  anno$GENENAME  <- NA_character_
  anno$SYMBOL    <- anno$SYMBOL_parsed
}

# Fill missing ENTREZID via ENSEMBL lookup
if (!is.null(anno_ens)) {
  anno_ens <- anno_ens[!duplicated(anno_ens$ENSEMBL), ]
  fill_idx <- which(is.na(anno$ENTREZID) & !is.na(anno$ENSEMBL))
  if (length(fill_idx) > 0) {
    match_idx <- match(anno$ENSEMBL[fill_idx], anno_ens$ENSEMBL)
    anno$ENTREZID[fill_idx] <- anno_ens$ENTREZID[match_idx]
    # Also fill SYMBOL if missing
    sym_fill <- is.na(anno$SYMBOL[fill_idx])
    anno$SYMBOL[fill_idx][sym_fill] <- anno_ens$SYMBOL[match_idx][sym_fill]
    # Also fill GENENAME if missing
    gn_fill <- is.na(anno$GENENAME[fill_idx])
    anno$GENENAME[fill_idx][gn_fill] <- anno_ens$GENENAME[match_idx][gn_fill]
  }
}

# Reorder columns; drop helper
anno <- anno[, c("gene_id", "SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME", "MSTRG")]

# id_type for downstream code compatibility — we use "gene_id" everywhere
id_type <- "gene_id"

cat("\n===== Final annotation coverage =====\n")
cat("  Total genes in dds  : ", nrow(dds), "\n", sep = "")
cat("  In annotation table : ", nrow(anno), "\n", sep = "")
cat("  With SYMBOL         : ", sum(!is.na(anno$SYMBOL)),
    " (", round(100*mean(!is.na(anno$SYMBOL)), 1), "%)\n", sep = "")
cat("  With ENTREZID       : ", sum(!is.na(anno$ENTREZID)),
    " (", round(100*mean(!is.na(anno$ENTREZID)), 1), "%)\n", sep = "")
cat("  With ENSEMBL        : ", sum(!is.na(anno$ENSEMBL)),
    " (", round(100*mean(!is.na(anno$ENSEMBL)), 1), "%)\n", sep = "")
cat("  Of which MSTRG-tagged with SYMBOL: ",
    sum(!is.na(anno$MSTRG) & !is.na(anno$SYMBOL)), "\n", sep = "")

# Save
saveRDS(anno, "gene_annotation.rds")
write.csv(anno, file.path(out_root, "tables", "gene_annotation.csv"),
          row.names = FALSE)
# =====================================================================
# ===== 8. DEG extraction (ashr shrinkage) - all 14 contrasts ========
# =====================================================================
extract_deg <- function(dds, contrast_vec, anno_df, id_type = "SYMBOL") {
  res     <- results(dds, contrast = contrast_vec, alpha = fdr)
  res_shr <- lfcShrink(dds, contrast = contrast_vec,
                       res = res, type = "ashr", quiet = TRUE)
  out <- as.data.frame(res_shr)
  out$gene_id <- rownames(out)
  
  # Merge annotation
  out <- merge(out, anno_df, by.x = "gene_id", by.y = id_type, all.x = TRUE)
  out <- out[order(out$padj, -abs(out$log2FoldChange)), ]
  
  # Reorder columns
  base_cols <- c("gene_id", "SYMBOL", "ENTREZID", "GENENAME",
                 "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
  if (id_type == "SYMBOL") base_cols <- setdiff(base_cols, "SYMBOL")
  if (id_type == "ENSEMBL") base_cols <- c(base_cols, "ENSEMBL")
  
  out[, intersect(base_cols, colnames(out))]
}

cat("\n===== Running DEG extraction for all 14 contrasts =====\n")
deg_results <- list()
for (nm in names(all_contrasts)) {
  cat("  Running: ", nm, "\n", sep = "")
  deg_results[[nm]] <- extract_deg(dds, all_contrasts[[nm]], anno, id_type)
}

# =====================================================================
# ===== 9. DEG summary + export =======================================
# =====================================================================
deg_summary <- do.call(rbind, lapply(names(deg_results), function(nm) {
  r <- deg_results[[nm]]
  ti <- contrast_tiers[contrast_tiers$contrast == nm, ]
  data.frame(
    contrast       = nm,
    tier           = ti$tier,
    control_group  = ti$control_group,
    total_tested   = sum(!is.na(r$padj)),
    sig_padj05     = sum(r$padj < padj_cutoff, na.rm = TRUE),
    up_lfc1        = sum(r$padj < padj_cutoff & r$log2FoldChange >  lfc_cutoff, na.rm = TRUE),
    dn_lfc1        = sum(r$padj < padj_cutoff & r$log2FoldChange < -lfc_cutoff, na.rm = TRUE),
    up_lfc2        = sum(r$padj < padj_cutoff & r$log2FoldChange >  2, na.rm = TRUE),
    dn_lfc2        = sum(r$padj < padj_cutoff & r$log2FoldChange < -2, na.rm = TRUE)
  )
}))
cat("\n===== DEG summary =====\n")
print(deg_summary, row.names = FALSE)
write.csv(deg_summary, file.path(out_root, "tables", "DEG_summary.csv"),
          row.names = FALSE)

# Per-contrast tables (full + up + down)
for (nm in names(deg_results)) {
  tier <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  r    <- deg_results[[nm]]
  
  # Full
  write.csv(r,
            file.path(out_root, "DEG_results",
                      paste0(tier, "__", nm, "_FULL.csv")),
            row.names = FALSE)
  # Up
  r_up <- r %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    dplyr::arrange(padj)
  write.csv(r_up,
            file.path(out_root, "DEG_results",
                      paste0(tier, "__", nm, "_UP.csv")),
            row.names = FALSE)
  # Down
  r_dn <- r %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    dplyr::arrange(padj)
  write.csv(r_dn,
            file.path(out_root, "DEG_results",
                      paste0(tier, "__", nm, "_DOWN.csv")),
            row.names = FALSE)
}

# =====================================================================
# ===== 10. GSEA .gct export (normalized counts) ======================
# =====================================================================
cat("\n===== Exporting GSEA-formatted files =====\n")

# Get ENTREZID for each gene
gsea_anno <- anno[match(rownames(normalized_counts_mat), anno[[id_type]]), ]
gsea_df <- data.frame(
  NAME        = gsea_anno$ENTREZID,
  Description = rownames(normalized_counts_mat),
  normalized_counts_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
# Remove rows without NAME (ENTREZID)
gsea_df <- gsea_df[!is.na(gsea_df$NAME) & gsea_df$NAME != "", ]
gsea_df <- gsea_df[!duplicated(gsea_df$NAME), ]
cat("GSEA matrix: ", nrow(gsea_df), " genes × ",
    ncol(normalized_counts_mat), " samples\n", sep = "")

# Write .txt (tab-separated, normalized counts for external GSEA)
write.table(gsea_df,
            file.path(out_root, "GSEA_export",
                      "normalized_counts_for_GSEA.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Also export .gct format header (Broad GSEA)
gct_file <- file.path(out_root, "GSEA_export",
                      "normalized_counts_for_GSEA.gct")
writeLines(c("#1.2",
             paste0(nrow(gsea_df), "\t", ncol(normalized_counts_mat))),
           con = gct_file)
write.table(gsea_df, gct_file, sep = "\t", quote = FALSE,
            row.names = FALSE, append = TRUE, col.names = TRUE)

# Sample .cls file (categorical, by group)
cls_file <- file.path(out_root, "GSEA_export", "groups_for_GSEA.cls")
n_samp <- nrow(colData); n_grp <- length(levels(colData$group))
writeLines(c(paste(n_samp, n_grp, "1", sep = " "),
             paste("#", paste(levels(colData$group), collapse = " ")),
             paste(as.integer(colData$group) - 1, collapse = " ")),
           con = cls_file)

cat("GSEA files exported to: ", file.path(out_root, "GSEA_export"), "\n", sep = "")

# =====================================================================
# ===== 11. EnhancedVolcano (all 14 contrasts) ========================
# =====================================================================
cat("\n===== Generating EnhancedVolcano plots =====\n")

volcano_legend <- c(
  "NS",
  expression(paste("Log"[2], "FC")),
  expression(italic(p) - value),
  expression(italic(p) ~ - value ~ and ~ Log[2] ~ FC)
)

for (nm in names(deg_results)) {
  tier <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  ctrl <- contrast_tiers$control_group[contrast_tiers$contrast == nm]
  r    <- deg_results[[nm]]
  
  # Need rownames for EnhancedVolcano
  rdf <- r
  rownames(rdf) <- make.unique(as.character(rdf$gene_id))
  
  fname <- file.path(out_root, "Volcano",
                     paste0("Volcano_", tier, "__", nm, ".png"))
  
  png(fname, width = 2400, height = 2600, res = 280)
  print(
    EnhancedVolcano(rdf,
                    lab = rownames(rdf),
                    x = "log2FoldChange",
                    y = "pvalue",
                    pCutoff = padj_cutoff,
                    FCcutoff = lfc_cutoff,
                    legendLabels = volcano_legend,
                    drawConnectors = TRUE,
                    widthConnectors = 0.2,
                    colConnectors = "grey30",
                    title = paste0(nm),
                    subtitle = paste0(tier, "  |  Control = ", ctrl),
                    caption = sprintf("padj < %.2f, |LFC| > %.1f", padj_cutoff, lfc_cutoff),
                    pointSize = 1.5, labSize = 3.0
    )
  )
  dev.off()
}
cat("  Saved ", length(deg_results), " EnhancedVolcano plots\n", sep = "")

# =====================================================================
# ===== 12. clusterProfiler GO + KEGG (all 14 contrasts) =============
# =====================================================================
cat("\n===== Running clusterProfiler GO + KEGG =====\n")

universe_entrez <- unique(na.omit(as.character(anno$ENTREZID)))

run_go <- function(deg_df, ont = "BP") {
  sig_genes <- deg_df %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::pull(ENTREZID) %>%
    na.omit() %>% unique() %>% as.character()
  if (length(sig_genes) < 10) return(NULL)
  
  tryCatch(
    enrichGO(gene = sig_genes,
             universe = universe_entrez,
             OrgDb = org.Mm.eg.db,
             keyType = "ENTREZID",
             ont = ont,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.10,
             readable = TRUE),
    error = function(e) { message("GO failed: ", e$message); NULL }
  )
}

run_kegg <- function(deg_df) {
  sig_genes <- deg_df %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::pull(ENTREZID) %>%
    na.omit() %>% unique() %>% as.character()
  if (length(sig_genes) < 10) return(NULL)
  
  res <- tryCatch(
    enrichKEGG(gene = sig_genes,
               organism = "mmu",
               universe = universe_entrez,
               keyType = "kegg",
               pvalueCutoff = 0.1,
               qvalueCutoff = 0.25),
    error = function(e) { message("KEGG failed: ", e$message); NULL }
  )
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    res <- setReadable(res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }
  res
}

go_results   <- list()
kegg_results <- list()

for (nm in names(deg_results)) {
  cat("  GO/KEGG: ", nm, "\n", sep = "")
  go_results[[nm]] <- list(
    BP = run_go(deg_results[[nm]], "BP"),
    CC = run_go(deg_results[[nm]], "CC"),
    MF = run_go(deg_results[[nm]], "MF")
  )
  kegg_results[[nm]] <- run_kegg(deg_results[[nm]])
}

# Export tables + dotplots
for (nm in names(deg_results)) {
  tier <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  
  # GO
  for (ont in c("BP", "CC", "MF")) {
    r <- go_results[[nm]][[ont]]
    if (is.null(r) || nrow(as.data.frame(r)) == 0) next
    write.csv(as.data.frame(r),
              file.path(out_root, "GO_KEGG",
                        paste0("GO_", ont, "__", tier, "__", nm, ".csv")),
              row.names = FALSE)
    p <- dotplot(r, showCategory = 15, font.size = 9) +
      ggtitle(paste0("GO ", ont, " — ", nm, "\n", tier)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(out_root, "GO_KEGG",
                     paste0("GO_", ont, "__", tier, "__", nm, ".png")),
           p, width = 8, height = 7, dpi = 300)
  }
  
  # KEGG
  r <- kegg_results[[nm]]
  if (!is.null(r) && nrow(as.data.frame(r)) > 0) {
    write.csv(as.data.frame(r),
              file.path(out_root, "GO_KEGG",
                        paste0("KEGG__", tier, "__", nm, ".csv")),
              row.names = FALSE)
    p <- dotplot(r, showCategory = 15, font.size = 9) +
      ggtitle(paste0("KEGG — ", nm, "\n", tier)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(out_root, "GO_KEGG",
                     paste0("KEGG__", tier, "__", nm, ".png")),
           p, width = 8, height = 7, dpi = 300)
  }
}

# =====================================================================
# ===== 13. Pathview — LFC overlay on KEGG pathways (Tier 4, 6) ======
# =====================================================================
cat("\n===== Running pathview for default contrasts =====\n")

run_pathview_for_contrast <- function(contrast_name, top_n_pathways = 5) {
  if (is.null(kegg_results[[contrast_name]])) {
    cat("  No KEGG result for ", contrast_name, "\n"); return(invisible(NULL))
  }
  r <- as.data.frame(kegg_results[[contrast_name]])
  if (nrow(r) == 0) return(invisible(NULL))
  
  tier <- contrast_tiers$tier[contrast_tiers$contrast == contrast_name]
  pv_dir <- file.path(out_root, "Pathview", paste0(tier, "__", contrast_name))
  dir.create(pv_dir, showWarnings = FALSE, recursive = TRUE)
  
  # LFC named vector (ENTREZID)
  d <- deg_results[[contrast_name]]
  fc <- d$log2FoldChange
  names(fc) <- d$ENTREZID
  fc <- fc[!is.na(names(fc)) & !duplicated(names(fc))]
  
  # Top N enriched pathways
  top_ids <- head(r$ID, top_n_pathways)
  
  # Run pathview (writes PNG to pv_dir)
  old_wd <- getwd()
  setwd(pv_dir)
  for (pid in top_ids) {
    tryCatch(
      pathview(gene.data = fc, pathway.id = sub("^mmu", "", pid),
               species = "mmu", limit = list(gene = 3, cpd = 1)),
      error = function(e) message("    pathview failed for ", pid, ": ", e$message)
    )
  }
  setwd(old_wd)
  cat("  Pathview done: ", contrast_name, " (", length(top_ids),
      " pathways)\n", sep = "")
}

for (cn in default_contrasts) run_pathview_for_contrast(cn, top_n_pathways = 5)

# =====================================================================
# ===== 14. gprofiler2 cross-validation (default contrasts) ==========
# =====================================================================
cat("\n===== Running gprofiler2 for default contrasts =====\n")

run_gprofiler2 <- function(contrast_name) {
  d <- deg_results[[contrast_name]]
  sig_genes <- d %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::pull(gene_id) %>%
    unique() %>% as.character()
  if (length(sig_genes) < 10) {
    cat("  gprofiler2 skipped (<10 sig genes): ", contrast_name, "\n")
    return(invisible(NULL))
  }
  
  tier <- contrast_tiers$tier[contrast_tiers$contrast == contrast_name]
  
  gres <- tryCatch(
    gost(query = sig_genes,
         organism = "mmusculus",
         user_threshold = 0.05,
         correction_method = "g_SCS",
         sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC","WP")),
    error = function(e) { message("gprofiler2 failed: ", e$message); NULL }
  )
  if (is.null(gres)) return(invisible(NULL))
  
  # Save table
  write.csv(gres$result,
            file.path(out_root, "gProfiler2",
                      paste0("gost__", tier, "__", contrast_name, ".csv")),
            row.names = FALSE)
  
  # Manhattan plot
  png(file.path(out_root, "gProfiler2",
                paste0("gost_manhattan__", tier, "__", contrast_name, ".png")),
      width = 2400, height = 1600, res = 200)
  print(gostplot(gres, capped = TRUE, interactive = FALSE))
  dev.off()
  
  invisible(gres)
}

gprofiler2_results <- list()
for (cn in default_contrasts) {
  gprofiler2_results[[cn]] <- run_gprofiler2(cn)
}

# =====================================================================
# ===== 15. limma + voom cross-validation (default contrasts) ========
# =====================================================================
cat("\n===== Running limma + voom for default contrasts =====\n")

run_limma_voom <- function(contrast_name, raw_counts, colData) {
  cv <- all_contrasts[[contrast_name]]
  level_treat <- cv[2]
  level_ctrl  <- cv[3]
  
  # Subset samples to only the two groups being compared
  keep_samp <- colData$group %in% c(level_treat, level_ctrl)
  cnt   <- raw_counts[, keep_samp, drop = FALSE]
  cdata <- colData[keep_samp, , drop = FALSE]
  cdata$group <- droplevels(cdata$group)
  
  # edgeR DGEList for voom
  dge <- DGEList(counts = cnt, group = cdata$group)
  keep <- filterByExpr(dge, group = cdata$group)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "TMM")
  
  design <- model.matrix(~ 0 + group, data = cdata)
  colnames(design) <- levels(cdata$group)
  
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  
  contrast_str <- paste0(level_treat, " - ", level_ctrl)
  cont_mat <- makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- contrasts.fit(fit, cont_mat)
  fit2 <- eBayes(fit2)
  
  tt <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH",
                 sort.by = "P")
  tt$gene_id <- rownames(tt)
  tt <- merge(tt, anno, by.x = "gene_id", by.y = id_type, all.x = TRUE)
  tt <- tt[order(tt$adj.P.Val), ]
  
  tier <- contrast_tiers$tier[contrast_tiers$contrast == contrast_name]
  write.csv(tt,
            file.path(out_root, "limma",
                      paste0("limma_voom__", tier, "__", contrast_name, ".csv")),
            row.names = FALSE)
  
  # Volcano (base R)
  png(file.path(out_root, "limma",
                paste0("limma_volcano__", tier, "__", contrast_name, ".png")),
      width = 2000, height = 2000, res = 250)
  par(mar = c(5, 5, 4, 2))
  with(tt, plot(logFC, -log10(P.Value), pch = 20, col = "grey70",
                main = paste0("limma+voom — ", contrast_name),
                xlab = "log2FC", ylab = "-log10(P)"))
  with(subset(tt, adj.P.Val < padj_cutoff & logFC >  lfc_cutoff),
       points(logFC, -log10(P.Value), pch = 20, col = "red"))
  with(subset(tt, adj.P.Val < padj_cutoff & logFC < -lfc_cutoff),
       points(logFC, -log10(P.Value), pch = 20, col = "blue"))
  abline(h = -log10(0.05), col = "black", lty = 2)
  abline(v = c(-lfc_cutoff, lfc_cutoff), col = "black", lty = 2)
  dev.off()
  
  cat("  limma done: ", contrast_name, " | ",
      sum(tt$adj.P.Val < padj_cutoff & abs(tt$logFC) > lfc_cutoff, na.rm = TRUE),
      " sig DEGs\n", sep = "")
  
  return(tt)
}

limma_results <- list()
for (cn in default_contrasts) {
  limma_results[[cn]] <- run_limma_voom(cn, counts_mat, colData)
}

# =====================================================================
# ===== 16. edgeR + glmQLFit cross-validation (default contrasts) ====
# =====================================================================
cat("\n===== Running edgeR + glmQLFit for default contrasts =====\n")

run_edgeR <- function(contrast_name, raw_counts, colData) {
  cv <- all_contrasts[[contrast_name]]
  level_treat <- cv[2]
  level_ctrl  <- cv[3]
  
  keep_samp <- colData$group %in% c(level_treat, level_ctrl)
  cnt   <- raw_counts[, keep_samp, drop = FALSE]
  cdata <- colData[keep_samp, , drop = FALSE]
  cdata$group <- droplevels(cdata$group)
  # Re-level so that ctrl is reference
  cdata$group <- relevel(cdata$group, ref = level_ctrl)
  
  dge <- DGEList(counts = cnt, group = cdata$group)
  keep <- filterByExpr(dge, group = cdata$group)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- calcNormFactors(dge, method = "TMM")
  
  design <- model.matrix(~ group, data = cdata)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  
  # Treatment with LFC threshold (more conservative)
  tr <- glmTreat(fit, coef = 2, lfc = lfc_cutoff)
  tab <- topTags(tr, n = Inf, adjust.method = "BH")$table
  tab$gene_id <- rownames(tab)
  tab <- merge(tab, anno, by.x = "gene_id", by.y = id_type, all.x = TRUE)
  tab <- tab[order(tab$FDR), ]
  
  tier <- contrast_tiers$tier[contrast_tiers$contrast == contrast_name]
  write.csv(tab,
            file.path(out_root, "edgeR",
                      paste0("edgeR_glmTreat__", tier, "__", contrast_name, ".csv")),
            row.names = FALSE)
  
  # Volcano (base R)
  png(file.path(out_root, "edgeR",
                paste0("edgeR_volcano__", tier, "__", contrast_name, ".png")),
      width = 2000, height = 2000, res = 250)
  par(mar = c(5, 5, 4, 2))
  with(tab, plot(logFC, -log10(PValue), pch = 20, col = "grey70",
                 main = paste0("edgeR glmTreat — ", contrast_name),
                 xlab = "log2FC", ylab = "-log10(P)"))
  with(subset(tab, FDR < padj_cutoff & logFC >  lfc_cutoff),
       points(logFC, -log10(PValue), pch = 20, col = "red"))
  with(subset(tab, FDR < padj_cutoff & logFC < -lfc_cutoff),
       points(logFC, -log10(PValue), pch = 20, col = "blue"))
  abline(h = -log10(0.05), col = "black", lty = 2)
  abline(v = c(-lfc_cutoff, lfc_cutoff), col = "black", lty = 2)
  dev.off()
  
  cat("  edgeR done: ", contrast_name, " | ",
      sum(tab$FDR < padj_cutoff, na.rm = TRUE), " sig DEGs (glmTreat)\n",
      sep = "")
  
  return(tab)
}

edger_results <- list()
for (cn in default_contrasts) {
  edger_results[[cn]] <- run_edgeR(cn, counts_mat, colData)
}

# DESeq2 vs limma vs edgeR overlap (Venn-like summary)
cat("\n===== Cross-platform DEG overlap (default contrasts) =====\n")
overlap_summary <- list()
for (cn in default_contrasts) {
  d_de  <- deg_results[[cn]] %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::pull(gene_id) %>% unique()
  d_li  <- limma_results[[cn]] %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < padj_cutoff,
                  abs(logFC) > lfc_cutoff) %>%
    dplyr::pull(gene_id) %>% unique()
  d_eg  <- edger_results[[cn]] %>%
    dplyr::filter(!is.na(FDR), FDR < padj_cutoff) %>%
    dplyr::pull(gene_id) %>% unique()
  
  all3 <- intersect(intersect(d_de, d_li), d_eg)
  overlap_summary[[cn]] <- data.frame(
    contrast = cn,
    DESeq2   = length(d_de),
    limma    = length(d_li),
    edgeR    = length(d_eg),
    all_three = length(all3)
  )
}
overlap_df <- do.call(rbind, overlap_summary)
print(overlap_df, row.names = FALSE)
write.csv(overlap_df,
          file.path(out_root, "tables", "DEG_overlap_3platforms.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 17. WGCNA (co-expression network) =============================
# =====================================================================
cat("\n===== Running WGCNA =====\n")
cat("(This may take 5-15 min depending on sample size and gene count)\n")

# Build a separate dds for WGCNA (no design model needed)
dds_wgcna <- DESeqDataSetFromMatrix(countData = counts_mat,
                                    colData   = colData,
                                    design    = ~ 1)
keep_w <- rowSums(counts(dds_wgcna) >= min_dds_count) >= min_dds_count_min_rowSums
dds_wgcna <- dds_wgcna[keep_w, ]

# VST + transpose (WGCNA expects samples as rows)
vsd_wgcna <- vst(dds_wgcna, blind = TRUE)
wgcna_expr <- t(assay(vsd_wgcna))

# Determine soft threshold
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(wgcna_expr, powerVector = powers,
                         dataIsExpr = TRUE, corFnc = cor,
                         networkType = "signed", verbose = 0)
saveRDS(sft, file.path(out_root, "WGCNA", "sft.rds"))

# Plot scale independence & mean connectivity
png(file.path(out_root, "WGCNA", "WGCNA_threshold_selection.png"),
    width = 2400, height = 1200, res = 200)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "red", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = 0.9, col = "red")
abline(h = 100, col = "blue", lty = 2)
dev.off()

# Auto-select soft threshold: lowest power with R^2 >= 0.80
fits <- data.frame(power = sft$fitIndices[, 1],
                   r2    = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
                   meank = sft$fitIndices[, 5])
above_80 <- fits$power[fits$r2 >= 0.80]
if (length(above_80) > 0) {
  selected_power <- min(above_80)
} else {
  selected_power <- 14  # fallback per WGCNA FAQ for n=20-30 signed
  cat("  WARNING: no power reached R^2 >= 0.80; using fallback power = 14\n")
}
cat("  Selected soft threshold power: ", selected_power, "\n", sep = "")

# Run blockwise modules
bwnet <- blockwiseModules(
  wgcna_expr,
  maxBlockSize    = 15000,
  TOMType         = "signed",
  power           = selected_power,
  networkType     = "signed",
  minModuleSize   = 30,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels    = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs        = FALSE,
  randomSeed      = wgcna_seed,
  verbose         = 0
)
saveRDS(bwnet, file.path(out_root, "WGCNA", "wgcna_blockwise.rds"))
cat("  Modules detected: ", length(unique(bwnet$colors)), "\n", sep = "")

# Module-trait correlation (group as binary trait)
module_eigengenes <- bwnet$MEs
group_design <- model.matrix(~ 0 + group, data = colData)
colnames(group_design) <- gsub("^group", "", colnames(group_design))

mod_trait_cor   <- cor(module_eigengenes, group_design, use = "p")
mod_trait_pval  <- corPvalueStudent(mod_trait_cor, nrow(wgcna_expr))

# Heatmap of module-trait
png(file.path(out_root, "WGCNA", "WGCNA_module_trait_heatmap.png"),
    width = 2400, height = 2000, res = 200)
text_mat <- paste0(signif(mod_trait_cor, 2), "\n(",
                   signif(mod_trait_pval, 2), ")")
dim(text_mat) <- dim(mod_trait_cor)
par(mar = c(10, 10, 4, 2))
labeledHeatmap(Matrix = mod_trait_cor,
               xLabels = colnames(mod_trait_cor),
               yLabels = colnames(mod_trait_cor),
               ySymbols = colnames(module_eigengenes),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = text_mat,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "WGCNA module-group correlation")
dev.off()

# Module-gene assignment table
module_df <- data.frame(
  gene_id = colnames(wgcna_expr),
  module_color  = labels2colors(bwnet$colors),
  module_number = bwnet$colors,
  stringsAsFactors = FALSE
)
module_df <- merge(module_df, anno, by.x = "gene_id", by.y = id_type,
                   all.x = TRUE)
write.csv(module_df,
          file.path(out_root, "WGCNA", "WGCNA_gene_module_assignments.csv"),
          row.names = FALSE)

cat("  WGCNA done\n")

# =====================================================================
# ===== 18. Save objects ==============================================
# =====================================================================
cat("\n===== Saving RDS objects =====\n")
saveRDS(dds,             "dds_full_8groups.rds")
saveRDS(vsd,             "vsd_full_8groups.rds")
saveRDS(deg_results,     "deg_results_full_8groups.rds")
saveRDS(anno,            "gene_annotation.rds")
saveRDS(contrast_tiers,  "contrast_tiers.rds")
saveRDS(go_results,      file.path(out_root, "GO_KEGG", "go_results.rds"))
saveRDS(kegg_results,    file.path(out_root, "GO_KEGG", "kegg_results.rds"))
saveRDS(limma_results,   file.path(out_root, "limma", "limma_results.rds"))
saveRDS(edger_results,   file.path(out_root, "edgeR", "edger_results.rds"))

cat("\n===== Pipeline complete =====\n")
cat("Output root: ", out_root, "/\n", sep = "")
cat("  QC/                 : Distribution plots, PCA, dispersion, sample distance\n")
cat("  DEG_results/        : 14 contrasts × {FULL, UP, DOWN}.csv\n")
cat("  Volcano/            : 14 EnhancedVolcano plots\n")
cat("  GO_KEGG/            : GO BP/CC/MF + KEGG (14 contrasts)\n")
cat("  Pathview/           : LFC overlay on KEGG pathways (default contrasts)\n")
cat("  gProfiler2/         : Manhattan plots + tables (default contrasts)\n")
cat("  limma/              : limma+voom cross-validation (default contrasts)\n")
cat("  edgeR/              : edgeR+glmTreat cross-validation (default contrasts)\n")
cat("  WGCNA/              : Soft threshold, modules, module-trait correlation\n")
cat("  GSEA_export/        : normalized_counts.gct + groups.cls\n")
cat("  tables/             : DEG summary + cross-platform overlap\n")
cat("\nRDS objects (working directory):\n")
cat("  dds_full_8groups.rds, vsd_full_8groups.rds, deg_results_full_8groups.rds\n")
cat("  gene_annotation.rds, contrast_tiers.rds\n")