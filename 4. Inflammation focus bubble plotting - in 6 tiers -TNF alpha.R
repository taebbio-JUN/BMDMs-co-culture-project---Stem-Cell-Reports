# =====================================================================
# Step 7: Figure 4 — 4-baseline integrated comparison
# Project: HN00273522 — BMDM × Macrophage × TNFα × Organoid co-culture
#
# Outputs:
#   Fig_4e_6way_Venn.{png,pdf}                     — 6-way circular Venn
#   Fig_4f_KEGG_bubble_matrix.{png,pdf,tiff}       — KEGG bubble matrix (14 contrasts)
#   Fig_4g_TopInflam_LFC_heatmap.{png,pdf,tiff}    — Top inflam LFC heatmap (14 contrasts)
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(RColorBrewer)
  library(clusterProfiler); library(enrichplot); library(DOSE)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

# =====================================================================
# ===== 1. Configuration ==============================================
# =====================================================================
out_root  <- "TNF-BMDMs-1k-5k-inflammation"
plot_font <- "Helvetica"

dir.create(out_root, showWarnings = FALSE)
dir.create(file.path(out_root, "tables"), showWarnings = FALSE)

padj_cutoff        <- 0.05
lfc_cutoff         <- 1.0
top_n_genes        <- 20
top_n_per_baseline <- 10

# =====================================================================
# ===== 2. Load Step 1 outputs ========================================
# =====================================================================
dds            <- readRDS("dds_full_8groups.rds")
deg_results    <- readRDS("deg_results_full_8groups.rds")
anno           <- readRDS("gene_annotation.rds")
contrast_tiers <- readRDS("contrast_tiers.rds")

cat("Loaded: ", nrow(dds), " genes × ", ncol(dds), " samples\n", sep = "")

# =====================================================================
# ===== 3. 4 baselines + contrast mapping =============================
# =====================================================================
baseline_map <- list(
  "vs Organoid (untreated)" = c(
    "T1_OrgBMDM1K_vs_Org",
    "T1_OrgBMDM5K_vs_Org",
    "T3_OrgTNFa_vs_Org",
    "T3_OrgTNFaBMDM1K_vs_Org",
    "T3_OrgTNFaBMDM5K_vs_Org"
  ),
  "vs Organoid + TNFα" = c(
    "T4_OrgTNFaBMDM1K_vs_OrgTNFa",
    "T4_OrgTNFaBMDM5K_vs_OrgTNFa"
  ),
  "vs Macrophage (untreated)" = c(
    "T2_OrgBMDM1K_vs_Mac",
    "T2_OrgBMDM5K_vs_Mac",
    "T5_MacTNFa_vs_Mac",
    "T5_OrgTNFaBMDM1K_vs_Mac",
    "T5_OrgTNFaBMDM5K_vs_Mac"
  ),
  "vs Macrophage + TNFα" = c(
    "T6_OrgTNFaBMDM1K_vs_MacTNFa",
    "T6_OrgTNFaBMDM5K_vs_MacTNFa"
  )
)
all_contrasts   <- unlist(baseline_map)
baseline_levels <- names(baseline_map)
cat("Total contrasts: ", length(all_contrasts), "\n", sep = "")

baseline_palette <- c(
  "vs Organoid (untreated)"   = "#1F4E79",   # dark blue
  "vs Organoid + TNFα"        = "#6FA8DC",   # light blue
  "vs Macrophage (untreated)" = "#A61C00",   # dark red
  "vs Macrophage + TNFα"      = "#E06666"    # light red
)

# Bubble plot contrast labels — 모두 "Org+BMDM ..." 통일
contrast_label_map <- c(
  T1_OrgBMDM1K_vs_Org         = "Org+BMDM 1K",
  T1_OrgBMDM5K_vs_Org         = "Org+BMDM 5K",
  T3_OrgTNFa_vs_Org           = "Org+TNFα",
  T3_OrgTNFaBMDM1K_vs_Org     = "Org+BMDM 1K+TNFα",
  T3_OrgTNFaBMDM5K_vs_Org     = "Org+BMDM 5K+TNFα",
  T4_OrgTNFaBMDM1K_vs_OrgTNFa = "Org+BMDM 1K+TNFα",
  T4_OrgTNFaBMDM5K_vs_OrgTNFa = "Org+BMDM 5K+TNFα",
  T2_OrgBMDM1K_vs_Mac         = "Org+BMDM 1K",
  T2_OrgBMDM5K_vs_Mac         = "Org+BMDM 5K",
  T5_MacTNFa_vs_Mac           = "Org+TNFα",
  T5_OrgTNFaBMDM1K_vs_Mac     = "Org+BMDM 1K+TNFα",
  T5_OrgTNFaBMDM5K_vs_Mac     = "Org+BMDM 5K+TNFα",
  T6_OrgTNFaBMDM1K_vs_MacTNFa = "Org+BMDM 1K+TNFα",
  T6_OrgTNFaBMDM5K_vs_MacTNFa = "Org+BMDM 5K+TNFα"
)

# =====================================================================
# ===== 4. Figure 4e — 6-way circular Venn diagram ===================
# =====================================================================
cat("\n===== Generating Fig 4e (6-way circular Venn) =====\n")

if (!requireNamespace("venn", quietly = TRUE)) install.packages("venn")
library(venn)

norm_counts <- counts(dds, normalized = TRUE)
group_vec   <- as.character(colData(dds)$group)

build_expressed_set <- function(group_name, threshold = 10, min_n = 2) {
  cols <- which(group_vec == group_name)
  if (length(cols) < 1) stop("No samples for group: ", group_name)
  expressed <- rowSums(norm_counts[, cols, drop = FALSE] >= threshold) >= min_n
  rownames(norm_counts)[expressed]
}

venn_groups <- c("Organoid","Macrophage","Organoid_TNFa","Macrophage_TNFa",
                 "Organoid_TNFa_BMDM_1K","Organoid_TNFa_BMDM_5K")
venn_labels <- c("Organoid","BMDM","Org+TNFα","BMDM+TNFα",
                 "Org+TNFα+BMDM 1K","Org+TNFα+BMDM 5K")

venn_sets <- list()
for (i in seq_along(venn_groups)) {
  venn_sets[[venn_labels[i]]] <- build_expressed_set(venn_groups[i])
  cat("  ", venn_groups[i], ": ", length(venn_sets[[venn_labels[i]]]),
      " expressed genes\n", sep = "")
}

# Save membership
all_g <- unique(unlist(venn_sets))
mem_df <- data.frame(gene_id = all_g, stringsAsFactors = FALSE)
for (lbl in venn_labels) mem_df[[lbl]] <- as.integer(all_g %in% venn_sets[[lbl]])
mem_df <- merge(mem_df, anno[, c("gene_id","SYMBOL","ENTREZID")],
                by = "gene_id", all.x = TRUE)
write.csv(mem_df,
          file.path(out_root, "tables", "Venn6_membership.csv"),
          row.names = FALSE)

venn_colors <- c("#1F4E79","#A61C00","#6FA8DC","#E06666","#33A02C","#FF7F00")

# PNG
png(file.path(out_root, "Fig_4e_6way_Venn.png"),
    width = 4000, height = 3500, res = 320, bg = "white")
par(mar = c(2, 2, 4, 2))
venn::venn(venn_sets, snames = venn_labels, ilabels = "counts",
           zcolor = venn_colors, opacity = 0.45, box = FALSE,
           ilcs = 1.0, sncs = 1.2, borders = TRUE, ellipse = TRUE)
title(main = "Expressed gene intersections — 6 key biological groups",
      cex.main = 1.8, font.main = 2, family = plot_font)
dev.off()

# PDF
pdf(file.path(out_root, "Fig_4e_6way_Venn.pdf"), width = 11, height = 9.5)
par(mar = c(2, 2, 4, 2))
venn::venn(venn_sets, snames = venn_labels, ilabels = "counts",
           zcolor = venn_colors, opacity = 0.45, box = FALSE,
           ilcs = 1.0, sncs = 1.2, borders = TRUE, ellipse = TRUE)
title(main = "Expressed gene intersections — 6 key biological groups",
      cex.main = 1.6, font.main = 2, family = plot_font)
dev.off()

# TIFF
tiff(file.path(out_root, "Fig_4e_6way_Venn.tiff"),
     width = 11, height = 9.5, units = "in",
     res = 300, bg = "white", compression = "lzw")
par(mar = c(2, 2, 4, 2))
venn::venn(venn_sets, snames = venn_labels, ilabels = "counts",
           zcolor = venn_colors, opacity = 0.45, box = FALSE,
           ilcs = 1.0, sncs = 1.2, borders = TRUE, ellipse = TRUE)
title(main = "Expressed gene intersections — 6 key biological groups",
      cex.main = 1.6, font.main = 2, family = plot_font)
dev.off()

cat("Fig 4e saved: PNG / PDF / TIFF\n")

# Intersection summary
cat("\n===== Venn intersection summary =====\n")
core_all6 <- Reduce(intersect, venn_sets)
cat("Core (in all 6 groups): ", length(core_all6), "\n", sep = "")
for (i in seq_along(venn_labels)) {
  uniq <- setdiff(venn_sets[[i]], unique(unlist(venn_sets[-i])))
  cat(sprintf("  %-22s unique : %d\n", venn_labels[i], length(uniq)))
}

# =====================================================================
# ===== 5. KEGG ORA for all 14 contrasts ==============================
# =====================================================================
cat("\n===== Running KEGG ORA for all 14 contrasts =====\n")

kegg_blacklist <- c(
  # Virus / Bacteria / Parasitic / Cancer / Cardiovascular / Neurodegenerative
  "Herpes simplex virus 1 infection","Epstein-Barr virus infection",
  "Kaposi sarcoma-associated herpesvirus infection",
  "Human papillomavirus infection","Human cytomegalovirus infection",
  "Human immunodeficiency virus 1 infection",
  "Human T-cell leukemia virus 1 infection",
  "Hepatitis B","Hepatitis C","Measles","Influenza A",
  "Coronavirus disease - COVID-19","Viral myocarditis",
  "Viral carcinogenesis","Viral life cycle - HIV-1",
  "Viral protein interaction with cytokine and cytokine receptor",
  "Staphylococcus aureus infection","Tuberculosis","Pertussis",
  "Legionellosis","Salmonella infection","Yersinia infection",
  "Vibrio cholerae infection","Pathogenic Escherichia coli infection",
  "Shigellosis","Bacterial invasion of epithelial cells",
  "Epithelial cell signaling in Helicobacter pylori infection",
  "Leishmaniasis","Malaria","Chagas disease",
  "African trypanosomiasis","Amoebiasis","Toxoplasmosis",
  "Pathways in cancer","Transcriptional misregulation in cancer",
  "MicroRNAs in cancer","Proteoglycans in cancer",
  "Chemical carcinogenesis - DNA adducts",
  "Chemical carcinogenesis - receptor activation",
  "Chemical carcinogenesis - reactive oxygen species",
  "Bladder cancer","Colorectal cancer","Pancreatic cancer",
  "Endometrial cancer","Glioma","Melanoma","Prostate cancer",
  "Thyroid cancer","Basal cell carcinoma","Renal cell carcinoma",
  "Small cell lung cancer","Non-small cell lung cancer",
  "Acute myeloid leukemia","Chronic myeloid leukemia",
  "Dilated cardiomyopathy","Hypertrophic cardiomyopathy",
  "Arrhythmogenic right ventricular cardiomyopathy",
  "Fluid shear stress and atherosclerosis","Lipid and atherosclerosis",
  "AGE-RAGE signaling pathway in diabetic complications",
  "Alzheimer disease","Parkinson disease","Huntington disease",
  "Amyotrophic lateral sclerosis","Prion disease",
  "Pathways of neurodegeneration - multiple diseases",
  "Rheumatoid arthritis","Systemic lupus erythematosus",
  "Type I diabetes mellitus","Autoimmune thyroid disease",
  "Graft-versus-host disease","Allograft rejection","Asthma",
  "Primary immunodeficiency",
  "Hematopoietic cell lineage",
  "Natural killer cell mediated cytotoxicity",
  "Neuroactive ligand-receptor interaction",
  "Maturity onset diabetes of the young",
  "Cytoskeleton in muscle cells"
)

filter_kegg <- function(res, blacklist) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(res)
  res_new <- res
  res_new@result <- res@result[!res@result$Description %in% blacklist, ]
  res_new
}

collapse_to_gene <- function(deg_df) {
  deg_df %>%
    dplyr::filter(!is.na(ENTREZID), !is.na(padj)) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::group_by(ENTREZID) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    as.data.frame()
}
deg_gene <- lapply(deg_results, collapse_to_gene)
universe_entrez <- unique(na.omit(as.character(anno$ENTREZID)))

run_kegg_ora <- function(deg_df) {
  sig <- deg_df %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::pull(ENTREZID) %>% na.omit() %>% unique() %>% as.character()
  if (length(sig) < 10) return(NULL)
  res <- tryCatch(
    enrichKEGG(gene = sig, organism = "mmu",
               universe = universe_entrez,
               pvalueCutoff = 0.1, qvalueCutoff = 0.25),
    error = function(e) NULL
  )
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    res <- setReadable(res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    res <- filter_kegg(res, kegg_blacklist)
  }
  res
}

kegg_list <- list()
for (cn in all_contrasts) {
  cat("  ", cn, "\n", sep = "")
  kegg_list[[cn]] <- run_kegg_ora(deg_gene[[cn]])
}

# =====================================================================
# ===== 6. Figure 4f — KEGG bubble matrix (data prep) ================
# =====================================================================
cat("\n===== Preparing Fig 4f data =====\n")

contrast_baseline_map <- list()
for (bl in names(baseline_map)) {
  for (cn in baseline_map[[bl]]) contrast_baseline_map[[cn]] <- bl
}

build_kegg_long_bubble <- function(kegg_list, contrast_baseline_map) {
  rows <- list()
  for (cn in names(contrast_baseline_map)) {
    r <- kegg_list[[cn]]
    if (is.null(r) || nrow(as.data.frame(r)) == 0) next
    df <- as.data.frame(r)
    df$Contrast <- cn
    df$Baseline <- contrast_baseline_map[[cn]]
    rows[[cn]] <- df[, c("ID","Description","Count","pvalue","p.adjust",
                         "Contrast","Baseline")]
  }
  do.call(rbind, rows)
}

kegg_bubble_long <- build_kegg_long_bubble(kegg_list, contrast_baseline_map)
kegg_bubble_long$Baseline <- factor(kegg_bubble_long$Baseline,
                                    levels = baseline_levels)
kegg_bubble_long$logp <- -log10(kegg_bubble_long$pvalue)

selected_pathways <- kegg_bubble_long %>%
  dplyr::group_by(Baseline, Description) %>%
  dplyr::summarise(best_logp = max(logp), .groups = "drop") %>%
  dplyr::group_by(Baseline) %>%
  dplyr::arrange(dplyr::desc(best_logp), .by_group = TRUE) %>%
  dplyr::slice_head(n = top_n_per_baseline) %>%
  dplyr::ungroup()

pathway_primary <- selected_pathways %>%
  dplyr::group_by(Description) %>%
  dplyr::slice_max(best_logp, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Baseline, dplyr::desc(best_logp))

kegg_bubble_sub <- kegg_bubble_long %>%
  dplyr::filter(Description %in% pathway_primary$Description)
kegg_bubble_sub$Description <- factor(kegg_bubble_sub$Description,
                                      levels = rev(pathway_primary$Description))

contrast_order <- unlist(baseline_map, use.names = FALSE)
kegg_bubble_sub$Contrast <- factor(kegg_bubble_sub$Contrast,
                                   levels = contrast_order)

# =====================================================================
# ===== Fig 4f — Bubble plot 빌드 + 저장 ==============================
# =====================================================================

# ----- 4f-A. 폰트 / 크기 (튜닝 포인트) -------------------------------
font_sz <- list(
  base       = 18, title      = 24, subtitle   = 16,
  axis_x     = 16, axis_y     = 16,
  legend_ttl = 16, legend_txt = 14
)
bubble_size_range <- c(4, 14)
bubble_stroke     <- 0.5
strip_size        <- 14
strip_margin      <- 4

# ----- 4f-B. Baseline 범례 dummy -------------------------------------
baseline_legend_data <- data.frame(
  Baseline = factor(baseline_levels, levels = baseline_levels)
)

# ----- 4f-C. 플롯 빌드 -----------------------------------------------
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
if (!has_ggh4x) {
  install.packages("ggh4x")
  has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
}

p_4f <- ggplot(kegg_bubble_sub,
               aes(x = Contrast, y = Description)) +
  geom_point(aes(size = Count, fill = logp),
             shape = 21, color = "black",
             stroke = bubble_stroke, alpha = 0.95) +
  geom_point(data = baseline_legend_data,
             aes(color = Baseline),
             x = NA, y = NA, size = 0, na.rm = TRUE,
             inherit.aes = FALSE, show.legend = TRUE)

if (has_ggh4x) {
  strip_fills <- unname(baseline_palette[baseline_levels])
  p_4f <- p_4f +
    ggh4x::facet_grid2(
      . ~ Baseline,
      scales = "free_x", space = "free_x",
      strip = ggh4x::strip_themed(
        background_x = lapply(strip_fills, function(col)
          element_rect(fill = col, color = "black", linewidth = 0.6)),
        text_x = lapply(strip_fills, function(col)
          element_text(color  = NA,
                       size   = strip_size,
                       margin = margin(t = strip_margin,
                                       b = strip_margin)))
      )
    )
} else {
  p_4f <- p_4f +
    facet_grid(. ~ Baseline, scales = "free_x", space = "free_x")
}

p_4f <- p_4f +
  scale_fill_gradient2(
    low = "#5865f2", mid = "#faffff", high = "#A61C00",
    midpoint = 5,
    name   = expression(-log[10]~italic(P)),
    limits = c(0, NA),
    guide  = guide_colorbar(barwidth = 1.4, barheight = 12, order = 1)
  ) +
  scale_size_continuous(
    range  = bubble_size_range, name = "Gene count",
    breaks = c(30, 60, 90), guide = guide_legend(order = 2)
  ) +
  scale_color_manual(
    values = baseline_palette, name = "Baseline", drop = FALSE,
    guide  = guide_legend(
      override.aes = list(size = 7, shape = 15, alpha = 1), order = 3
    )
  ) +
  scale_x_discrete(labels = contrast_label_map) +
  labs(x = NULL, y = NULL,
       title = "KEGG pathway enrichment — bubble matrix (4-baseline × 14 contrasts)",
       subtitle = paste0("Top ", top_n_per_baseline,
                         " pathways per baseline (union); bubble size = gene count, color = -log10(p)")) +
  theme_bw(base_size = font_sz$base, base_family = plot_font) +
  theme(
    text             = element_text(family = plot_font),
    plot.title       = element_text(face = "bold", hjust = 0.5,
                                    size = font_sz$title, margin = margin(b = 6)),
    plot.subtitle    = element_text(hjust = 0.5, size = font_sz$subtitle,
                                    color = "grey30", margin = margin(b = 14)),
    axis.text.x      = element_text(size = font_sz$axis_x, face = "bold",
                                    color = "black", angle = 45,
                                    hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = font_sz$axis_y, color = "black"),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.spacing.x  = unit(0.5, "lines"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = font_sz$legend_ttl),
    legend.text      = element_text(size = font_sz$legend_txt),
    legend.box       = "vertical",
    legend.spacing.y = unit(0.4, "cm"),
    plot.margin      = margin(15, 15, 15, 15)
  )

# ----- 4f-D. 화면 출력 -----------------------------------------------
while (!is.null(dev.list())) dev.off()
print(p_4f)

# ----- 4f-E. 저장 (PNG / PDF / TIFF) ---------------------------------
fig_w_4f <- 20     # ★ 직접 입력 (inches)
fig_h_4f <- 16     # ★ 직접 입력 (inches)

ggsave(file.path(out_root, "Fig_4f_KEGG_bubble_matrix.png"),
       p_4f, width = fig_w_4f, height = fig_h_4f,
       dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_root, "Fig_4f_KEGG_bubble_matrix.pdf"),
       p_4f, width = fig_w_4f, height = fig_h_4f,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(out_root, "Fig_4f_KEGG_bubble_matrix.tiff"),
       p_4f, width = fig_w_4f, height = fig_h_4f,
       device = "tiff", dpi = 300, bg = "white",
       compression = "lzw", limitsize = FALSE)

cat("\nFig 4f saved: PNG / PDF / TIFF (",
    fig_w_4f, " × ", fig_h_4f, " in)\n", sep = "")

write.csv(kegg_bubble_sub,
          file.path(out_root, "tables", "KEGG_bubble_matrix_long.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 7. Figure 4g — Top inflammation gene LFC heatmap =============
# =====================================================================
cat("\n===== Preparing Fig 4g data =====\n")

# ----- 7-A. Inflammation panel ---------------------------------------
inflam_panel_expanded <- c(
  # Cytokines
  "Tnf","Il1a","Il1b","Il6","Il12a","Il12b","Il18","Il33","Il17a","Il17f",
  "Il10","Tgfb1","Tgfb2","Tgfb3","Ifna1","Ifnb1","Ifng",
  "Csf1","Csf2","Csf3","Osm","Lif","Clcf1",
  # Cytokine receptors
  "Tnfrsf1a","Tnfrsf1b","Il6st","Il6ra","Il1r1","Il1r2",
  "Il10ra","Il10rb","Il12rb1","Il12rb2","Il18r1","Il18rap","Il23r",
  "Csf1r","Csf2ra","Csf2rb","Csf3r","Ifnar1","Ifnar2","Ifngr1","Ifngr2",
  "Stat1","Stat2","Stat3","Stat4","Stat5a","Stat5b","Stat6",
  "Jak1","Jak2","Jak3","Tyk2",
  # Chemokines
  "Ccl2","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9",
  "Ccl12","Ccl17","Ccl19","Ccl20","Ccl22","Ccl24","Ccl25",
  "Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcl10","Cxcl11",
  "Cxcl12","Cxcl13","Cxcl16",
  "Ccr1","Ccr2","Ccr5","Ccr6","Ccr7","Cxcr2","Cxcr3","Cxcr4","Cmklr1",
  # Acute phase / alarmins
  "Saa3","Crp","Hp","Lbp","Lcn2","Ptx3","Cp","Orm1","Orm2",
  "S100a4","S100a6","S100a8","S100a9","Hmgb1",
  # Vascular / adhesion
  "Icam1","Vcam1","Sele","Selp","Ptgs2",
  # TLR / NLR
  "Tlr2","Tlr3","Tlr4","Tlr7","Tlr9","Myd88","Ticam1","Tirap",
  "Nlrp3","Aim2","Casp1","Casp4","Pycard","Gsdmd","Gsdme",
  # NF-κB
  "Nfkb1","Nfkb2","Rela","Relb","Nfkbia","Nfkbid","Nfkbie",
  "Ikbkg","Tnfaip3","Tnip1","Tnip3","Traf1","Traf3",
  # IFN / ISG
  "Irf1","Irf3","Irf5","Irf7","Irf8","Irf9",
  "Mx1","Mx2","Isg15","Rsad2","Oas1a","Oas2","Oas3",
  "Ifit1","Ifit2","Ifit3","Usp18","Bst2","Gbp2","Gbp3","Gbp4",
  # SOCS
  "Socs1","Socs2","Socs3","Cish",
  # Macrophage polarization
  "Nos2","Arg1","Cd86","Cd80","Cd163","Cd68","Adgre1","Mrc1","Mafb",
  "Marco","Msr1","Trem2",
  # Complement
  "C3","C1qa","C1qb","C1qc","C3ar1","C5ar1",
  # ECM
  "Mmp9","Mmp12","Mmp13","Mmp3","Mmp8",
  # Apoptotic
  "Fas","Fasl","Casp3","Casp8","Bcl2a1a","Bcl2a1b","Bcl2a1d",
  # Anti-inflammatory feedback
  "Il1rn","Hmox1","Klf2","Klf4"
)
inflam_panel_expanded <- unique(inflam_panel_expanded)
cat("Inflammation panel: ", length(inflam_panel_expanded), " genes\n", sep = "")

# ----- 7-B. Collect LFC, rank, select top ----------------------------
lfc_per_gene <- list()
for (cn in all_contrasts) {
  d <- deg_results[[cn]] %>%
    dplyr::filter(SYMBOL %in% inflam_panel_expanded, !is.na(padj))
  lfc_per_gene[[cn]] <- data.frame(
    SYMBOL = d$SYMBOL, LFC = d$log2FoldChange,
    padj = d$padj, Contrast = cn,
    stringsAsFactors = FALSE
  )
}
lfc_long <- do.call(rbind, lfc_per_gene)
lfc_long <- lfc_long %>% dplyr::distinct(SYMBOL, Contrast, .keep_all = TRUE)

gene_rank <- lfc_long %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    max_abs_lfc = max(abs(LFC), na.rm = TRUE),
    min_padj    = min(padj, na.rm = TRUE),
    n_sig       = sum(padj < padj_cutoff & abs(LFC) > lfc_cutoff, na.rm = TRUE),
    rank_score  = max_abs_lfc * (-log10(pmax(min_padj, 1e-300))),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(rank_score))

top_genes <- head(gene_rank$SYMBOL, top_n_genes)
cat("\nTop ", top_n_genes, " inflammation genes selected\n", sep = "")
print(head(gene_rank %>% dplyr::filter(SYMBOL %in% top_genes), 15))

# ----- 7-C. Build LFC / padj matrices --------------------------------
ordered_contrasts <- unlist(baseline_map, use.names = FALSE)

lfc_mat <- matrix(NA, nrow = length(top_genes),
                  ncol = length(ordered_contrasts),
                  dimnames = list(top_genes, ordered_contrasts))
padj_mat <- lfc_mat
for (i in seq_along(top_genes)) {
  for (j in seq_along(ordered_contrasts)) {
    d <- deg_results[[ordered_contrasts[j]]] %>%
      dplyr::filter(SYMBOL == top_genes[i])
    if (nrow(d) > 0) {
      lfc_mat[i, j]  <- d$log2FoldChange[1]
      padj_mat[i, j] <- d$padj[1]
    }
  }
}
keep_rows <- rowSums(!is.na(lfc_mat)) > 0
lfc_mat  <- lfc_mat[keep_rows, , drop = FALSE]
padj_mat <- padj_mat[keep_rows, , drop = FALSE]

# ----- 7-D. Column labels (vs ... 제거, 중복 방지용 공백) ------------
col_labels <- c(
  T1_OrgBMDM1K_vs_Org         = "Org+BMDM 1K",
  T1_OrgBMDM5K_vs_Org         = "Org+BMDM 5K",
  T3_OrgTNFa_vs_Org           = "Org+TNFα",
  T3_OrgTNFaBMDM1K_vs_Org     = "Org+BMDM 1K+TNFα",
  T3_OrgTNFaBMDM5K_vs_Org     = "Org+BMDM 5K+TNFα",
  T4_OrgTNFaBMDM1K_vs_OrgTNFa = "Org+BMDM 1K+TNFα ",
  T4_OrgTNFaBMDM5K_vs_OrgTNFa = "Org+BMDM 5K+TNFα ",
  T2_OrgBMDM1K_vs_Mac         = "Org+BMDM 1K ",
  T2_OrgBMDM5K_vs_Mac         = "Org+BMDM 5K ",
  T5_MacTNFa_vs_Mac           = "Org+TNFα ",
  T5_OrgTNFaBMDM1K_vs_Mac     = "Org+BMDM 1K+TNFα  ",
  T5_OrgTNFaBMDM5K_vs_Mac     = "Org+BMDM 5K+TNFα  ",
  T6_OrgTNFaBMDM1K_vs_MacTNFa = "Org+BMDM 1K+TNFα   ",
  T6_OrgTNFaBMDM5K_vs_MacTNFa = "Org+BMDM 5K+TNFα   "
)
colnames(lfc_mat)  <- col_labels[ordered_contrasts]
colnames(padj_mat) <- col_labels[ordered_contrasts]

# ----- 7-E. Column annotation ----------------------------------------
col_baselines <- c(
  rep("vs Organoid (untreated)",   5),
  rep("vs Organoid + TNFα",        2),
  rep("vs Macrophage (untreated)", 5),
  rep("vs Macrophage + TNFα",      2)
)
col_group <- data.frame(
  Baseline = factor(col_baselines, levels = baseline_levels),
  row.names = colnames(lfc_mat)
)
ann_color <- list(Baseline = baseline_palette)

lfc_max <- max(abs(lfc_mat), na.rm = TRUE)
lfc_max <- min(lfc_max, 12)

# =====================================================================
# ===== Fig 4g — ComplexHeatmap 빌드 + 저장 ===========================
# =====================================================================

# ----- 4g-A. 폰트 / 셀 크기 (튜닝 포인트) ----------------------------
hm_params_4g <- list(
  fontsize        = 18,    # ★ 전체 base (title, legend)
  fontsize_row    = 14,    # ★ gene symbol
  fontsize_col    = 12,    # ★ column label
  cellwidth       = 32,    # ★ 셀 너비
  cellheight      = 24,    # ★ 셀 높이
  legend_title    = 14     # ★ Log2(FC) 라벨 크기
)

# ----- 4g-B. ComplexHeatmap 빌드 -------------------------------------
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
}
library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun_4g <- colorRamp2(
  c(-lfc_max, -lfc_max/2, 0, lfc_max/2, lfc_max),
  c("#5865f2", "#a8b0f5", "#e8e8e8", "#F0A0A0", "#F05650")
)

top_anno_4g <- HeatmapAnnotation(
  Baseline = col_group$Baseline,
  col      = list(Baseline = ann_color$Baseline),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Baseline = list(
      title     = "Baseline",
      title_gp  = gpar(fontsize = hm_params_4g$fontsize, fontface = "bold"),
      labels_gp = gpar(fontsize = hm_params_4g$fontsize_col)
    )
  )
)

ht_4g <- Heatmap(
  lfc_mat,
  name              = "Log2(FC)",
  col               = col_fun_4g,
  cluster_rows      = FALSE, cluster_columns = FALSE,
  show_row_names    = TRUE, show_column_names = TRUE,
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = hm_params_4g$fontsize_row),
  column_names_gp   = gpar(fontsize = hm_params_4g$fontsize_col),
  column_names_rot  = 45,
  top_annotation    = top_anno_4g,
  column_split      = col_group$Baseline,
  column_title      = NULL,
  border            = TRUE,
  rect_gp           = gpar(col = "black", lwd = 0.5),
  width             = unit(hm_params_4g$cellwidth  * ncol(lfc_mat), "pt"),
  height            = unit(hm_params_4g$cellheight * nrow(lfc_mat), "pt"),
  heatmap_legend_param = list(
    title          = "Log2(FC)",
    title_gp       = gpar(fontsize = hm_params_4g$legend_title, fontface = "bold"),
    labels_gp      = gpar(fontsize = hm_params_4g$fontsize_col),
    title_position = "leftcenter-rot",
    legend_height  = unit(5, "cm"),
    grid_width     = unit(0.5, "cm")
  )
)

# ----- 4g-C. draw 함수 ----------------------------------------------
draw_ht_4g <- function() {
  draw(ht_4g,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       column_title    = paste0("Top ", top_n_genes,
                                " inflammation genes — LFC (4-baseline comprehensive)"),
       column_title_gp = gpar(fontsize = hm_params_4g$fontsize + 2,
                              fontface = "bold"))
}

# ----- 4g-D. 화면 출력 -----------------------------------------------
while (!is.null(dev.list())) dev.off()
draw_ht_4g()

# ----- 4g-E. 저장 (PNG / PDF / TIFF) ---------------------------------
fig_w_4g <- 18     # ★ 직접 입력 (inches)
fig_h_4g <- 12     # ★ 직접 입력 (inches)

png(file.path(out_root, "Fig_4g_TopInflam_LFC_heatmap.png"),
    width = fig_w_4g, height = fig_h_4g,
    units = "in", res = 300, bg = "white")
draw_ht_4g()
dev.off()

pdf(file.path(out_root, "Fig_4g_TopInflam_LFC_heatmap.pdf"),
    width = fig_w_4g, height = fig_h_4g, bg = "white")
draw_ht_4g()
dev.off()

tiff(file.path(out_root, "Fig_4g_TopInflam_LFC_heatmap.tiff"),
     width = fig_w_4g, height = fig_h_4g,
     units = "in", res = 300, bg = "white",
     compression = "lzw")
draw_ht_4g()
dev.off()

cat("\nFig 4g saved: PNG / PDF / TIFF (",
    fig_w_4g, " × ", fig_h_4g, " in)\n", sep = "")

# ----- 4g-F. Tables export -------------------------------------------
lfc_df_out  <- data.frame(SYMBOL = rownames(lfc_mat),  lfc_mat,
                          check.names = FALSE)
padj_df_out <- data.frame(SYMBOL = rownames(padj_mat), padj_mat,
                          check.names = FALSE)
write.csv(lfc_df_out,
          file.path(out_root, "tables", "Top_inflam_LFC_AllBaselines.csv"),
          row.names = FALSE)
write.csv(padj_df_out,
          file.path(out_root, "tables", "Top_inflam_padj_AllBaselines.csv"),
          row.names = FALSE)
write.csv(gene_rank,
          file.path(out_root, "tables", "Inflam_gene_ranking_AllBaselines.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 8. Summary ====================================================
# =====================================================================
cat("\n===== Step 7 complete =====\n")
cat("Output root: ", out_root, "/\n\n", sep = "")
cat("  Fig_4e_6way_Venn.{png,pdf,tiff}\n")
cat("  Fig_4f_KEGG_bubble_matrix.{png,pdf,tiff}\n")
cat("  Fig_4g_TopInflam_LFC_heatmap.{png,pdf,tiff}\n")
cat("\nTables/:\n")
cat("  Venn6_membership.csv\n")
cat("  KEGG_bubble_matrix_long.csv\n")
cat("  Inflam_gene_ranking_AllBaselines.csv\n")
cat("  Top_inflam_LFC_AllBaselines.csv\n")
cat("  Top_inflam_padj_AllBaselines.csv\n")

# =====================================================================
# Fig 4g — Top inflammation gene LFC heatmap (print 가능)
# Fig 1g와 동일한 패턴
# =====================================================================

# ===== 4g-3. 폰트 / 셀 크기 — 튜닝 포인트 ============================
hm_params_4g <- list(
  fontsize             = 18,    # 전체 base
  fontsize_row         = 16,    # gene symbol
  fontsize_col         = 14,    # column label
  cellwidth            = 32,    # 셀 너비
  cellheight           = 32,    # 셀 높이
  legend_title         = 14,    # Log2(FC) 라벨
  anno_legend_title    = 20,    # ★ 추가 — "Baseline" 제목
  anno_legend_labels   = 16     # ★ 추가 — "vs Organoid ..." 라벨
)

# ===== 4g-4. Heatmap with ComplexHeatmap ============================
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
}
library(ComplexHeatmap)
library(circlize)
library(grid)

# Color function
col_fun_4g <- colorRamp2(
  c(-lfc_max, -lfc_max/2, 0, lfc_max/2, lfc_max),
  c("#5865f2", "#a8b0f5", "#e8e8e8", "#F0A0A0", "#F05650")
)

# Top annotation (Baseline 색 막대)
top_anno_4g <- HeatmapAnnotation(
  Baseline = col_group$Baseline,
  col      = list(Baseline = ann_color$Baseline),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Baseline = list(
      title     = "Baseline",
      title_gp  = gpar(fontsize = hm_params_4g$anno_legend_title,   # ★ 변경
                       fontface = "bold"),
      labels_gp = gpar(fontsize = hm_params_4g$anno_legend_labels,  # ★ 변경
                       fontface = "bold")
    )
  )
)

# Heatmap 빌드
ht_4g <- Heatmap(
  lfc_mat,
  name              = "Log2(FC)",
  col               = col_fun_4g,
  cluster_rows      = FALSE, cluster_columns = FALSE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = hm_params_4g$fontsize_row,
                           fontface = "bold"),
  column_names_gp   = gpar(fontsize = hm_params_4g$fontsize_col,
                           fontface = "bold"),
  column_names_rot  = 45,
  top_annotation    = top_anno_4g,
  column_split      = col_group$Baseline,
  column_title      = NULL,
  border            = TRUE,
  border_gp         = gpar(col = "black", lwd = 1.5),
  rect_gp           = gpar(col = "black", lwd = 1.2),
  width             = unit(hm_params_4g$cellwidth  * ncol(lfc_mat), "pt"),
  height            = unit(hm_params_4g$cellheight * nrow(lfc_mat), "pt"),
  heatmap_legend_param = list(
    title          = "Log2(FC)",
    title_gp       = gpar(fontsize = hm_params_4g$legend_title,
                          fontface = "bold"),
    labels_gp      = gpar(fontsize = hm_params_4g$fontsize_col,
                          fontface = "bold"),
    title_position = "leftcenter-rot",
    legend_height  = unit(5, "cm"),
    grid_width     = unit(0.5, "cm")
  )
)

# ===== 4g-5. draw 함수 (저장 시 재사용) ==============================
draw_ht_4g <- function() {
  draw(ht_4g,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       column_title    = paste0("Top ", top_n_genes,
                                " inflammation genes — LFC (4-baseline comprehensive)"),
       column_title_gp = gpar(fontsize = hm_params_4g$fontsize + 2,
                              fontface = "bold"))
}

# ===== 4g-6. 화면 출력 ===============================================
while (!is.null(dev.list())) dev.off()
draw_ht_4g()

# ===== 4g-7. 저장 (PNG / PDF / TIFF) =================================
fig_h_4g <- 12     # ★ 직접 입력 (inches)
fig_w_4g <- 18     # ★ 직접 입력 (inches) — 14 contrasts라 길게

# PNG
png(file.path(out_root, "Fig_4g_TopInflam_LFC_heatmap.png"),
    width = fig_w_4g, height = fig_h_4g,
    units = "in", res = 300, bg = "white")
draw_ht_4g()
dev.off()

# PDF (vector)
pdf(file.path(out_root, "Fig_4g_TopInflam_LFC_heatmap.pdf"),
    width = fig_w_4g, height = fig_h_4g, bg = "white")
draw_ht_4g()
dev.off()

# TIFF (출판용, 300 dpi + LZW)
tiff(file.path(out_root, "Fig_4g_TopInflam_LFC_heatmap.tiff"),
     width = fig_w_4g, height = fig_h_4g,
     units = "in", res = 300, bg = "white",
     compression = "lzw")
draw_ht_4g()
dev.off()

cat("\nFig 4g saved: PNG / PDF / TIFF (",
    fig_w_4g, " × ", fig_h_4g, " in)\n", sep = "")

# ===== 4g-8. Tables export ===========================================
lfc_df_out  <- data.frame(SYMBOL = rownames(lfc_mat),  lfc_mat,
                          check.names = FALSE)
padj_df_out <- data.frame(SYMBOL = rownames(padj_mat), padj_mat,
                          check.names = FALSE)
write.csv(lfc_df_out,
          file.path(out_root, "tables", "Top_inflam_LFC_AllBaselines.csv"),
          row.names = FALSE)
write.csv(padj_df_out,
          file.path(out_root, "tables", "Top_inflam_padj_AllBaselines.csv"),
          row.names = FALSE)
write.csv(gene_rank,
          file.path(out_root, "tables", "Inflam_gene_ranking_AllBaselines.csv"),
          row.names = FALSE)

# =====================================================================
# Fig 4f — KEGG bubble matrix (4-baseline × 14 contrasts)
# Fig 1g와 동일한 저장 패턴
# =====================================================================

# ===== 4f-3. 폰트 / 크기 — 튜닝 포인트 ===============================
font_sz <- list(
  title       = 24,   # ★ 메인 타이틀
  subtitle    = 16,   # ★ 서브타이틀
  axis_x      = 16,   # ★ 하단 contrast 라벨
  axis_y      = 22,   # ★ 좌측 pathway 이름
  legend_ttl  = 16,   # ★ 범례 제목
  legend_txt  = 16,   # ★ 범례 값
  base        = 18    # ★ theme_bw base
)
bubble_size_range <- c(4, 14)
bubble_stroke     <- 0.5
strip_size        <- 14    # ★ 상단 색 막대 두께 (글자는 투명)
strip_margin      <- 4

# ===== 4f-4. Bubble plot 빌드 ========================================
baseline_legend_data <- data.frame(
  Baseline = factor(baseline_levels, levels = baseline_levels)
)
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)

p_4f <- ggplot(kegg_bubble_sub,
               aes(x = Contrast, y = Description)) +
  geom_point(aes(size = Count, fill = logp),
             shape = 21, color = "black",
             stroke = bubble_stroke, alpha = 0.95) +
  geom_point(data = baseline_legend_data,
             aes(color = Baseline),
             x = NA, y = NA, size = 0, na.rm = TRUE,
             inherit.aes = FALSE, show.legend = TRUE)

if (has_ggh4x) {
  strip_fills <- unname(baseline_palette[baseline_levels])
  
  p_4f <- p_4f +
    ggh4x::facet_grid2(
      . ~ Baseline,
      scales = "free_x", space = "free_x",
      strip = ggh4x::strip_themed(
        background_x = lapply(strip_fills, function(col)
          element_rect(fill = col, color = "black", linewidth = 0.6)),
        text_x = lapply(strip_fills, function(col)
          element_text(color  = NA,
                       size   = strip_size,
                       margin = margin(t = strip_margin,
                                       b = strip_margin)))
      )
    )
} else {
  p_4f <- p_4f +
    facet_grid(. ~ Baseline, scales = "free_x", space = "free_x")
}

p_4f <- p_4f +
  scale_fill_gradient2(
    low = "#5865f2", mid = "#faffff", high = "#A61C00",
    midpoint = 5,
    name   = expression(-log[10]~italic(P)),
    limits = c(0, NA),
    guide  = guide_colorbar(barwidth = 1.4, barheight = 12, order = 1)
  ) +
  scale_size_continuous(
    range  = bubble_size_range,
    name   = "Gene count",
    breaks = c(30, 60, 90),
    guide  = guide_legend(order = 2)
  ) +
  scale_color_manual(
    values = baseline_palette,
    name   = "Baseline",
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(size = 7, shape = 15, alpha = 1),
      order = 3
    )
  ) +
  scale_x_discrete(labels = contrast_label_map) +
  labs(x = NULL, y = NULL,
       title = "KEGG pathway enrichment — bubble matrix (4-baseline × 14 contrasts)",
       subtitle = paste0("Top ", top_n_per_baseline,
                         " pathways per baseline (union); bubble size = gene count, color = -log10(p)")) +
  theme_bw(base_size = font_sz$base, base_family = plot_font) +
  theme(
    text             = element_text(family = plot_font),
    plot.title       = element_text(face = "bold", hjust = 0.5,
                                    size = font_sz$title,
                                    margin = margin(b = 6)),
    plot.subtitle    = element_text(hjust = 0.5,
                                    size = font_sz$subtitle,
                                    color = "grey30",
                                    margin = margin(b = 14)),
    axis.text.x      = element_text(size = font_sz$axis_x,
                                    face = "bold", color = "black",
                                    angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = font_sz$axis_y, color = "black"),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.spacing.x  = unit(0.5, "lines"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = font_sz$legend_ttl),
    legend.text      = element_text(size = font_sz$legend_txt),
    legend.box       = "vertical",
    legend.spacing.y = unit(0.4, "cm"),
    plot.margin      = margin(15, 15, 15, 15)
  )

# ===== 4f-5. 화면 출력 ===============================================
while (!is.null(dev.list())) dev.off()
print(p_4f)

# ===== 4f-6. 저장 (PNG / PDF / TIFF) =================================
fig_h_4f <- 16     # ★ 직접 입력 (inches)
fig_w_4f <- 20     # ★ 직접 입력 (inches)

# PNG
ggsave(file.path(out_root, "Fig_4f_KEGG_bubble_matrix.png"),
       p_4f, width = fig_w_4f, height = fig_h_4f,
       dpi = 300, bg = "white", limitsize = FALSE)

# PDF (vector)
ggsave(file.path(out_root, "Fig_4f_KEGG_bubble_matrix.pdf"),
       p_4f, width = fig_w_4f, height = fig_h_4f,
       device = cairo_pdf, limitsize = FALSE)

# TIFF (출판용, 300 dpi + LZW)
ggsave(file.path(out_root, "Fig_4f_KEGG_bubble_matrix.tiff"),
       p_4f, width = fig_w_4f, height = fig_h_4f,
       device = "tiff", dpi = 300, bg = "white",
       compression = "lzw", limitsize = FALSE)

cat("\nFig 4f saved: PNG / PDF / TIFF\n")

# ===== 4f-7. Tables export ===========================================
write.csv(kegg_bubble_sub,
          file.path(out_root, "tables", "KEGG_bubble_matrix_long.csv"),
          row.names = FALSE)