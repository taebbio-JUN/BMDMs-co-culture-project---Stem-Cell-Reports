# =====================================================================
# Step 8: Figure 5 — Stem cell / lineage analysis (4-baseline framework)
# Project: HN00273522 — BMDM x Macrophage x TNFa x Organoid co-culture
#
# Gene panel: 8-panel stem cell + lineage (fate decision focus)
#   ISC_markers, Wnt_signaling, Notch_signaling, Hippo_YAP,
#   Enterocyte, Goblet, Paneth, Enteroendocrine
#
# Enrichment: KEGG + GO (BP, MF, CC)
#
# Two layouts:
#   Layout A (basal, no TNFa) — 4 contrasts, 2 panel (vs Org / vs Mac)
#     Main figure
#   Layout B (full 14 contrasts) — 4 baseline panels
#     Supplementary figure
#
# Outputs:
#   Fig_5e_6way_Venn_SC.{png,pdf}
#   Layout A:
#     Fig_5f_KEGG_bubble_basal.{png,pdf}
#     Fig_5g_GO_BP_bubble_basal.{png,pdf}
#     Fig_5h_GO_MF_bubble_basal.{png,pdf}
#     Fig_5i_GO_CC_bubble_basal.{png,pdf}
#     Fig_5j_StemCell_LFC_heatmap_basal.{png,pdf}
#   Layout B:
#     Fig_S5f_KEGG_bubble_full.{png,pdf}
#     Fig_S5g_GO_BP_bubble_full.{png,pdf}
#     Fig_S5h_GO_MF_bubble_full.{png,pdf}
#     Fig_S5i_GO_CC_bubble_full.{png,pdf}
#     Fig_S5j_StemCell_LFC_heatmap_full.{png,pdf}
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
  library(clusterProfiler); library(enrichplot); library(DOSE)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

# =====================================================================
# ===== 1. Configuration =============================================
# =====================================================================
out_root  <- "Step10_StemCell_4baseline"
plot_font <- "Helvetica"

dir.create(out_root, showWarnings = FALSE)
dir.create(file.path(out_root, "tables"), showWarnings = FALSE)

padj_cutoff        <- 0.05
lfc_cutoff         <- 1.0
top_n_per_panel    <- 5         # heatmap: 5 genes per gene-panel x 8 = 40 top
top_n_per_baseline <- 10        # bubble plot pathways per baseline (Layout B)
top_n_basal        <- 12        # bubble plot pathways per baseline (Layout A)

# =====================================================================
# ===== 2. Load Step 1 outputs =======================================
# =====================================================================
dds            <- readRDS("dds_full_8groups.rds")
deg_results    <- readRDS("deg_results_full_8groups.rds")
anno           <- readRDS("gene_annotation.rds")
contrast_tiers <- readRDS("contrast_tiers.rds")

cat("Loaded: ", nrow(dds), " genes x ", ncol(dds), " samples\n", sep = "")

# =====================================================================
# ===== 3. Define baselines and contrasts ============================
# =====================================================================
baseline_map_full <- list(
  "vs Organoid (untreated)" = c(
    "T1_OrgBMDM1K_vs_Org",
    "T1_OrgBMDM5K_vs_Org",
    "T3_OrgTNFa_vs_Org",
    "T3_OrgTNFaBMDM1K_vs_Org",
    "T3_OrgTNFaBMDM5K_vs_Org"
  ),
  "vs Organoid + TNFa" = c(
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
  "vs Macrophage + TNFa" = c(
    "T6_OrgTNFaBMDM1K_vs_MacTNFa",
    "T6_OrgTNFaBMDM5K_vs_MacTNFa"
  )
)
all_contrasts_full <- unlist(baseline_map_full)
baseline_levels_full <- names(baseline_map_full)

# Layout A — basal, no TNFa: only 4 contrasts, 2 baselines
baseline_map_basal <- list(
  "vs Organoid (untreated)" = c(
    "T1_OrgBMDM1K_vs_Org",
    "T1_OrgBMDM5K_vs_Org"
  ),
  "vs Macrophage (untreated)" = c(
    "T2_OrgBMDM1K_vs_Mac",
    "T2_OrgBMDM5K_vs_Mac"
  )
)
all_contrasts_basal   <- unlist(baseline_map_basal)
baseline_levels_basal <- names(baseline_map_basal)

cat("Layout A (basal): ", length(all_contrasts_basal), " contrasts\n", sep = "")
cat("Layout B (full) : ", length(all_contrasts_full),  " contrasts\n", sep = "")

# Baseline color palette
baseline_palette_full <- c(
  "vs Organoid (untreated)"   = "#1F4E79",
  "vs Organoid + TNFa"        = "#6FA8DC",
  "vs Macrophage (untreated)" = "#A61C00",
  "vs Macrophage + TNFa"      = "#E06666"
)
baseline_palette_basal <- baseline_palette_full[c("vs Organoid (untreated)",
                                                  "vs Macrophage (untreated)")]

# =====================================================================
# ===== 4. Stem cell + lineage gene panel ============================
# =====================================================================
layered_panel <- list(
  ISC_markers = c(
    "Lgr5","Ascl2","Olfm4","Smoc2","Mex3a","Tnfrsf19",
    "Bmi1","Hopx","Tert","Lrig1","Msi1","Clu","Ly6a","Krt19",
    "Sox9","Prom1"
  ),
  Wnt_signaling = c(
    "Wnt3","Wnt3a","Wnt5a","Wnt2b",
    "Ctnnb1","Tcf7l2","Lef1","Apc","Gsk3b","Axin2",
    "Lrp5","Lrp6","Fzd5","Fzd7","Porcn",
    "Rnf43","Znrf3","Dkk1","Sfrp1","Notum","Nkd1",
    "Myc","Ccnd1","Cd44","Sox4","Sp5"
  ),
  Notch_signaling = c(
    "Notch1","Notch2","Notch3",
    "Dll1","Dll4","Jag1","Jag2",
    "Hes1","Hes5","Hey1","Heyl",
    "Rbpj","Maml1","Numb","Dlk1","Lfng","Mfng","Atoh1"
  ),
  Hippo_YAP = c(
    "Yap1","Wwtr1",
    "Tead1","Tead2","Tead3","Tead4",
    "Lats1","Lats2","Stk3","Stk4","Sav1","Mob1a","Mob1b","Nf2",
    "Amotl1","Amotl2","Ccn1","Ccn2",
    "Ctgf","Cyr61","Ankrd1","Birc5"
  ),
  Enterocyte = c(
    "Vil1","Alpi","Sis","Epcam",
    "Anpep","Mep1a","Lct","Maoa",
    "Apoa1","Apoa4","Apob","Fabp1","Fabp2",
    "Slc5a1","Slc15a1","Slc2a5","Slc7a8"
  ),
  Goblet = c(
    "Muc2","Muc3","Muc4","Muc13",
    "Tff3","Agr2","Fcgbp","Clca1","Zg16","Reg4",
    "Spdef","Klf4","Gfi1"
  ),
  Paneth = c(
    "Lyz1",
    "Defa1","Defa3","Defa4","Defa5","Defa17","Defa21","Defa22",
    "Mmp7","Pla2g2a","Ang4","Itln1",
    "Reg3g","Reg3b"
  ),
  Enteroendocrine = c(
    "Chga","Chgb",
    "Gcg","Pyy","Cck","Gip",
    "Sct","Sst","Nts","Ghrl",
    "Neurog3","Neurod1","Pax4","Pax6"
  )
)

panel_order <- names(layered_panel)
sc_all_genes <- unique(unlist(layered_panel))
cat("Stem cell panel: ", length(sc_all_genes), " unique genes across ",
    length(panel_order), " sub-panels\n", sep = "")

# Panel color palette
panel_palette <- c(
  ISC_markers     = "#2E8B57",
  Wnt_signaling   = "#4A90E2",
  Notch_signaling = "#8B4789",
  Hippo_YAP       = "#1F3A6E",
  Enterocyte      = "#8C510A",
  Goblet          = "#35978F",
  Paneth          = "#FDAE61",
  Enteroendocrine = "#B58900"
)

# Canonical first-panel assignment (for duplicates across panels)
gene_panel_unique <- do.call(rbind, lapply(panel_order, function(p) {
  data.frame(Panel = p, SYMBOL = layered_panel[[p]], stringsAsFactors = FALSE)
})) %>%
  dplyr::mutate(Panel = factor(Panel, levels = panel_order)) %>%
  dplyr::arrange(Panel) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)


# =====================================================================
# ===== 6. Run KEGG + GO ORA for all 14 contrasts ====================
# =====================================================================
cat("\n===== Running KEGG + GO ORA (BP, MF, CC) for all 14 contrasts =====\n")

# Blacklist (same as Step 7, IBD kept, irrelevant terms removed)
kegg_blacklist <- c(
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

filter_enrich <- function(res, blacklist) {
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

# Generic enrichment runners (KEGG / GO BP / GO MF / GO CC)
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
    res <- filter_enrich(res, kegg_blacklist)
  }
  res
}

run_go_ora <- function(deg_df, ont) {
  sig <- deg_df %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff,
                  abs(log2FoldChange) > lfc_cutoff) %>%
    dplyr::pull(ENTREZID) %>% na.omit() %>% unique() %>% as.character()
  if (length(sig) < 10) return(NULL)
  res <- tryCatch(
    enrichGO(gene = sig, OrgDb = org.Mm.eg.db,
             keyType = "ENTREZID", ont = ont,
             universe = universe_entrez,
             pvalueCutoff = 0.1, qvalueCutoff = 0.25,
             readable = TRUE),
    error = function(e) NULL
  )
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    # GO does not need disease blacklist; optional simplify
    res <- tryCatch(
      simplify(res, cutoff = 0.7, by = "p.adjust", select_fun = min),
      error = function(e) res
    )
  }
  res
}

# Run all enrichments
kegg_list  <- list()
go_bp_list <- list()
go_mf_list <- list()
go_cc_list <- list()

for (cn in all_contrasts_full) {
  cat("  ", cn, "\n", sep = "")
  kegg_list[[cn]]  <- run_kegg_ora(deg_gene[[cn]])
  go_bp_list[[cn]] <- run_go_ora(deg_gene[[cn]], ont = "BP")
  go_mf_list[[cn]] <- run_go_ora(deg_gene[[cn]], ont = "MF")
  go_cc_list[[cn]] <- run_go_ora(deg_gene[[cn]], ont = "CC")
}

# =====================================================================
# ===== 7. Contrast labels and strip labels (shared) =================
# =====================================================================
contrast_label_map <- c(
  T1_OrgBMDM1K_vs_Org         = "BMDM 1K",
  T1_OrgBMDM5K_vs_Org         = "BMDM 5K",
  T3_OrgTNFa_vs_Org           = "TNFa",
  T3_OrgTNFaBMDM1K_vs_Org     = "BMDM 1K+TNFa",
  T3_OrgTNFaBMDM5K_vs_Org     = "BMDM 5K+TNFa",
  T4_OrgTNFaBMDM1K_vs_OrgTNFa = "BMDM 1K+TNFa",
  T4_OrgTNFaBMDM5K_vs_OrgTNFa = "BMDM 5K+TNFa",
  T2_OrgBMDM1K_vs_Mac         = "Org+BMDM 1K",
  T2_OrgBMDM5K_vs_Mac         = "Org+BMDM 5K",
  T5_MacTNFa_vs_Mac           = "TNFa",
  T5_OrgTNFaBMDM1K_vs_Mac     = "Org+BMDM 1K+TNFa",
  T5_OrgTNFaBMDM5K_vs_Mac     = "Org+BMDM 5K+TNFa",
  T6_OrgTNFaBMDM1K_vs_MacTNFa = "Org+BMDM 1K+TNFa",
  T6_OrgTNFaBMDM5K_vs_MacTNFa = "Org+BMDM 5K+TNFa"
)

baseline_strip_labels <- c(
  "vs Organoid (untreated)"   = "vs Organoid (untreated)",
  "vs Organoid + TNFa"        = "vs Organoid + TNFa",
  "vs Macrophage (untreated)" = "vs Macrophage (untreated)",
  "vs Macrophage + TNFa"      = "vs Macrophage + TNFa"
)

# Check ggh4x for colored strips
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
if (!has_ggh4x) {
  install.packages("ggh4x")
  has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)
}

# =====================================================================
# ===== 8. Bubble plot builder (re-usable function) ==================
# =====================================================================
build_bubble_long <- function(enrich_list, baseline_map) {
  # Build contrast -> baseline map
  c2b <- list()
  for (bl in names(baseline_map)) {
    for (cn in baseline_map[[bl]]) c2b[[cn]] <- bl
  }
  rows <- list()
  for (cn in names(c2b)) {
    r <- enrich_list[[cn]]
    if (is.null(r) || nrow(as.data.frame(r)) == 0) next
    df <- as.data.frame(r)
    df$Contrast <- cn
    df$Baseline <- c2b[[cn]]
    rows[[cn]] <- df[, c("ID","Description","Count","pvalue","p.adjust",
                         "Contrast","Baseline")]
  }
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  out$logp <- -log10(out$pvalue)
  out
}

make_bubble_plot <- function(enrich_list,
                             baseline_map,
                             baseline_levels,
                             baseline_palette,
                             top_n,
                             title_main,
                             title_sub) {
  long <- build_bubble_long(enrich_list, baseline_map)
  if (is.null(long)) {
    message("No enrichment data; returning NULL")
    return(NULL)
  }
  long$Baseline <- factor(long$Baseline, levels = baseline_levels)
  
  # Top N per baseline
  selected <- long %>%
    dplyr::group_by(Baseline, Description) %>%
    dplyr::summarise(best_logp = max(logp), .groups = "drop") %>%
    dplyr::group_by(Baseline) %>%
    dplyr::arrange(dplyr::desc(best_logp), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  
  pathway_primary <- selected %>%
    dplyr::group_by(Description) %>%
    dplyr::slice_max(best_logp, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Baseline, dplyr::desc(best_logp))
  
  sub <- long %>%
    dplyr::filter(Description %in% pathway_primary$Description)
  sub$Description <- factor(sub$Description,
                            levels = rev(pathway_primary$Description))
  
  contrast_order <- unlist(baseline_map, use.names = FALSE)
  sub$Contrast <- factor(sub$Contrast, levels = contrast_order)
  
  # Build plot
  p <- ggplot(sub, aes(x = Contrast, y = Description)) +
    geom_point(aes(size = Count, fill = logp),
               shape = 21, color = "black", stroke = 0.4, alpha = 0.95)
  
  if (has_ggh4x) {
    strip_fills <- unname(baseline_palette[baseline_levels])
    p <- p + ggh4x::facet_grid2(
      . ~ Baseline,
      scales = "free_x", space = "free_x",
      labeller = labeller(Baseline = baseline_strip_labels),
      strip = ggh4x::strip_themed(
        background_x = lapply(strip_fills, function(col)
          element_rect(fill = col, color = "black", linewidth = 0.8)),
        text_x = lapply(strip_fills, function(col) {
          text_col <- ifelse(col %in% c("#1F4E79","#A61C00"), "white", "black")
          element_text(color = text_col, face = "bold", size = 11,
                       family = plot_font, lineheight = 0.95)
        })
      )
    )
  } else {
    p <- p + facet_grid(. ~ Baseline, scales = "free_x", space = "free_x",
                        labeller = labeller(Baseline = baseline_strip_labels))
  }
  
  p <- p +
    scale_fill_gradient2(
      low = "#5865f2", mid = "#faffff", high = "#ff3333",
      midpoint = 5,
      name = expression(-log[10]~italic(P)),
      limits = c(0, NA),
      guide = guide_colorbar(barwidth = 1, barheight = 8)
    ) +
    scale_size_continuous(range = c(2, 8), name = "Gene count",
                          breaks = c(30, 60, 90)) +
    scale_x_discrete(labels = contrast_label_map) +
    labs(x = NULL, y = NULL,
         title = title_main, subtitle = title_sub) +
    theme_bw(base_size = 11, base_family = plot_font) +
    theme(
      text                = element_text(family = plot_font),
      plot.title          = element_text(face = "bold", hjust = 0.5,
                                         size = 14, margin = margin(b = 4)),
      plot.subtitle       = element_text(hjust = 0.5, size = 10,
                                         color = "grey30",
                                         margin = margin(b = 12)),
      axis.text.x         = element_text(size = 8.5, face = "bold",
                                         color = "black",
                                         angle = 45, hjust = 1, vjust = 1),
      axis.text.y         = element_text(size = 9, color = "black"),
      axis.ticks          = element_line(color = "black", linewidth = 0.4),
      panel.grid.major    = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor    = element_blank(),
      panel.border        = element_rect(color = "black", fill = NA,
                                         linewidth = 0.6),
      panel.spacing.x     = unit(0.4, "lines"),
      strip.placement     = "outside",
      legend.position     = "right",
      legend.title        = element_text(face = "bold", size = 10),
      legend.text         = element_text(size = 9),
      legend.box          = "vertical",
      legend.spacing.y    = unit(0.3, "cm"),
      plot.margin         = margin(15, 15, 15, 15)
    )
  
  list(plot = p, data = sub, n_pathways = length(unique(sub$Description)))
}

# =====================================================================
# ===== 9. Generate bubble plots for both layouts =====================
# =====================================================================
generate_and_save_bubble <- function(enrich_list, layout_label,
                                     enrich_type, fig_id) {
  if (layout_label == "basal") {
    bmap <- baseline_map_basal
    bl   <- baseline_levels_basal
    pal  <- baseline_palette_basal
    top  <- top_n_basal
    sub_label <- "basal (no TNFa) - 4 contrasts"
    fig_prefix <- "Fig_5"
  } else {
    bmap <- baseline_map_full
    bl   <- baseline_levels_full
    pal  <- baseline_palette_full
    top  <- top_n_per_baseline
    sub_label <- "full 14 contrasts"
    fig_prefix <- "Fig_S5"
  }
  
  title_main <- paste0(enrich_type, " enrichment - ", sub_label)
  title_sub  <- paste0("Top ", top,
                       " per baseline (union); bubble size = gene count, color = -log10(p)")
  
  result <- make_bubble_plot(
    enrich_list, bmap, bl, pal, top, title_main, title_sub
  )
  
  if (is.null(result)) {
    cat("  [", layout_label, " / ", enrich_type, "] no data\n", sep = "")
    return(invisible(NULL))
  }
  
  fname <- paste0(fig_prefix, fig_id, "_",
                  gsub(" ", "_", enrich_type),
                  "_bubble_", layout_label)
  fig_w <- if (layout_label == "basal") 9  else 13
  fig_h <- max(8, result$n_pathways * 0.32 + 4)
  
  ggsave(file.path(out_root, paste0(fname, ".png")),
         result$plot, width = fig_w, height = fig_h,
         dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(out_root, paste0(fname, ".pdf")),
         result$plot, width = fig_w, height = fig_h,
         device = cairo_pdf, limitsize = FALSE)
  
  # Save table
  write.csv(result$data,
            file.path(out_root, "tables",
                      paste0(gsub(" ", "_", enrich_type), "_bubble_",
                             layout_label, ".csv")),
            row.names = FALSE)
  
  cat("  [", layout_label, " / ", enrich_type, "] ",
      result$n_pathways, " pathways -> ", fname, ".png\n", sep = "")
}

cat("\n===== Layout A (basal, no TNFa) =====\n")
generate_and_save_bubble(kegg_list,  "basal", "KEGG",  "f")
generate_and_save_bubble(go_bp_list, "basal", "GO BP", "g")
generate_and_save_bubble(go_mf_list, "basal", "GO MF", "h")
generate_and_save_bubble(go_cc_list, "basal", "GO CC", "i")

cat("\n===== Layout B (full 14 contrasts) =====\n")
generate_and_save_bubble(kegg_list,  "full", "KEGG",  "f")
generate_and_save_bubble(go_bp_list, "full", "GO BP", "g")
generate_and_save_bubble(go_mf_list, "full", "GO MF", "h")
generate_and_save_bubble(go_cc_list, "full", "GO CC", "i")

# =====================================================================
# ===== 10. Stem cell LFC heatmap (panel-balanced top 5 x 8 panels) ==
# =====================================================================
cat("\n===== Generating Fig 5j/S5j (Stem cell LFC heatmaps) =====\n")

# Helper: collect LFC for all panel genes across given contrasts
collect_lfc_long <- function(contrasts_set, gene_set) {
  rows <- list()
  for (cn in contrasts_set) {
    d <- deg_results[[cn]] %>%
      dplyr::filter(SYMBOL %in% gene_set, !is.na(padj))
    rows[[cn]] <- data.frame(
      SYMBOL    = d$SYMBOL,
      LFC       = d$log2FoldChange,
      padj      = d$padj,
      Contrast  = cn,
      stringsAsFactors = FALSE
    )
  }
  long <- do.call(rbind, rows) %>%
    dplyr::distinct(SYMBOL, Contrast, .keep_all = TRUE) %>%
    dplyr::left_join(gene_panel_unique, by = "SYMBOL") %>%
    dplyr::filter(!is.na(Panel))
  long
}

# Panel-balanced top-N selection: top N per panel by combined rank score
select_panel_balanced_top <- function(long_df, n_per_panel) {
  ranked <- long_df %>%
    dplyr::group_by(SYMBOL, Panel) %>%
    dplyr::summarise(
      max_abs_lfc = max(abs(LFC), na.rm = TRUE),
      min_padj    = min(padj, na.rm = TRUE),
      n_sig       = sum(padj < padj_cutoff & abs(LFC) > lfc_cutoff,
                        na.rm = TRUE),
      rank_score  = max_abs_lfc * (-log10(pmax(min_padj, 1e-300))),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(rank_score)) %>%
    dplyr::arrange(Panel, dplyr::desc(rank_score))
  
  top <- ranked %>%
    dplyr::group_by(Panel) %>%
    dplyr::slice_head(n = n_per_panel) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Panel, dplyr::desc(rank_score))
  
  list(top = top, ranked = ranked)
}

# Build LFC matrix from long
build_lfc_matrix <- function(top_genes, contrasts_set) {
  lfc_mat  <- matrix(NA, nrow = length(top_genes), ncol = length(contrasts_set),
                     dimnames = list(top_genes, contrasts_set))
  padj_mat <- lfc_mat
  for (i in seq_along(top_genes)) {
    for (j in seq_along(contrasts_set)) {
      d <- deg_results[[contrasts_set[j]]] %>%
        dplyr::filter(SYMBOL == top_genes[i])
      if (nrow(d) > 0) {
        lfc_mat[i, j]  <- d$log2FoldChange[1]
        padj_mat[i, j] <- d$padj[1]
      }
    }
  }
  list(lfc = lfc_mat, padj = padj_mat)
}

# Heatmap column labels
col_labels_full <- c(
  T1_OrgBMDM1K_vs_Org         = "BMDM 1K\nvs Org",
  T1_OrgBMDM5K_vs_Org         = "BMDM 5K\nvs Org",
  T3_OrgTNFa_vs_Org           = "TNFa\nvs Org",
  T3_OrgTNFaBMDM1K_vs_Org     = "TNFa+BMDM 1K\nvs Org",
  T3_OrgTNFaBMDM5K_vs_Org     = "TNFa+BMDM 5K\nvs Org",
  T4_OrgTNFaBMDM1K_vs_OrgTNFa = "BMDM 1K+TNFa\nvs Org+TNFa",
  T4_OrgTNFaBMDM5K_vs_OrgTNFa = "BMDM 5K+TNFa\nvs Org+TNFa",
  T2_OrgBMDM1K_vs_Mac         = "Org+BMDM 1K\nvs Mac",
  T2_OrgBMDM5K_vs_Mac         = "Org+BMDM 5K\nvs Mac",
  T5_MacTNFa_vs_Mac           = "Mac+TNFa\nvs Mac",
  T5_OrgTNFaBMDM1K_vs_Mac     = "Org+TNFa+BMDM 1K\nvs Mac",
  T5_OrgTNFaBMDM5K_vs_Mac     = "Org+TNFa+BMDM 5K\nvs Mac",
  T6_OrgTNFaBMDM1K_vs_MacTNFa = "Org+BMDM 1K+TNFa\nvs Mac+TNFa",
  T6_OrgTNFaBMDM5K_vs_MacTNFa = "Org+BMDM 5K+TNFa\nvs Mac+TNFa"
)

# =====================================================================
# ===== 10. Stem cell LFC heatmap — 2분할 × 2 layout = 4 heatmap ======
# =====================================================================
cat("\n===== Generating split heatmaps (4 total, print 가능) =====\n")

# ----- Heatmap 분할 정의 ---------------------------------------------
heatmap_split <- list(
  Heatmap1_StemSignaling  = c("ISC_markers","Wnt_signaling","Notch_signaling"),
  Heatmap2_LineageOutcome = c("Hippo_YAP","Enterocyte","Goblet","Paneth","Enteroendocrine")
)

# ----- 출력 폴더 -----------------------------------------------------
supply_root <- "supply_stemness"
dir.create(supply_root, showWarnings = FALSE)
dir.create(file.path(supply_root, "tables"), showWarnings = FALSE)

# ----- 폰트 / 셀 크기 — 튜닝 포인트 ----------------------------------
hm_params <- list(
  fontsize_row = 11,    # ★ gene symbol
  fontsize_col = 10,    # ★ column label
  cellwidth    = 30,    # ★ 셀 너비
  cellheight   = 24,    # ★ 셀 높이
  fontsize     = 12     # ★ title / legend
)

# ----- Column label sets ---------------------------------------------
col_labels_basal <- c(
  T1_OrgBMDM1K_vs_Org = "BMDM 1K\nvs Org",
  T1_OrgBMDM5K_vs_Org = "BMDM 5K\nvs Org",
  T2_OrgBMDM1K_vs_Mac = "Org+BMDM 1K\nvs Mac",
  T2_OrgBMDM5K_vs_Mac = "Org+BMDM 5K\nvs Mac"
)
# col_labels_full 은 기존 코드에 이미 정의되어 있음

# =====================================================================
# ----- Heatmap 빌드 함수 (객체 반환, 자동 화면 출력/저장 X) ----------
# =====================================================================
build_split_heatmap <- function(panels_subset,
                                contrasts_set,
                                baseline_map,
                                baseline_levels,
                                baseline_palette,
                                layout_label,
                                split_name,
                                col_label_set,
                                gaps_col_pattern = NULL,
                                fs = hm_params) {
  
  # 1) 해당 split의 gene set
  panel_genes_subset <- unique(unlist(layered_panel[panels_subset]))
  
  # 2) LFC 수집
  long <- collect_lfc_long(contrasts_set, panel_genes_subset)
  if (nrow(long) == 0) {
    cat("  [", split_name, " / ", layout_label, "] no data\n", sep = "")
    return(NULL)
  }
  
  # 3) Panel-balanced top N
  sel <- select_panel_balanced_top(long, top_n_per_panel)
  top_genes_ordered <- sel$top$SYMBOL
  
  # 4) LFC matrix
  mats <- build_lfc_matrix(top_genes_ordered, contrasts_set)
  lfc_mat  <- mats$lfc
  padj_mat <- mats$padj
  
  keep_rows <- rowSums(!is.na(lfc_mat)) > 0
  lfc_mat  <- lfc_mat[keep_rows, , drop = FALSE]
  padj_mat <- padj_mat[keep_rows, , drop = FALSE]
  top_keep <- sel$top %>% dplyr::filter(SYMBOL %in% rownames(lfc_mat))
  
  # 5) Column labels
  colnames(lfc_mat) <- col_label_set[colnames(lfc_mat)]
  
  # 6) Column annotation (명시적 순서 매칭)
  c2b <- list()
  for (bl in names(baseline_map)) {
    for (cn in baseline_map[[bl]]) c2b[[cn]] <- bl
  }
  col_baselines <- unlist(c2b[contrasts_set], use.names = FALSE)
  col_group <- data.frame(
    Baseline = factor(col_baselines, levels = baseline_levels)
  )
  rownames(col_group) <- colnames(lfc_mat)
  
  # 7) Row annotation
  row_panels <- top_keep$Panel
  row_group <- data.frame(
    Panel = factor(row_panels, levels = panels_subset)
  )
  rownames(row_group) <- rownames(lfc_mat)
  
  ann_color <- list(
    Baseline = baseline_palette[baseline_levels],
    Panel    = panel_palette[panels_subset]
  )
  
  # 8) Color scale
  lfc_max <- max(abs(lfc_mat), na.rm = TRUE)
  lfc_max <- min(lfc_max, 12)
  breaks  <- seq(-lfc_max, lfc_max, length.out = 101)
  
  # 9) Row gaps (panel 경계)
  panel_counts <- table(factor(row_panels, levels = panels_subset))
  panel_cumsum <- cumsum(as.integer(panel_counts))
  gaps_row <- as.integer(unlist(lapply(head(panel_cumsum, -1),
                                       function(x) rep(x, 2))))
  
  gaps_col <- if (is.null(gaps_col_pattern)) NULL else gaps_col_pattern
  
  # 10) pheatmap 빌드 — ComplexHeatmap 마스킹 방지를 위해 namespace 명시
  p <- pheatmap::pheatmap(
    lfc_mat,
    color = colorRampPalette(c("#5865f2","#e8e8e8","#ff3232"))(100),
    breaks = breaks,
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_col = col_group,
    annotation_row = row_group,
    annotation_colors = ann_color,
    annotation_names_col = FALSE,
    annotation_names_row = FALSE,
    gaps_col = gaps_col,
    gaps_row = gaps_row,
    display_numbers = FALSE,
    fontsize        = fs$fontsize,
    fontsize_row    = fs$fontsize_row,
    fontsize_col    = fs$fontsize_col,
    cellwidth       = fs$cellwidth,
    cellheight      = fs$cellheight,
    angle_col       = "45",
    border_color    = "black",
    lwd             = 1,
    main            = paste0(split_name, " — LFC (", layout_label, ")"),
    silent          = TRUE       # ★ 자동 화면 출력 OFF
  )
  
  list(plot = p, lfc = lfc_mat, padj = padj_mat,
       row_panels = row_panels, ranked = sel$ranked,
       fig_name = paste0(split_name, "_", layout_label))
}

# =====================================================================
# ----- 저장 함수 -----------------------------------------------------
# =====================================================================
save_split_heatmap <- function(res, fig_w = NULL, fig_h = NULL) {
  if (is.null(res)) return(invisible(NULL))
  
  if (is.null(fig_h)) fig_h <- max(8, nrow(res$lfc) * 0.35 + 4)
  if (is.null(fig_w)) fig_w <- max(8, ncol(res$lfc) * 0.9 + 4)
  
  fn <- res$fig_name
  
  # PNG
  ggsave(file.path(supply_root, paste0(fn, ".png")),
         res$plot$gtable, width = fig_w, height = fig_h,
         dpi = 300, bg = "white", limitsize = FALSE)
  # PDF
  ggsave(file.path(supply_root, paste0(fn, ".pdf")),
         res$plot$gtable, width = fig_w, height = fig_h,
         device = cairo_pdf, limitsize = FALSE)
  # TIFF
  ggsave(file.path(supply_root, paste0(fn, ".tiff")),
         res$plot$gtable, width = fig_w, height = fig_h,
         device = "tiff", dpi = 300, bg = "white",
         compression = "lzw", limitsize = FALSE)
  
  # Tables
  lfc_df_out  <- data.frame(SYMBOL = rownames(res$lfc),
                            Panel = res$row_panels,
                            res$lfc, check.names = FALSE)
  padj_df_out <- data.frame(SYMBOL = rownames(res$padj),
                            Panel = res$row_panels,
                            res$padj, check.names = FALSE)
  write.csv(lfc_df_out,
            file.path(supply_root, "tables", paste0(fn, "_LFC.csv")),
            row.names = FALSE)
  write.csv(padj_df_out,
            file.path(supply_root, "tables", paste0(fn, "_padj.csv")),
            row.names = FALSE)
  write.csv(res$ranked,
            file.path(supply_root, "tables", paste0(fn, "_ranking.csv")),
            row.names = FALSE)
  
  cat("  Saved: ", fn, ".{png,pdf,tiff} (",
      fig_w, " × ", fig_h, " in)\n", sep = "")
}

# =====================================================================
# ===== Heatmap 1 (basal) — ISC + Wnt + Notch, no TNFα ===============
# =====================================================================
cat("\n--- Heatmap 1 (basal) ---\n")
res_hm1_basal <- build_split_heatmap(
  panels_subset    = heatmap_split$Heatmap1_StemSignaling,
  contrasts_set    = all_contrasts_basal,
  baseline_map     = baseline_map_basal,
  baseline_levels  = baseline_levels_basal,
  baseline_palette = baseline_palette_basal,
  layout_label     = "basal",
  split_name       = "Heatmap1_StemSignaling",
  col_label_set    = col_labels_basal,
  gaps_col_pattern = c(2, 2, 2)
)

while (!is.null(dev.list())) dev.off()
grid::grid.newpage(); grid::grid.draw(res_hm1_basal$plot$gtable)

# ▼ 마음에 들면 저장 ▼
# save_split_heatmap(res_hm1_basal, fig_w = 10, fig_h = 14)


# =====================================================================
# ===== Heatmap 1 (full) — ISC + Wnt + Notch, 14 contrasts ===========
# =====================================================================
cat("\n--- Heatmap 1 (full) ---\n")
res_hm1_full <- build_split_heatmap(
  panels_subset    = heatmap_split$Heatmap1_StemSignaling,
  contrasts_set    = all_contrasts_full,
  baseline_map     = baseline_map_full,
  baseline_levels  = baseline_levels_full,
  baseline_palette = baseline_palette_full,
  layout_label     = "full",
  split_name       = "Heatmap1_StemSignaling",
  col_label_set    = col_labels_full,
  gaps_col_pattern = c(5, 5, 5, 7, 7, 7, 12, 12, 12)
)

while (!is.null(dev.list())) dev.off()
grid::grid.newpage(); grid::grid.draw(res_hm1_full$plot$gtable)

# ▼ 마음에 들면 저장 ▼
# save_split_heatmap(res_hm1_full, fig_w = 18, fig_h = 14)


# =====================================================================
# ===== Heatmap 2 (basal) — Hippo + Ent + Gob + Pan + EE =============
# =====================================================================
cat("\n--- Heatmap 2 (basal) ---\n")
res_hm2_basal <- build_split_heatmap(
  panels_subset    = heatmap_split$Heatmap2_LineageOutcome,
  contrasts_set    = all_contrasts_basal,
  baseline_map     = baseline_map_basal,
  baseline_levels  = baseline_levels_basal,
  baseline_palette = baseline_palette_basal,
  layout_label     = "basal",
  split_name       = "Heatmap2_LineageOutcome",
  col_label_set    = col_labels_basal,
  gaps_col_pattern = c(2, 2, 2)
)

while (!is.null(dev.list())) dev.off()
grid::grid.newpage(); grid::grid.draw(res_hm2_basal$plot$gtable)

# ▼ 마음에 들면 저장 ▼
# save_split_heatmap(res_hm2_basal, fig_w = 10, fig_h = 16)


# =====================================================================
# ===== Heatmap 2 (full) — Lineage panels, 14 contrasts ==============
# =====================================================================
cat("\n--- Heatmap 2 (full) ---\n")
res_hm2_full <- build_split_heatmap(
  panels_subset    = heatmap_split$Heatmap2_LineageOutcome,
  contrasts_set    = all_contrasts_full,
  baseline_map     = baseline_map_full,
  baseline_levels  = baseline_levels_full,
  baseline_palette = baseline_palette_full,
  layout_label     = "full",
  split_name       = "Heatmap2_LineageOutcome",
  col_label_set    = col_labels_full,
  gaps_col_pattern = c(5, 5, 5, 7, 7, 7, 12, 12, 12)
)

while (!is.null(dev.list())) dev.off()
grid::grid.newpage(); grid::grid.draw(res_hm2_full$plot$gtable)

# ▼ 마음에 들면 저장 ▼
# save_split_heatmap(res_hm2_full, fig_w = 18, fig_h = 16)


# =====================================================================
# ===== 11. Summary ===================================================
# =====================================================================
cat("\n===== Step 8 complete (stem cell + lineage analysis) =====\n")
cat("Output roots:\n")
cat("  ", out_root,    "/    — bubble plots + Venn\n", sep = "")
cat("  ", supply_root, "/    — split heatmaps (4 total)\n\n", sep = "")
cat("Layout A (basal, no TNFα) - main figure:\n")
cat("  Fig_5f_KEGG_bubble_basal.{png,pdf}\n")
cat("  Fig_5g_GO_BP_bubble_basal.{png,pdf}\n")
cat("  Fig_5h_GO_MF_bubble_basal.{png,pdf}\n")
cat("  Fig_5i_GO_CC_bubble_basal.{png,pdf}\n\n")
cat("Layout B (full 14 contrasts) - supplementary:\n")
cat("  Fig_S5f_KEGG_bubble_full.{png,pdf}\n")
cat("  Fig_S5g_GO_BP_bubble_full.{png,pdf}\n")
cat("  Fig_S5h_GO_MF_bubble_full.{png,pdf}\n")
cat("  Fig_S5i_GO_CC_bubble_full.{png,pdf}\n\n")
cat("Stem cell heatmaps (split 2 × layout 2 = 4 total):\n")
cat("  Heatmap1_StemSignaling_basal.{png,pdf,tiff}   — ISC + Wnt + Notch (no TNFα)\n")
cat("  Heatmap1_StemSignaling_full.{png,pdf,tiff}    — ISC + Wnt + Notch (14 contrasts)\n")
cat("  Heatmap2_LineageOutcome_basal.{png,pdf,tiff}  — Hippo + Ent + Gob + Pan + EE (no TNFα)\n")
cat("  Heatmap2_LineageOutcome_full.{png,pdf,tiff}   — Hippo + Ent + Gob + Pan + EE (14 contrasts)\n\n")
cat("Tables in supply_stemness/tables/:\n")
cat("  Heatmap{1,2}_*_{basal,full}_{LFC,padj,ranking}.csv\n\n")
cat("Gene panel split:\n")
cat("  Heatmap 1 (Stem + Signaling):\n")
for (p in heatmap_split$Heatmap1_StemSignaling) {
  cat("    ", p, " (", length(layered_panel[[p]]), " genes)\n", sep = "")
}
cat("  Heatmap 2 (Lineage Outcome):\n")
for (p in heatmap_split$Heatmap2_LineageOutcome) {
  cat("    ", p, " (", length(layered_panel[[p]]), " genes)\n", sep = "")
}
cat("\nMechanistic narrative:\n")
cat("  Heatmap 1 - upstream stem cell signaling input\n")
cat("    → Does BMDM co-culture modulate ISC pool & Wnt/Notch fate-decision input?\n")
cat("  Heatmap 2 - downstream lineage commitment outcome\n")
cat("    → Does the upstream change translate into altered differentiation outcomes?\n")
cat("  Compare basal vs full layouts:\n")
cat("    - basal layout: BMDM effect alone\n")
cat("    - full layout: BMDM effect under TNFα-conditioned environment\n")