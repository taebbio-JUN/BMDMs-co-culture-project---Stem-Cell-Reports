# =====================================================================
# Step 4: Recreate Nature 2025 Figure 1e–g style figures
# Project: HN00273522 — BMDM × Macrophage × TNFα × Organoid co-culture
# Reference: doi.org/10.1038/s41467-025-59639-9 Figure 1
#
# Outputs:
#   Fig_1e_4way_Venn.png          — 4-way Venn of expressed genes
#   Fig_1f_KEGG_paired_bars.png   — Tier 1 + Tier 2 KEGG paired bars
#   Fig_1g_TopInflam_LFC_heatmap.png — Top 10 inflam gene LFC heatmap
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
  library(clusterProfiler); library(enrichplot); library(DOSE)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

# ===== 1. Configuration ==============================================
out_root <- "Step2_Fig1_recreate"
plot_font <- "Helvetica"

dir.create(out_root, showWarnings = FALSE)
dir.create(file.path(out_root, "tables"), showWarnings = FALSE)

padj_cutoff <- 0.05
lfc_cutoff  <- 1.0

# 8-group color scheme
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

# ===== 2. Load Step 1 outputs ========================================
dds            <- readRDS("dds_full_8groups.rds")
vsd            <- readRDS("vsd_full_8groups.rds")
deg_results    <- readRDS("deg_results_full_8groups.rds")
anno           <- readRDS("gene_annotation.rds")
contrast_tiers <- readRDS("contrast_tiers.rds")

cat("Loaded: ", nrow(dds), " genes × ", ncol(dds), " samples\n", sep = "")

# =====================================================================
# ===== 3. Define "expressed" gene set per group (Fig 1e) ============
# =====================================================================
# Criteria: normalized count ≥ 10 in ≥ 2 out of 3 replicates
# (Applied AFTER DESeq2 pre-filter, so noise is minimal)

cat("\n===== Defining expressed gene sets per group =====\n")

# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE)
cat("Normalized count matrix: ", nrow(norm_counts), " genes × ",
    ncol(norm_counts), " samples\n", sep = "")

# Target groups for Venn
venn_groups <- c("Organoid", "Macrophage",
                 "Organoid_BMDM_1K", "Organoid_BMDM_5K")
venn_labels <- c("Intestinal Organoid",
                 "BMDMs",
                 "Intestinal Organoid +\nBMDMs (1K)",
                 "Intestinal Organoid +\nBMDMs (5K)")

# Build per-group expressed gene set
group_vec <- as.character(colData(dds)$group)

build_expressed_set <- function(group_name, threshold = 10, min_n = 2) {
  cols <- which(group_vec == group_name)
  if (length(cols) < 1) {
    stop("No samples found for group: ", group_name)
  }
  # ≥ threshold count in ≥ min_n samples
  expressed <- rowSums(norm_counts[, cols, drop = FALSE] >= threshold) >= min_n
  rownames(norm_counts)[expressed]
}

venn_sets <- list()
for (i in seq_along(venn_groups)) {
  venn_sets[[venn_labels[i]]] <- build_expressed_set(venn_groups[i])
  cat("  ", venn_groups[i], ": ", length(venn_sets[[venn_labels[i]]]),
      " expressed genes\n", sep = "")
}

# Save expression status per gene (long format for reference)
all_expressed <- unique(unlist(venn_sets))
venn_membership <- data.frame(
  gene_id = all_expressed,
  Organoid                            = all_expressed %in% venn_sets[[venn_labels[1]]],
  BMDMs                               = all_expressed %in% venn_sets[[venn_labels[2]]],
  Organoid_BMDM_1K                    = all_expressed %in% venn_sets[[venn_labels[3]]],
  Organoid_BMDM_5K                    = all_expressed %in% venn_sets[[venn_labels[4]]]
)
# Annotate with SYMBOL
venn_membership <- merge(venn_membership, anno[, c("gene_id","SYMBOL","ENTREZID")],
                         by = "gene_id", all.x = TRUE)
write.csv(venn_membership,
          file.path(out_root, "tables", "Venn_membership_per_gene.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 4. Figure 1e — 4-way Venn diagram =============================
# =====================================================================
cat("\n===== Generating Figure 1e (4-way Venn) =====\n")

# Try ggvenn first; fallback to VennDiagram if needed
has_ggvenn <- requireNamespace("ggvenn", quietly = TRUE)
has_VD     <- requireNamespace("VennDiagram", quietly = TRUE)

if (!has_ggvenn && !has_VD) {
  cat("[NOTE] Installing ggvenn...\n")
  install.packages("ggvenn")
  has_ggvenn <- requireNamespace("ggvenn", quietly = TRUE)
}

venn_palette <- c("#A6CEE3","#FB9A99","#B2DF8A","#FDBF6F")  # 4 groups

if (has_ggvenn) {
  library(ggvenn)
  
  p_venn <- ggvenn(
    venn_sets,
    fill_color  = venn_palette,
    fill_alpha  = 0.55,
    stroke_color = "black",
    stroke_size = 0.5,
    set_name_size = 5,
    text_size   = 4.5,
    show_percentage = FALSE
  ) +
    ggtitle("Expressed genes — 4-way Venn") +
    theme(plot.title = element_text(family = plot_font, face = "bold",
                                    hjust = 0.5, size = 16))
  
  ggsave(file.path(out_root, "Fig_1e_4way_Venn.png"), p_venn,
         width = 8.5, height = 8.5, dpi = 300, bg = "white")
  ggsave(file.path(out_root, "Fig_1e_4way_Venn.pdf"), p_venn,
         width = 8.5, height = 8.5, device = cairo_pdf)
  
} else if (has_VD) {
  # Fallback: VennDiagram package
  library(VennDiagram)
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  png(file.path(out_root, "Fig_1e_4way_Venn.png"),
      width = 2400, height = 2400, res = 300, bg = "white")
  grid.draw(venn.diagram(
    x = venn_sets,
    category.names = venn_labels,
    filename = NULL,
    fill = venn_palette,
    alpha = 0.55,
    cex = 1.3,
    cat.cex = 1.1,
    cat.fontface = "bold",
    cat.col = "black",
    margin = 0.1,
    main = "Expressed genes — 4-way Venn",
    main.cex = 1.5,
    main.fontface = "bold"
  ))
  dev.off()
}

# Print intersection sizes for narrative
cat("\n===== Venn intersection summary =====\n")
core_all4 <- Reduce(intersect, venn_sets)
cat("Core (all 4 groups)              : ", length(core_all4), "\n", sep = "")
cat("Organoid-unique                  : ",
    length(setdiff(venn_sets[[venn_labels[1]]],
                   unique(unlist(venn_sets[-1])))), "\n", sep = "")
cat("BMDMs-unique                     : ",
    length(setdiff(venn_sets[[venn_labels[2]]],
                   unique(unlist(venn_sets[-2])))), "\n", sep = "")
cat("Org+BMDM 1K unique               : ",
    length(setdiff(venn_sets[[venn_labels[3]]],
                   unique(unlist(venn_sets[-3])))), "\n", sep = "")
cat("Org+BMDM 5K unique               : ",
    length(setdiff(venn_sets[[venn_labels[4]]],
                   unique(unlist(venn_sets[-4])))), "\n", sep = "")
cat("Shared with Macrophage           : ",
    length(intersect(venn_sets[[venn_labels[2]]],
                     intersect(venn_sets[[venn_labels[3]]],
                               venn_sets[[venn_labels[4]]]))),
    "\n", sep = "")

# =====================================================================
# ===== 5. KEGG ORA for Tier 1 + Tier 2 (Fig 1f preparation) =========
# =====================================================================
cat("\n===== Running KEGG ORA for Tier 1 + Tier 2 =====\n")

# Comprehensive disease pathway BLACKLIST (virus, bacteria, parasitic, cancer,
# cardio, neurodegenerative — all non-physiological disease contexts)
kegg_blacklist <- c(
  # ---- Virus ----
  "Herpes simplex virus 1 infection","Epstein-Barr virus infection",
  "Kaposi sarcoma-associated herpesvirus infection",
  "Human papillomavirus infection","Human cytomegalovirus infection",
  "Human immunodeficiency virus 1 infection",
  "Human T-cell leukemia virus 1 infection",
  "Hepatitis B","Hepatitis C","Measles","Influenza A",
  "Coronavirus disease - COVID-19","Viral myocarditis",
  "Viral carcinogenesis","Viral life cycle - HIV-1",
  "Viral protein interaction with cytokine and cytokine receptor",
  # ---- Bacterial ----
  "Staphylococcus aureus infection","Tuberculosis","Pertussis",
  "Legionellosis","Salmonella infection","Yersinia infection",
  "Vibrio cholerae infection","Pathogenic Escherichia coli infection",
  "Shigellosis","Bacterial invasion of epithelial cells",
  "Epithelial cell signaling in Helicobacter pylori infection",
  # ---- Parasitic ----
  "Leishmaniasis","Malaria","Chagas disease",
  "African trypanosomiasis","Amoebiasis","Toxoplasmosis",
  # ---- Cancer (mostly not relevant to organoid-BMDM co-culture) ----
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
  # ---- Cardiovascular ----
  "Dilated cardiomyopathy","Hypertrophic cardiomyopathy",
  "Arrhythmogenic right ventricular cardiomyopathy",
  "Fluid shear stress and atherosclerosis","Lipid and atherosclerosis",
  "AGE-RAGE signaling pathway in diabetic complications",
  # ---- Neurodegenerative ----
  "Alzheimer disease","Parkinson disease","Huntington disease",
  "Amyotrophic lateral sclerosis","Prion disease",
  "Pathways of neurodegeneration - multiple diseases",
  # ---- Autoimmune (not relevant) ----
  "Rheumatoid arthritis","Systemic lupus erythematosus",
  "Type I diabetes mellitus","Autoimmune thyroid disease",
  "Graft-versus-host disease","Allograft rejection","Asthma",
  "Primary immunodeficiency","Inflammatory bowel disease"
  # ↑ IBD를 뺄지 말지는 manuscript narrative에 따라 결정. 
  # organoid+BMDM co-culture는 IBD model이므로 보통은 *keep*하는 게 자연스러우나,
  # 사용자가 "disease pathway 제거"를 강조했으므로 일단 제외.
  # IBD를 다시 살리고 싶으면 위 줄 끝에 # 주석 처리하시면 됩니다.
)

# Apply blacklist
filter_kegg <- function(res, blacklist) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(res)
  res_new <- res
  res_new@result <- res@result[!res@result$Description %in% blacklist, ]
  res_new
}

# Collapse DEG to gene level (ENTREZID-based)
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

# Tier 1 + Tier 2 contrasts only
fig1f_contrasts <- c("T1_OrgBMDM1K_vs_Org", "T1_OrgBMDM5K_vs_Org",
                     "T2_OrgBMDM1K_vs_Mac", "T2_OrgBMDM5K_vs_Mac")

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
for (cn in fig1f_contrasts) {
  cat("  ", cn, "\n", sep = "")
  kegg_list[[cn]] <- run_kegg_ora(deg_gene[[cn]])
}

# =====================================================================
# ===== 6. Figure 1f — Stacked 2-panel KEGG bar plot ================
# =====================================================================
cat("\n===== Generating Figure 1f (stacked 2-panel KEGG bars) =====\n")

# Build long-format data
build_kegg_long <- function(kegg_list, contrasts_t1, contrasts_t2) {
  rows <- list()
  for (cn in c(contrasts_t1, contrasts_t2)) {
    r <- kegg_list[[cn]]
    if (is.null(r) || nrow(as.data.frame(r)) == 0) next
    df <- as.data.frame(r)
    df$Contrast <- cn
    df$Comparison <- ifelse(cn %in% contrasts_t1, "vs Organoid", "vs Macrophage")
    df$BMDM_dose  <- ifelse(grepl("BMDM1K", cn), "1K", "5K")
    rows[[cn]] <- df[, c("ID","Description","Count","pvalue","p.adjust",
                         "Contrast","Comparison","BMDM_dose")]
  }
  do.call(rbind, rows)
}

kegg_long <- build_kegg_long(kegg_list,
                             c("T1_OrgBMDM1K_vs_Org","T1_OrgBMDM5K_vs_Org"),
                             c("T2_OrgBMDM1K_vs_Mac","T2_OrgBMDM5K_vs_Mac"))

# Aggregate within Comparison (combine 1K and 5K → take more significant)
kegg_panel <- kegg_long %>%
  dplyr::group_by(ID, Description, Comparison) %>%
  dplyr::summarise(
    pvalue   = min(pvalue, na.rm = TRUE),
    p.adjust = min(p.adjust, na.rm = TRUE),
    Count    = max(Count, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(logp = -log10(pvalue))

# Order: within each Comparison panel, by -log10(p) descending
# Build per-panel data frames separately to allow independent ordering
panel_top <- kegg_panel %>%
  dplyr::filter(Comparison == "vs Organoid") %>%
  dplyr::arrange(dplyr::desc(logp))
panel_bot <- kegg_panel %>%
  dplyr::filter(Comparison == "vs Macrophage") %>%
  dplyr::arrange(dplyr::desc(logp))

# Optional: limit to top N per panel for clarity (default = all)
top_n_per_panel <- 20    # ← 조정 가능
panel_top <- head(panel_top, top_n_per_panel)
panel_bot <- head(panel_bot, top_n_per_panel)

# Combine, with Comparison as factor for facet order
kegg_combined <- dplyr::bind_rows(panel_top, panel_bot)
kegg_combined$Comparison <- factor(kegg_combined$Comparison,
                                   levels = c("vs Organoid","vs Macrophage"))

# Within-panel ordering: each panel orders independently
# Trick: create a unique panel-specific factor
kegg_combined <- kegg_combined %>%
  dplyr::group_by(Comparison) %>%
  dplyr::arrange(Comparison, logp) %>%   # ascending so highest is at top
  dplyr::mutate(row_id = paste(Comparison, Description, sep = "__")) %>%
  dplyr::ungroup()

kegg_combined$row_id <- factor(kegg_combined$row_id,
                               levels = kegg_combined$row_id)

# Color palette: pink (vs Organoid), teal (vs Macrophage)
comparison_palette <- c(
  "vs Organoid"   = "#F4A6A6",
  "vs Macrophage" = "#A8D5C2"
)

# Build plot
p_1f <- ggplot(kegg_combined,
               aes(x = logp, y = row_id, fill = Comparison)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_point(aes(size = Count), shape = 21,
             color = "black", fill = "white", stroke = 0.5) +
  # Use Description as y-axis labels (strip the row_id prefix)
  scale_y_discrete(labels = function(x) sub("^[^_]+__", "", x)) +
  facet_grid(Comparison ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_fill_manual(values = comparison_palette,
                    name = "MaugOs vs") +
  scale_size_continuous(range = c(2, 6), name = "Gene count",
                        breaks = c(30, 60, 90)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(-log[10]~"(p-value)"),
       y = NULL,
       title = "KEGG pathway enrichment — paired comparison") +
  theme_bw(base_size = 11, base_family = plot_font) +
  theme(
    text                = element_text(family = plot_font),
    plot.title          = element_text(family = plot_font, face = "bold",
                                       hjust = 0.5, size = 13),
    axis.title.x        = element_text(face = "bold", size = 11),
    axis.text.y         = element_text(family = plot_font, size = 9.5,
                                       face = "plain", color = "black"),
    axis.text.x         = element_text(family = plot_font, size = 9,
                                       color = "black"),
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_blank(),
    panel.border        = element_rect(color = "black", fill = NA,
                                       linewidth = 0.6),
    strip.placement     = "outside",
    strip.background.y  = element_rect(fill = NA, color = NA),
    strip.text.y.left   = element_blank(),   # ★ category strip 제거
    panel.spacing.y     = unit(0.4, "lines"),
    legend.position     = "top",
    legend.box          = "horizontal",
    legend.title        = element_text(face = "bold", size = 10),
    legend.text         = element_text(size = 9),
    plot.margin         = margin(15, 15, 15, 15)
  )

# Auto-size based on # pathways per panel
n_top <- nrow(panel_top)
n_bot <- nrow(panel_bot)
fig_h <- max(7, (n_top + n_bot) * 0.30 + 2)

ggsave(file.path(out_root, "Fig_1f_KEGG_paired_bars.png"), p_1f,
       width = 9, height = fig_h, dpi = 300, bg = "white")
ggsave(file.path(out_root, "Fig_1f_KEGG_paired_bars.pdf"), p_1f,
       width = 9, height = fig_h, device = cairo_pdf)

# Save aggregated KEGG table
write.csv(kegg_combined,
          file.path(out_root, "tables", "KEGG_2panel_Tier1_Tier2.csv"),
          row.names = FALSE)

cat("Pathways shown — vs Organoid: ", n_top,
    " | vs Macrophage: ", n_bot, "\n", sep = "")
# =====================================================================
# ===== 7. Top 10 inflammation gene selection + LFC heatmap (Fig 1g) =
# =====================================================================
cat("\n===== Generating Figure 1g (Top 10 inflam LFC heatmap) =====\n")

# Inflammation panel — Cytokine + Chemokine + Acute_Inflammation
# =========================================================
# Curated gene panels for BMDM × intestinal organoid co-culture
#  (mouse symbols, MGI convention: Title-case)
# =========================================================

# ---- 1. Inflammation / immune signaling ------------------
cytokines <- c(
  "Tnf","Il1a","Il1b","Il6","Il12a","Il12b","Il18","Il33","Il17a",
  "Il10","Tgfb1","Ifna1","Ifnb1","Ifng",
  "Csf1","Csf2","Csf3","Osm","Lif","Clcf1"
)

cytokine_receptors <- c(
  "Tnfrsf1a","Tnfrsf1b","Il6st","Il1r1","Il10ra",
  "Stat1","Stat3"
)

chemokines <- c(
  "Ccl2","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9",
  "Ccl12","Ccl17","Ccl19","Ccl20","Ccl22",
  "Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcl10","Cxcl11",
  "Cxcl12","Cxcl16"
)

chemokine_receptors <- c(
  "Ccr1","Ccr2","Ccr5","Ccr7","Cxcr3","Cxcr4","Cmklr1"
)

acute_phase_alarmins <- c(
  "Saa3","Lcn2","Ptx3","Cp","Orm1","Orm2",
  "S100a8","S100a9","S100a4","S100a6","Hmgb1",
  "Crp","Hp","Lbp","Mt1","Mt2"
)

vascular_adhesion_cox <- c(
  "Icam1","Vcam1","Sele","Selp","Ptgs2"
)

prr_inflammasome <- c(
  "Tlr2","Tlr4","Tlr7","Tlr9",
  "Nlrp3","Aim2","Casp1","Casp4","Pycard","Gsdmd"
)

nfkb_pathway <- c(
  "Nfkb1","Nfkb2","Rela","Nfkbia","Nfkbid","Tnfaip3","Tnip1"
)

ifn_isg <- c(
  "Irf3","Irf7","Mx1","Mx2","Isg15","Rsad2",
  "Oas1a","Oas2","Ifit1","Ifit2","Ifit3"
)

socs_feedback <- c("Socs1","Socs2","Socs3")

macrophage_polarization <- c(
  "Nos2","Arg1","Cd86","Cd80","Cd163","Mrc1","Mafb"
)

complement <- c(
  "C3","C1qa","C1qb","C1qc","C3ar1","C5ar1"
)

ecm_remodeling <- c("Mmp9","Mmp12","Mmp13")


# ---- 2. Proteostasis / stress response -------------------
hsp_chaperones <- c(
  "Hsp90aa1","Hspa1a","Hspa1b","Hspa5","Hspa8","Hsf1"
)

upr_er_stress <- c(
  "Ddit3","Atf3","Atf4","Atf6","Xbp1"
)


# ---- 3. Redox / Nrf2 axis --------------------------------
nrf2_antioxidant <- c(
  "Nfe2l2","Hmox1","Nqo1","Gclc",
  "Sod1","Sod2","Cat","Gpx1"
)


# ---- 4. DNA damage / cell death --------------------------
dna_damage_p53 <- c("Gadd45a","Trp53")

apoptosis <- c(
  # executioner / initiator caspases
  "Casp3","Casp7","Casp8","Casp9",
  # BH3-only / pro-apoptotic Bcl2
  "Bax","Bak1","Bid","Bbc3","Pmaip1","Bcl2l11",
  # anti-apoptotic
  "Bcl2","Bcl2l1","Mcl1","Birc5","Cflar",
  # apoptosome / mitochondrial release
  "Apaf1","Cycs","Diablo"
)


# ---- 5. Lysosome / autophagy -----------------------------
cathepsins <- c("Ctsa","Ctsb","Ctsd","Ctsl","Ctss","Ctsk")

lysosome_autophagy <- c(
  "Lamp1","Lamp2","Lamp3",
  "Tfeb","Tfe3","Mitf",
  "Hexa","Hexb","Gba","Lipa",
  "Sqstm1","Map1lc3a","Map1lc3b","Atp6v0d2"
)


# =========================================================
# Combined master panel (deduplicated, order preserved)
# =========================================================
gene_panel <- unique(c(
  cytokines, cytokine_receptors, chemokines, chemokine_receptors,
  acute_phase_alarmins, vascular_adhesion_cox, prr_inflammasome,
  nfkb_pathway, ifn_isg, socs_feedback, macrophage_polarization,
  complement, ecm_remodeling,
  hsp_chaperones, upr_er_stress,
  nrf2_antioxidant,
  dna_damage_p53, apoptosis,
  cathepsins, lysosome_autophagy
))

# Optional: keep category labels for heatmap row-splitting / GSVA modules
panel_annot <- data.frame(
  gene = c(cytokines, cytokine_receptors, chemokines, chemokine_receptors,
           acute_phase_alarmins, vascular_adhesion_cox, prr_inflammasome,
           nfkb_pathway, ifn_isg, socs_feedback, macrophage_polarization,
           complement, ecm_remodeling,
           hsp_chaperones, upr_er_stress,
           nrf2_antioxidant,
           dna_damage_p53, apoptosis,
           cathepsins, lysosome_autophagy),
  module = c(
    rep("Cytokines", length(cytokines)),
    rep("Cytokine_Receptors", length(cytokine_receptors)),
    rep("Chemokines", length(chemokines)),
    rep("Chemokine_Receptors", length(chemokine_receptors)),
    rep("Acute_Phase_Alarmins", length(acute_phase_alarmins)),
    rep("Vascular_COX", length(vascular_adhesion_cox)),
    rep("PRR_Inflammasome", length(prr_inflammasome)),
    rep("NFkB", length(nfkb_pathway)),
    rep("IFN_ISG", length(ifn_isg)),
    rep("SOCS", length(socs_feedback)),
    rep("Mac_Polarization", length(macrophage_polarization)),
    rep("Complement", length(complement)),
    rep("ECM", length(ecm_remodeling)),
    rep("HSP", length(hsp_chaperones)),
    rep("UPR", length(upr_er_stress)),
    rep("Nrf2_Antioxidant", length(nrf2_antioxidant)),
    rep("DDR_p53", length(dna_damage_p53)),
    rep("Apoptosis", length(apoptosis)),
    rep("Cathepsins", length(cathepsins)),
    rep("Lysosome_Autophagy", length(lysosome_autophagy))
  ),
  stringsAsFactors = FALSE
)
# Resolve gene-level duplicates by keeping first category
panel_annot <- panel_annot[!duplicated(panel_annot$gene), ]

length(gene_panel)            # final unique gene count
table(panel_annot$module)     # module sizes

# Collect LFC from 4 contrasts
lfc_per_gene <- list()
for (cn in fig1f_contrasts) {
  d <- deg_results[[cn]] %>%
    dplyr::filter(SYMBOL %in% inflam_panel, !is.na(padj))
  lfc_per_gene[[cn]] <- data.frame(
    SYMBOL    = d$SYMBOL,
    LFC       = d$log2FoldChange,
    padj      = d$padj,
    Contrast  = cn,
    stringsAsFactors = FALSE
  )
}

lfc_long <- do.call(rbind, lfc_per_gene)
lfc_long <- lfc_long %>% dplyr::distinct(SYMBOL, Contrast, .keep_all = TRUE)

# Rank by combined criterion: max |LFC| × significance
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

top10_genes <- head(gene_rank$SYMBOL, 10)
cat("\nTop 10 inflammation genes selected (auto):\n")
print(gene_rank %>% dplyr::filter(SYMBOL %in% top10_genes))

# Build LFC matrix: row = top 10 gene, col = 4 contrasts
lfc_mat <- matrix(NA, nrow = length(top10_genes), ncol = length(fig1f_contrasts),
                  dimnames = list(top10_genes, fig1f_contrasts))
padj_mat <- lfc_mat

for (i in seq_along(top10_genes)) {
  for (j in seq_along(fig1f_contrasts)) {
    d <- deg_results[[fig1f_contrasts[j]]] %>%
      dplyr::filter(SYMBOL == top10_genes[i])
    if (nrow(d) > 0) {
      lfc_mat[i, j]  <- d$log2FoldChange[1]
      padj_mat[i, j] <- d$padj[1]
    }
  }
}

# Drop rows with NA in all contrasts
keep_rows <- rowSums(!is.na(lfc_mat)) > 0
lfc_mat  <- lfc_mat[keep_rows, , drop = FALSE]
padj_mat <- padj_mat[keep_rows, , drop = FALSE]

# Significance asterisks
sig_mat <- ifelse(is.na(padj_mat), "",
                  ifelse(padj_mat < 0.001, "***",
                         ifelse(padj_mat < 0.01,  "**",
                                ifelse(padj_mat < 0.05,  "*", ""))))

# Column labels (cleaner)
col_labels <- c(
  T1_OrgBMDM1K_vs_Org = "BMDM 1K\nvs Organoid",
  T1_OrgBMDM5K_vs_Org = "BMDM 5K\nvs Organoid",
  T2_OrgBMDM1K_vs_Mac = "BMDM 1K\nvs Macrophage",
  T2_OrgBMDM5K_vs_Mac = "BMDM 5K\nvs Macrophage"
)
colnames(lfc_mat) <- col_labels[colnames(lfc_mat)]
colnames(sig_mat) <- col_labels[colnames(sig_mat)]

# Column group annotation (vs Organoid / vs Macrophage)
col_group <- data.frame(
  Comparison = factor(c("vs Organoid","vs Organoid",
                        "vs Macrophage","vs Macrophage"),
                      levels = c("vs Organoid","vs Macrophage")),
  row.names = colnames(lfc_mat)
)

ann_color <- list(
  Comparison = c("vs Organoid" = "#F4A6A6",
                 "vs Macrophage" = "#A8D5C2")
)

# Symmetric color scale around 0
lfc_max <- max(abs(lfc_mat), na.rm = TRUE)
lfc_max <- min(lfc_max, 8)  # cap for visibility
breaks  <- seq(-lfc_max, lfc_max, length.out = 101)

pheatmap(lfc_mat,
         color = colorRampPalette(c("#0055FF","#7FB2FF","white","#F8B4BA","#C73E1D"))(100),
         breaks = breaks,
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = col_group,
         annotation_colors = ann_color,
         gaps_col = 2,
         display_numbers = sig_mat,
         number_color = "black", fontsize_number = 11,
         fontsize_row = 11, fontsize_col = 10,
         cellwidth = 50, cellheight = 24,
         angle_col = 0,
         main = "Top 10 inflammation genes — LFC across 4 contrasts\n(* p<.05, ** p<.01, *** p<.001)",
         filename = file.path(out_root, "Fig_1g_TopInflam_LFC_heatmap.png"),
         width = 9, height = max(6, nrow(lfc_mat) * 0.32 + 2.5))

# Save LFC table
lfc_df_out <- data.frame(
  SYMBOL = rownames(lfc_mat),
  lfc_mat,
  check.names = FALSE
)
padj_df_out <- data.frame(
  SYMBOL = rownames(padj_mat),
  padj_mat,
  check.names = FALSE
)
write.csv(lfc_df_out,
          file.path(out_root, "tables", "Top10_inflam_LFC.csv"),
          row.names = FALSE)
write.csv(padj_df_out,
          file.path(out_root, "tables", "Top10_inflam_padj.csv"),
          row.names = FALSE)
write.csv(gene_rank,
          file.path(out_root, "tables", "Inflam_gene_ranking.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 8. Summary ====================================================
# =====================================================================
cat("\n===== Step 4 complete =====\n")
cat("Output root: ", out_root, "/\n\n", sep = "")
cat("  Fig_1e_4way_Venn.{png,pdf}              — 4-way Venn diagram\n")
cat("  Fig_1f_KEGG_paired_bars.{png,pdf}       — Paired KEGG bars (Tier 1 + 2)\n")
cat("  Fig_1g_TopInflam_LFC_heatmap.png        — Top 10 inflammation LFC heatmap\n")
cat("\nTables/:\n")
cat("  Venn_membership_per_gene.csv\n")
cat("  KEGG_paired_Tier1_Tier2.csv\n")
cat("  Inflam_gene_ranking.csv\n")
cat("  Top10_inflam_LFC.csv\n")
cat("  Top10_inflam_padj.csv\n")
cat("\nKEGG blacklist applied: ", length(kegg_blacklist),
    " disease/pathogen pathways removed\n", sep = "")




# =====================================================================
# Fig 1f — KEGG bubble matrix (2 baseline × 4 contrasts)
# Fig 4f와 동일한 규격: 색 막대 strip + 우측 색 범례 + 회색 격자 제거
# =====================================================================

# ===== 1. 폰트 / 크기 / 색상 ==========================================
font_sz <- list(
  base       = 18, title      = 22, subtitle   = 18,
  axis_x     = 16, axis_y     = 18,
  legend_ttl = 18, legend_txt = 16
)
bubble_size_range <- c(4, 14)
bubble_stroke     <- 0.5
strip_size        <- 14     # 막대 두께
strip_margin      <- 5
top_n_per_baseline <- 12    # 베이스라인당 top N pathway

# 2 baseline 정의 + 팔레트 + contrast 매핑
baseline_map_1f <- list(
  "vs Organoid"   = c("T1_OrgBMDM1K_vs_Org",
                                  "T1_OrgBMDM5K_vs_Org"),
  "vs Macrophage" = c("T2_OrgBMDM1K_vs_Mac",
                                  "T2_OrgBMDM5K_vs_Mac")
)
baseline_levels_1f  <- names(baseline_map_1f)
baseline_palette_1f <- c(
  "vs Organoid"   = "#1F4E79",   # dark blue
  "vs Macrophage" = "#A61C00"    # dark red
)

contrast_label_map_1f <- c(
  T1_OrgBMDM1K_vs_Org = "Org+BMDM 1K",
  T1_OrgBMDM5K_vs_Org = "Org+BMDM 5K",
  T2_OrgBMDM1K_vs_Mac = "Org+BMDM 1K",
  T2_OrgBMDM5K_vs_Mac = "Org+BMDM 5K"
)

# ===== 2. Long-format build (Fig 4f와 동일 로직) =====================
contrast_baseline_map_1f <- list()
for (bl in names(baseline_map_1f)) {
  for (cn in baseline_map_1f[[bl]]) contrast_baseline_map_1f[[cn]] <- bl
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

kegg_bubble_long_1f <- build_kegg_long_bubble(kegg_list,
                                              contrast_baseline_map_1f)
kegg_bubble_long_1f$Baseline <- factor(kegg_bubble_long_1f$Baseline,
                                       levels = baseline_levels_1f)
kegg_bubble_long_1f$logp <- -log10(kegg_bubble_long_1f$pvalue)

# ===== 3. Pathway 선정 — 베이스라인별 top N union =====================
selected_paths_1f <- kegg_bubble_long_1f %>%
  dplyr::group_by(Baseline, Description) %>%
  dplyr::summarise(best_logp = max(logp), .groups = "drop") %>%
  dplyr::group_by(Baseline) %>%
  dplyr::arrange(dplyr::desc(best_logp), .by_group = TRUE) %>%
  dplyr::slice_head(n = top_n_per_baseline) %>%
  dplyr::ungroup()

# 각 pathway는 best_logp 가장 큰 baseline에 primary 배정
pathway_primary_1f <- selected_paths_1f %>%
  dplyr::group_by(Description) %>%
  dplyr::slice_max(best_logp, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Baseline, dplyr::desc(best_logp))

kegg_bubble_sub_1f <- kegg_bubble_long_1f %>%
  dplyr::filter(Description %in% pathway_primary_1f$Description)

kegg_bubble_sub_1f$Description <- factor(
  kegg_bubble_sub_1f$Description,
  levels = rev(pathway_primary_1f$Description)
)

contrast_order_1f <- unlist(baseline_map_1f, use.names = FALSE)
kegg_bubble_sub_1f$Contrast <- factor(kegg_bubble_sub_1f$Contrast,
                                      levels = contrast_order_1f)

# ===== 4. 플롯 빌드 ===================================================
baseline_legend_data_1f <- data.frame(
  Baseline = factor(baseline_levels_1f, levels = baseline_levels_1f)
)
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)

p_1f <- ggplot(kegg_bubble_sub_1f,
               aes(x = Contrast, y = Description)) +
  geom_point(aes(size = Count, fill = logp),
             shape = 21, color = "black",
             stroke = bubble_stroke, alpha = 0.95) +
  geom_point(data = baseline_legend_data_1f,
             aes(color = Baseline),
             x = NA, y = NA, size = 0, na.rm = TRUE,
             inherit.aes = FALSE, show.legend = TRUE)

if (has_ggh4x) {
  strip_fills_1f <- unname(baseline_palette_1f[baseline_levels_1f])
  
  p_1f <- p_1f +
    ggh4x::facet_grid2(
      . ~ Baseline,
      scales = "free_x", space = "free_x",
      strip = ggh4x::strip_themed(
        background_x = lapply(strip_fills_1f, function(col)
          element_rect(fill = col, color = "black", linewidth = 0.6)),
        text_x = lapply(strip_fills_1f, function(col)
          element_text(color  = NA,             # 글자 투명
                       size   = strip_size,
                       margin = margin(t = strip_margin,
                                       b = strip_margin)))
      )
    )
} else {
  p_1f <- p_1f +
    facet_grid(. ~ Baseline, scales = "free_x", space = "free_x")
}

p_1f <- p_1f +
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
    values = baseline_palette_1f,
    name   = "Baseline",
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(size = 7, shape = 15, alpha = 1),
      order = 3
    )
  ) +
  scale_x_discrete(labels = contrast_label_map_1f) +
  labs(x = NULL, y = NULL,
       title = "KEGG pathway enrichment — bubble matrix (Tier 1 + Tier 2)",
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
    panel.grid       = element_blank(),                    # 회색 격자 제거
    panel.border     = element_rect(color = "black", fill = NA,
                                    linewidth = 0.7),
    panel.spacing.x  = unit(0.5, "lines"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = font_sz$legend_ttl),
    legend.text      = element_text(size = font_sz$legend_txt),
    legend.box       = "vertical",
    legend.spacing.y = unit(0.4, "cm"),
    plot.margin      = margin(15, 15, 15, 15)
  )

# ===== 5. 출력 ========================================================
while (!is.null(dev.list())) dev.off()
print(p_1f)

# ===== 6. 저장 (PNG / PDF / TIFF) =====================================
n_path_1f <- length(unique(kegg_bubble_sub_1f$Description))
fig_w_1f  <- 11
fig_h_1f  <- max(8, n_path_1f * 0.40 + 4)

# PNG
ggsave(file.path(out_root, "Fig_1f_KEGG_bubble_matrix.png"),
       p_1f, width = fig_w_1f, height = fig_h_1f,
       dpi = 300, bg = "white", limitsize = FALSE)

# PDF (vector)
ggsave(file.path(out_root, "Fig_1f_KEGG_bubble_matrix.pdf"),
       p_1f, width = fig_w_1f, height = fig_h_1f,
       device = cairo_pdf, limitsize = FALSE)

# TIFF (출판용 — 300 dpi + LZW 압축 + white background)
ggsave(file.path(out_root, "Fig_1f_KEGG_bubble_matrix.tiff"),
       p_1f, width = fig_w_1f, height = fig_h_1f,
       device = "tiff", dpi = 300, bg = "white",
       compression = "lzw", limitsize = FALSE)

write.csv(kegg_bubble_sub_1f,
          file.path(out_root, "tables", "KEGG_bubble_matrix_Tier1_Tier2.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 7. Top 10 inflammation gene LFC heatmap (Fig 1g) — print 가능
# =====================================================================
cat("\n===== Generating Figure 1g (Top 10 inflam LFC heatmap) =====\n")

# (gene_panel, panel_annot 정의 부분은 기존 그대로 — 생략)
# ...

# ===== 7-1. LFC 수집 (★ inflam_panel → gene_panel 로 수정) ==========
lfc_per_gene <- list()
for (cn in fig1f_contrasts) {
  d <- deg_results[[cn]] %>%
    dplyr::filter(SYMBOL %in% gene_panel, !is.na(padj))   # ← gene_panel
  lfc_per_gene[[cn]] <- data.frame(
    SYMBOL    = d$SYMBOL,
    LFC       = d$log2FoldChange,
    padj      = d$padj,
    Contrast  = cn,
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

top10_genes <- head(gene_rank$SYMBOL, 10)
cat("\nTop 10 inflammation genes selected (auto):\n")
print(gene_rank %>% dplyr::filter(SYMBOL %in% top10_genes))

# ===== 7-2. LFC / padj matrix (기존 그대로) =========================
lfc_mat <- matrix(NA, nrow = length(top10_genes), ncol = length(fig1f_contrasts),
                  dimnames = list(top10_genes, fig1f_contrasts))
padj_mat <- lfc_mat
for (i in seq_along(top10_genes)) {
  for (j in seq_along(fig1f_contrasts)) {
    d <- deg_results[[fig1f_contrasts[j]]] %>%
      dplyr::filter(SYMBOL == top10_genes[i])
    if (nrow(d) > 0) {
      lfc_mat[i, j]  <- d$log2FoldChange[1]
      padj_mat[i, j] <- d$padj[1]
    }
  }
}
keep_rows <- rowSums(!is.na(lfc_mat)) > 0
lfc_mat  <- lfc_mat[keep_rows, , drop = FALSE]
padj_mat <- padj_mat[keep_rows, , drop = FALSE]

sig_mat <- ifelse(is.na(padj_mat), "",
                  ifelse(padj_mat < 0.001, "***",
                         ifelse(padj_mat < 0.01,  "**",
                                ifelse(padj_mat < 0.05,  "*", ""))))

# ★ 컬럼명 단순화 (vs ... 제거) — 중복 방지 위해 dose만 표시
col_labels <- c(
  T1_OrgBMDM1K_vs_Org = "Org+BMDM 1K ",       # 끝에 공백 → 다른 컬럼명과 구별
  T1_OrgBMDM5K_vs_Org = "Org+BMDM 5K ",
  T2_OrgBMDM1K_vs_Mac = "Org+BMDM 1K",
  T2_OrgBMDM5K_vs_Mac = "Org+BMDM 5K"
)
colnames(lfc_mat) <- col_labels[colnames(lfc_mat)]
colnames(sig_mat) <- col_labels[colnames(sig_mat)]

col_group <- data.frame(
  Panel = factor(c("vs Organoid","vs Organoid",            # ★ Comparison → Panel
                   "vs Macrophage","vs Macrophage"),
                 levels = c("vs Organoid","vs Macrophage")),
  row.names = colnames(lfc_mat)
)
ann_color <- list(
  Panel = c("vs Organoid"   = "#1F4E79",                   # ★ list 키도 Panel
            "vs Macrophage" = "#A61C00")
)

lfc_max <- max(abs(lfc_mat), na.rm = TRUE)
lfc_max <- min(lfc_max, 8)
breaks  <- seq(-lfc_max, lfc_max, length.out = 101)
# ===== Legend label 위치 (color bar 위에 "Log2FC" 라벨용) ===========
legend_breaks_1g <- pretty(c(-lfc_max, lfc_max), n = 5)
legend_breaks_1g <- legend_breaks_1g[legend_breaks_1g >= -lfc_max &
                                       legend_breaks_1g <=  lfc_max]
legend_labels_1g <- c(as.character(legend_breaks_1g), "Log2FC")
legend_breaks_1g <- c(legend_breaks_1g, lfc_max)


# ===== 7-3. 폰트 / 셀 크기 — 튜닝 포인트 =============================
hm_params_1g <- list(
  fontsize        = 20,    # ★ 전체 base (title, legend)
  fontsize_row    = 18,    # ★ gene symbol
  fontsize_col    = 18,    # ★ column label
  fontsize_number = 16,    # ★ display_numbers (asterisks)
  cellwidth       = 40,    # ★ 셀 너비
  cellheight      = 30     # ★ 셀 높이
)

# ===== 7-4. Heatmap with ComplexHeatmap (기존 그대로) ===============
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
}
library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun_1g <- colorRamp2(
  c(-lfc_max, -lfc_max/2, 0, lfc_max/2, lfc_max),
  c("#0055FF", "#7FB2FF", "white", "#F8B4BA", "#C73E1D")
)

top_anno_1g <- HeatmapAnnotation(
  Panel = col_group$Panel,
  col   = list(Panel = ann_color$Panel),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Panel = list(title = "Panel",
                 title_gp  = gpar(fontsize = hm_params_1g$fontsize, fontface = "bold"),
                 labels_gp = gpar(fontsize = hm_params_1g$fontsize_col))
  )
)

ht_1g <- Heatmap(
  lfc_mat,
  name              = "Log2(FC)",
  col               = col_fun_1g,
  cluster_rows      = FALSE, cluster_columns = FALSE,
  show_row_names    = TRUE,
  show_column_names = TRUE,
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = hm_params_1g$fontsize_row),
  column_names_gp   = gpar(fontsize = hm_params_1g$fontsize_col),
  column_names_rot  = 45,
  top_annotation    = top_anno_1g,
  column_split      = col_group$Panel,
  column_title      = NULL,
  border            = TRUE,
  rect_gp           = gpar(col = "black", lwd = 0.5),
  width             = unit(hm_params_1g$cellwidth  * ncol(lfc_mat), "pt"),
  height            = unit(hm_params_1g$cellheight * nrow(lfc_mat), "pt"),
  heatmap_legend_param = list(
    title          = "Log2(FC)",
    title_gp       = gpar(fontsize = 14, fontface = "bold"),
    labels_gp      = gpar(fontsize = hm_params_1g$fontsize_col),
    title_position = "leftcenter-rot",
    legend_height  = unit(5, "cm"),
    grid_width     = unit(0.5, "cm")
  )
)

# ===== 7-5. draw 함수 (저장 시 재사용) ==============================
draw_ht_1g <- function() {
  draw(ht_1g,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       column_title    = "Top 10 inflammation genes — LFC across 4 contrasts",
       column_title_gp = gpar(fontsize = hm_params_1g$fontsize + 2,
                              fontface = "bold"))
}

# ===== 7-6. 화면 출력 ===============================================
while (!is.null(dev.list())) dev.off()
draw_ht_1g()

# ===== 7-7. 저장 (PNG / PDF / TIFF) — ComplexHeatmap 방식 ===========
fig_h_1g <- max(7, nrow(lfc_mat) * 0.45 + 3)
fig_w_1g <- 10

# PNG
png(file.path(out_root, "Fig_1g_TopInflam_LFC_heatmap.png"),
    width = fig_w_1g, height = fig_h_1g,
    units = "in", res = 300, bg = "white")
draw_ht_1g()
dev.off()

# PDF (vector)
pdf(file.path(out_root, "Fig_1g_TopInflam_LFC_heatmap.pdf"),
    width = fig_w_1g, height = fig_h_1g, bg = "white")
draw_ht_1g()
dev.off()

# TIFF (출판용, 300 dpi + LZW)
tiff(file.path(out_root, "Fig_1g_TopInflam_LFC_heatmap.tiff"),
     width = fig_w_1g, height = fig_h_1g,
     units = "in", res = 300, bg = "white",
     compression = "lzw")
draw_ht_1g()
dev.off()

cat("\nFig 1g saved: PNG / PDF / TIFF\n")

# ===== 7-8. Tables export ===========================================
lfc_df_out  <- data.frame(SYMBOL = rownames(lfc_mat), lfc_mat,
                          check.names = FALSE)
padj_df_out <- data.frame(SYMBOL = rownames(padj_mat), padj_mat,
                          check.names = FALSE)
write.csv(lfc_df_out,
          file.path(out_root, "tables", "Top10_inflam_LFC.csv"),
          row.names = FALSE)
write.csv(padj_df_out,
          file.path(out_root, "tables", "Top10_inflam_padj.csv"),
          row.names = FALSE)
write.csv(gene_rank,
          file.path(out_root, "tables", "Inflam_gene_ranking.csv"),
          row.names = FALSE)