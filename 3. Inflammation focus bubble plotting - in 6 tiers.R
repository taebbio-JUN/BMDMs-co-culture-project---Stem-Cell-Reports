# =====================================================================
# Step 3: INFLAMMATION SUB-PANEL ACTIVATION ANALYSIS (6-Tier)
# Project: HN00273522 — BMDM × Macrophage × TNFα × Intestinal Organoid
#
# Super-categories (2):
#   ① Inflammation                  : 7 sub-panels
#   ② Inflammatory_Apoptosis_Resolution : 2 sub-panels
#
# Visualizations:
#   • LFC-based Heatmap (per tier, separate folder)
#   • Sub-panel colored Volcano (per contrast)
#   • KEGG enrichment (per contrast)
#   • GSVA supplementary
# =====================================================================

# ===== 0. Packages ===================================================
suppressPackageStartupMessages({
  library(DESeq2);   library(ggplot2); library(ggrepel)
  library(dplyr);    library(tidyr);   library(tibble)
  library(pheatmap); library(RColorBrewer); library(patchwork)
  library(clusterProfiler); library(enrichplot); library(DOSE)
  library(org.Mm.eg.db);   library(AnnotationDbi)
  library(GSVA)     # for GSVA supplementary
})

# ===== 1. Configuration ==============================================
out_root <- "Step3_Inflammation_Activation"

padj_cutoff <- 0.05
lfc_cutoff  <- 1.0
n_each_side <- 50
p_cap       <- 300
plot_font   <- "Helvetica"

# ===== 2. Load Step 1 outputs ========================================
dds            <- readRDS("dds_full_8groups.rds")
vsd            <- readRDS("vsd_full_8groups.rds")
deg_results    <- readRDS("deg_results_full_8groups.rds")
anno           <- readRDS("transcript_annotation.rds")
contrast_tiers <- readRDS("contrast_tiers.rds")

tier_names <- unique(contrast_tiers$tier)

# 8-group color scheme
group_colors_full <- c(
  Organoid              = "#A6CEE3",
  Organoid_TNFa         = "#1F78B4",
  Organoid_BMDM_1K      = "#B2DF8A",
  Organoid_TNFa_BMDM_1K = "#33A02C",
  Organoid_BMDM_5K      = "#FDBF6F",
  Organoid_TNFa_BMDM_5K = "#FF7F00",
  Macrophage            = "#FB9A99",
  Macrophage_TNFa       = "#E31A1C"
)

# ===== 3. Folder structure ===========================================
dir.create(out_root, showWarnings = FALSE)
for (tn in tier_names) {
  dir.create(file.path(out_root, tn),                showWarnings = FALSE)
  dir.create(file.path(out_root, tn, "heatmap"),     showWarnings = FALSE)
  dir.create(file.path(out_root, tn, "volcano"),     showWarnings = FALSE)
  dir.create(file.path(out_root, tn, "enrich"),      showWarnings = FALSE)
}
dir.create(file.path(out_root, "tables"),            showWarnings = FALSE)
dir.create(file.path(out_root, "GSVA_supplementary"), showWarnings = FALSE)
dir.create(file.path(out_root, "GSVA_supplementary", "heatmap"), showWarnings = FALSE)

# =====================================================================
# ===== 4. Sub-panel definitions (9 panels, 2 super-categories) =======
# =====================================================================

sub_panels <- list(
  
  # ===== Super-category ① INFLAMMATION =====
  
  Acute_Inflammation = c(
    # Acute-phase proteins (liver/local secreted)
    "Saa3","Crp","Hp","Lcn2","Lbp","Ptx3","Cp","Orm1","Orm2",
    # Alarmins / DAMP
    "S100a4","S100a6","S100a8","S100a9","Hmgb1",
    # Early-response cytokines
    "Il6","Tnf","Il1b","Il1a",
    # Early chemokines (first-wave)
    "Cxcl1","Cxcl2","Ccl2",
    # Vascular adhesion (acute leukocyte extravasation)
    "Icam1","Vcam1","Sele","Selp",
    # COX-2 / prostaglandin (acute)
    "Ptgs2"
  ),
  
  TLR_NFkB_Signaling = c(
    # TLR family
    "Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tlr13",
    # TLR adapters / co-receptors
    "Cd14","Ly96","Cd180","Lbp","Myd88","Ticam1","Ticam2","Tirap",
    "Irak1","Irak3","Irak4",
    # C-type lectin (PRR)
    "Card9","Clec4a1","Clec4a2","Clec4a3","Clec4e","Clec4n","Clec5a",
    "Clec7a","Clec12a",
    # NF-κB family
    "Nfkb1","Nfkb2","Rela","Relb","Rel","Nfkbib",
    # NF-κB inhibitors / regulators (feedback)
    "Nfkbia","Nfkbid","Nfkbie","Ikbkg","Ikbkb",
    # TRAF / receptor adapters
    "Traf1","Traf2","Traf3","Traf6","Traf3ip3","Tifab",
    "Tnfrsf1a","Tnfrsf1b"
  ),
  
  Cytokine_Storm_Core = c(
    # Pro-inflammatory ligands (core)
    "Tnf","Il1a","Il1b","Il6","Il12a","Il12b","Il18","Il33","Il17a","Il17f",
    # Regulatory ligands (homeostatic balance — keep here)
    "Il10","Tgfb1","Tgfb2","Tgfb3",
    # Type I/II IFN ligands
    "Ifna1","Ifnb1","Ifng","Ifnl2","Ifnl3",
    # Myeloid CSF / gp130 cytokines
    "Csf1","Csf2","Csf3","Clcf1","Osm","Lif",
    # Cytokine receptors
    "Tnfrsf1a","Tnfrsf1b","Il6st","Il6ra","Il1r1","Il1r2","Il10ra","Il10rb",
    "Il12rb1","Il12rb2","Il18r1","Il18rap","Il23r","Il1rl1","Il1rl2",
    "Il7r","Il2rg",
    "Csf1r","Csf2ra","Csf2rb","Csf2rb2","Csf3r",
    "Ifnar1","Ifnar2","Ifngr1","Ifngr2",
    # Canonical JAK-STAT signaling
    "Stat1","Stat2","Stat3","Stat4","Stat5a","Stat5b","Stat6",
    "Jak1","Jak2","Jak3","Tyk2"
  ),
  
  Chemokine_Recruitment = c(
    # CC family
    "Ccl2","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9",
    "Ccl11","Ccl12","Ccl17","Ccl19","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25",
    # CXC family
    "Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcl10","Cxcl11","Cxcl12",
    "Cxcl13","Cxcl14","Cxcl16",
    # CX3C / others
    "Cx3cl1","Xcl1","Ppbp","Pf4",
    # Receptors
    "Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr9","Ccrl2",
    "Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Cx3cr1","Cmklr1"
  ),
  
  Interferon_Response = c(
    # IFN regulatory factors
    "Irf1","Irf3","Irf5","Irf7","Irf8","Irf9",
    # Cytosolic nucleic acid sensors
    "Ddx58","Ifih1","Dhx58","Mavs","Sting1","Tmem173","Cgas","Zbp1",
    # Classical ISGs
    "Ifit1","Ifit2","Ifit3","Ifi27","Ifi35","Ifi44","Ifi47","Rtp4",
    "Mx1","Mx2","Isg15","Rsad2","Usp18","Ly6e","Cmpk2","Herc6",
    "Samd9l","Xaf1",
    "Oas1a","Oas1g","Oas2","Oas3","Oasl1","Oasl2",
    # IFN-γ-induced GTPases
    "Gbp2","Gbp3","Gbp4","Gbp7","Gbp9","Igtp","Iigp1","Tgtp1","Tgtp2",
    # ISG15 conjugation system
    "Uba7","Ube2l6","Bst2",
    # Slfn family
    "Slfn1","Slfn2","Slfn4","Slfn5","Slfn8",
    # AIM2-like / IFI20x family
    "Ifi202b","Ifi203","Ifi204","Ifi205","Mnda","Mndal","Pyhin1",
    # IFN-induced TRIM ubiquitin ligases
    "Trim15","Trim21","Trim25","Trim30a","Trim30c","Trim30d",
    "Trim34a","Trim34b","Trim40","Trim47"
  ),
  
  Macrophage_Induced_Inflammation = c(
    # M1 polarization markers
    "Nos2","Cd86","Cd80","Cd40",
    # Macrophage identity
    "Cd68","Emr1","Adgre1","Csf1r","Aif1","Mafb","Mrc1","Cd14","Marco","Msr1",
    # M1 effector enzymes
    "Mmp9","Mmp12","Mmp13","Mmp3","Mmp8","Mmp27","Ptgs2","Tbxas1",
    # NADPH oxidase (ROS production)
    "Cybb","Cyba","Ncf1","Ncf2","Ncf4","Nox1","Duoxa2","Hvcn1",
    # iNOS / NO pathway
    "Gch1","Ddah1","Ddah2","Arg2",
    # M1-associated effectors
    "Slpi","Cd274","Havcr2","Slc11a1",
    # Macrophage-specific signaling
    "Irg1","Acod1","Ido1","Ido2"
  ),
  
  Lysosomal_Phagocytic = c(
    # Cathepsins (lysosomal proteases)
    "Ctsa","Ctsb","Ctsd","Ctsg","Ctsk","Ctsl","Ctss","Ctsz",
    # LAMPs
    "Lamp1","Lamp2","Lamp3","Laptm5",
    # Lysosomal biogenesis (TFEB axis)
    "Tfeb","Tfe3","Mitf",
    # Lysosomal hydrolases
    "Hexa","Hexb","Gba","Lipa","Lgmn","Acp2","Napsa",
    # V-ATPase
    "Atp6v0d1","Atp6v0d2","Atp6v1h",
    # Autophagy-lysosome interface
    "Sqstm1","Map1lc3a","Map1lc3b","Becn1","Atg5","Atg7","Atg12","Gabarapl1",
    # Phagosome maturation
    "Lyz1","Lyz2","Mpo","Elane","Prtn3","Ngp","Srgn",
    # Fc receptors (phagocytic uptake)
    "Fcgr1","Fcgr2b","Fcgr3","Fcgr4","Fcer1g",
    # Scavenger receptors
    "Cd36","Olr1","Scarf1","Stab1","Stab2"
  ),
  
  # ===== Super-category ② INFLAMMATORY_APOPTOSIS_RESOLUTION =====
  
  Inflammatory_Cell_Death = c(
    # Apoptosis — Tumor suppressor + caspases
    "Trp53","Casp3","Casp6","Casp7","Casp8","Casp9","Casp2","Casp10",
    # Apoptosis — Pro-apoptotic Bcl-2 family
    "Bax","Bak1","Bid","Bbc3","Pmaip1","Bcl2l11","Bik","Bmf","Hrk",
    # Apoptosis — Anti-apoptotic Bcl-2 family
    "Bcl2","Bcl2l1","Mcl1","Bcl2a1a","Bcl2a1b","Bcl2a1d","Bcl2l2",
    # Apoptosis — Mitochondrial / adapter
    "Apaf1","Cycs","Diablo","Aifm1","Endog",
    # IAP / FLIP / receptor adapters
    "Birc2","Birc3","Birc5","Xiap","Cflar","Fadd","Fas","Fasl","Tnfrsf10b",
    # Necroptosis
    "Ripk1","Ripk3","Mlkl","Zbp1","Faim",
    # Pyroptosis
    "Casp1","Casp4","Gsdmd","Gsdme","Pycard","Aim2",
    "Nlrp1a","Nlrp3","Nlrp6","Nlrc4","Mefv",
    # Ferroptosis — core
    "Gpx4","Acsl4","Lpcat3","Alox15","Aifm2","Slc7a11","Slc3a2",
    # Ferroptosis — iron handling
    "Tfrc","Fth1","Ftl1","Slc40a1","Ncoa4","Pcbp1","Pcbp2",
    "Steap4","Trpm2","Trpm7"
  ),
  
  Resolution_Anti_inflammatory = c(
    # Anti-inflammatory cytokines & receptors
    "Il10","Il10ra","Il10rb","Tgfb1","Tgfb2","Tgfb3",
    "Tgfbr1","Tgfbr2","Tgfbr3","Il1rn","Il13","Il4",
    # SOCS family (JAK-STAT negative feedback)
    "Socs1","Socs2","Socs3","Socs4","Socs5","Socs6","Socs7","Cish",
    # A20 / ubiquitin-editing (NF-κB termination)
    "Tnfaip3","Tnip1","Tnip2","Tnip3","Tnfaip8","Tnfaip8l2","Cyld",
    # Anti-inflammatory transcription factors
    "Klf2","Klf4","Foxp3","Nr4a1","Nr4a2","Nr4a3",
    # Resolution lipid mediators (enzymes)
    "Alox5","Ptgr1","Ptgr2",
    # Pro-resolving / efferocytosis
    "Mertk","Axl","Tyro3","Gas6","Anxa1","Lxn",
    # IL-1 / TLR pathway feedback
    "Irak3","Sigirr","Il1r2",
    # Heme oxygenase (anti-inflammatory)
    "Hmox1","Mt1","Mt2",
    # M2 macrophage / alternative activation markers
    "Arg1","Mrc1","Chil3","Retnla","Cd163","Mgl2"
  )
)

# Super-category map
super_cat_map <- c(
  Acute_Inflammation              = "Inflammation",
  TLR_NFkB_Signaling              = "Inflammation",
  Cytokine_Storm_Core             = "Inflammation",
  Chemokine_Recruitment           = "Inflammation",
  Interferon_Response             = "Inflammation",
  Macrophage_Induced_Inflammation = "Inflammation",
  Lysosomal_Phagocytic            = "Inflammation",
  Inflammatory_Cell_Death         = "Inflammatory_Apoptosis_Resolution",
  Resolution_Anti_inflammatory    = "Inflammatory_Apoptosis_Resolution"
)
sub_panel_order <- names(sub_panels)
super_cat_order <- c("Inflammation", "Inflammatory_Apoptosis_Resolution")

# Build long table — gene × panel (first-listed panel priority)
gene_panel_long <- do.call(rbind, lapply(sub_panel_order, function(p) {
  data.frame(SubPanel = p, SuperCat = super_cat_map[[p]],
             SYMBOL = sub_panels[[p]], stringsAsFactors = FALSE)
}))
gene_panel_unique <- gene_panel_long %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# Save panel definitions
panel_def_df <- gene_panel_long %>%
  dplyr::group_by(SuperCat, SubPanel) %>%
  dplyr::summarise(n_genes = dplyr::n(),
                   genes = paste(SYMBOL, collapse = ", "),
                   .groups = "drop")
write.csv(panel_def_df,
          file.path(out_root, "tables", "panel_definitions.csv"),
          row.names = FALSE)

cat("\n===== Sub-panel summary =====\n")
panel_summary <- gene_panel_long %>%
  dplyr::group_by(SuperCat, SubPanel) %>%
  dplyr::summarise(n_total = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(n_unique_after_dedup = sapply(SubPanel, function(p) {
    sum(gene_panel_unique$SubPanel == p)
  }))
print(panel_summary)

# Sub-panel color palette
sub_panel_palette <- c(
  Acute_Inflammation              = "#D62728",
  TLR_NFkB_Signaling              = "#E08214",
  Cytokine_Storm_Core             = "#F46D43",
  Chemokine_Recruitment           = "#FDAE61",
  Interferon_Response             = "#762A83",
  Macrophage_Induced_Inflammation = "#9970AB",
  Lysosomal_Phagocytic            = "#C2A5CF",
  Inflammatory_Cell_Death         = "#1B7837",
  Resolution_Anti_inflammatory    = "#5AAE61"
)
super_cat_palette <- c(
  Inflammation                       = "#C73E1D",
  Inflammatory_Apoptosis_Resolution  = "#1B7837"
)

# =====================================================================
# ===== 5. Transcript → Gene collapse =================================
# =====================================================================
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

# Annotate with sub-panel
build_annotated_deg <- function(deg_df, contrast_name) {
  deg_df %>%
    dplyr::filter(!is.na(SYMBOL), !is.na(log2FoldChange), !is.na(padj)) %>%
    dplyr::mutate(
      contrast_full  = contrast_name,
      sig            = padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff,
      neg_log10_padj = -log10(pmax(padj, 10^(-p_cap))),
      dot_class = dplyr::case_when(
        is.na(padj) | padj >= padj_cutoff                  ~ "ns",
        padj < padj_cutoff & log2FoldChange >=  lfc_cutoff ~ "sig_up",
        padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "sig_down",
        padj < padj_cutoff & log2FoldChange > 0            ~ "sub_up",
        padj < padj_cutoff & log2FoldChange < 0            ~ "sub_down",
        TRUE                                               ~ "ns"
      )
    ) %>%
    dplyr::left_join(gene_panel_unique, by = "SYMBOL") %>%
    dplyr::mutate(
      SubPanel = ifelse(is.na(SubPanel), "Other", SubPanel),
      SuperCat = ifelse(is.na(SuperCat), "Other", SuperCat)
    )
}
deg_annotated <- lapply(names(deg_gene), function(nm)
  build_annotated_deg(deg_gene[[nm]], nm))
names(deg_annotated) <- names(deg_gene)

# =====================================================================
# ===== 6. Build per-tier LFC matrix (genes × contrasts) ==============
# =====================================================================
# Returns: list of LFC matrices, one per tier
#         row = gene (SYMBOL), col = contrast (within tier)
#         NA / low-expression genes are REMOVED (not shown as white)

build_tier_lfc_matrix <- function(panel_genes, tier_name, deg_annotated_list,
                                  contrast_tiers) {
  contrasts_in_tier <- contrast_tiers$contrast[contrast_tiers$tier == tier_name]
  if (length(contrasts_in_tier) == 0) return(NULL)
  
  # Build LFC matrix
  lfc_mat <- sapply(contrasts_in_tier, function(cn) {
    d <- deg_annotated_list[[cn]]
    lfc_vec <- d$log2FoldChange[match(panel_genes, d$SYMBOL)]
    names(lfc_vec) <- panel_genes
    lfc_vec
  })
  rownames(lfc_mat) <- panel_genes
  
  # Build padj matrix (for significance asterisks)
  padj_mat <- sapply(contrasts_in_tier, function(cn) {
    d <- deg_annotated_list[[cn]]
    pj <- d$padj[match(panel_genes, d$SYMBOL)]
    names(pj) <- panel_genes
    pj
  })
  rownames(padj_mat) <- panel_genes
  
  # REMOVE genes that have NA in all contrasts (not present in filtered dds)
  keep <- apply(lfc_mat, 1, function(x) !all(is.na(x)))
  lfc_mat  <- lfc_mat[keep, , drop = FALSE]
  padj_mat <- padj_mat[keep, , drop = FALSE]
  
  # If a gene has NA in *some* (but not all) contrasts, also drop it
  # to ensure clean heatmap (no white cells in middle of rows)
  complete_rows <- complete.cases(lfc_mat)
  lfc_mat  <- lfc_mat[complete_rows, , drop = FALSE]
  padj_mat <- padj_mat[complete_rows, , drop = FALSE]
  
  list(lfc = lfc_mat, padj = padj_mat,
       contrasts = contrasts_in_tier,
       genes = rownames(lfc_mat))
}

# =====================================================================
# ===== 7. Heatmap: per-tier, per-super-category ======================
# =====================================================================
draw_lfc_heatmap <- function(tier_name, super_cat, sub_panels_in_cat,
                             deg_annotated_list, contrast_tiers,
                             file_suffix) {
  
  # Combine all sub-panel genes within super-category
  panel_genes_all <- unique(unlist(sub_panels[sub_panels_in_cat]))
  
  mat_list <- build_tier_lfc_matrix(panel_genes_all, tier_name,
                                    deg_annotated_list, contrast_tiers)
  if (is.null(mat_list)) return(invisible(NULL))
  lfc_mat  <- mat_list$lfc
  padj_mat <- mat_list$padj
  
  if (nrow(lfc_mat) < 2) {
    cat("  Skipping ", tier_name, " × ", super_cat, ": <2 genes available\n",
        sep = "")
    return(invisible(NULL))
  }
  
  # Assign sub-panel to each gene (first-listed priority)
  gene_to_subpanel <- setNames(gene_panel_unique$SubPanel,
                               gene_panel_unique$SYMBOL)
  sub_vec <- gene_to_subpanel[rownames(lfc_mat)]
  sub_vec_f <- factor(sub_vec, levels = intersect(sub_panels_in_cat,
                                                  unique(sub_vec)))
  
  # Order rows by sub-panel
  ord <- order(sub_vec_f)
  lfc_mat  <- lfc_mat[ord, , drop = FALSE]
  padj_mat <- padj_mat[ord, , drop = FALSE]
  sub_vec_f <- sub_vec_f[ord]
  
  row_groups <- data.frame(SubPanel = sub_vec_f, row.names = rownames(lfc_mat))
  
  # Gap rows between sub-panels
  gaps_row <- cumsum(table(droplevels(sub_vec_f)))
  gaps_row <- as.numeric(gaps_row[-length(gaps_row)])
  
  # Annotation colors
  panels_present <- levels(droplevels(sub_vec_f))
  ann_colors <- list(SubPanel = sub_panel_palette[panels_present])
  
  # Significance asterisks (* < 0.05, ** < 0.01, *** < 0.001)
  sig_mat <- ifelse(is.na(padj_mat), "",
                    ifelse(padj_mat < 0.001, "***",
                           ifelse(padj_mat < 0.01,  "**",
                                  ifelse(padj_mat < 0.05,  "*", ""))))
  
  # Color breaks — symmetric around 0
  lfc_max <- max(abs(lfc_mat), na.rm = TRUE)
  lfc_max <- min(lfc_max, 6)  # cap at ±6 for visibility
  breaks  <- seq(-lfc_max, lfc_max, length.out = 101)
  
  # Get control group for title
  ctrl_group <- unique(contrast_tiers$control_group[contrast_tiers$tier == tier_name])
  
  # Title
  ttl <- paste0("[", super_cat, "] LFC Heatmap — ", tier_name,
                "\n(Control = ", ctrl_group, ",  * p<.05, ** p<.01, *** p<.001)")
  
  # Auto-size height
  h <- max(7, nrow(lfc_mat) * 0.13 + 3)
  w <- max(7, ncol(lfc_mat) * 1.3 + 5)
  
  pheatmap(lfc_mat,
           color = colorRampPalette(c("#0055FF","#7FB2FF","white","#F8B4BA","#C73E1D"))(100),
           breaks = breaks,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = row_groups,
           annotation_colors = ann_colors,
           gaps_row = gaps_row,
           display_numbers = sig_mat,
           number_color = "black",
           fontsize_number = 7,
           fontsize_row = 7, fontsize_col = 9,
           cellwidth = 36, cellheight = 9,
           angle_col = 45,
           main = ttl,
           filename = file.path(out_root, tier_name, "heatmap",
                                paste0("Heatmap_LFC_", file_suffix, ".png")),
           width = w, height = h)
  
  # Save matrix as CSV
  out_df <- data.frame(SYMBOL = rownames(lfc_mat),
                       SubPanel = as.character(sub_vec_f),
                       lfc_mat, check.names = FALSE)
  write.csv(out_df,
            file.path(out_root, tier_name, "heatmap",
                      paste0("LFC_matrix_", file_suffix, ".csv")),
            row.names = FALSE)
  
  invisible(list(lfc = lfc_mat, n_genes = nrow(lfc_mat)))
}

cat("\n===== Generating per-tier LFC heatmaps =====\n")
heatmap_log <- data.frame()
for (tn in tier_names) {
  # Super-category 1: Inflammation (7 sub-panels)
  r1 <- draw_lfc_heatmap(
    tier_name        = tn,
    super_cat        = "Inflammation",
    sub_panels_in_cat = c("Acute_Inflammation","TLR_NFkB_Signaling",
                          "Cytokine_Storm_Core","Chemokine_Recruitment",
                          "Interferon_Response",
                          "Macrophage_Induced_Inflammation",
                          "Lysosomal_Phagocytic"),
    deg_annotated_list = deg_annotated,
    contrast_tiers = contrast_tiers,
    file_suffix = "Inflammation"
  )
  
  # Super-category 2: Inflammatory_Apoptosis_Resolution (2 sub-panels)
  r2 <- draw_lfc_heatmap(
    tier_name        = tn,
    super_cat        = "Inflammatory_Apoptosis_Resolution",
    sub_panels_in_cat = c("Inflammatory_Cell_Death",
                          "Resolution_Anti_inflammatory"),
    deg_annotated_list = deg_annotated,
    contrast_tiers = contrast_tiers,
    file_suffix = "Apoptosis_Resolution"
  )
  
  n1 <- if (is.null(r1)) 0 else r1$n_genes
  n2 <- if (is.null(r2)) 0 else r2$n_genes
  heatmap_log <- rbind(heatmap_log,
                       data.frame(tier = tn,
                                  Inflammation_n_genes = n1,
                                  Apop_Resolution_n_genes = n2))
  cat("  ", tn, " : Inflammation=", n1,
      " genes, Apop/Resolution=", n2, " genes\n", sep = "")
}
write.csv(heatmap_log,
          file.path(out_root, "tables", "Heatmap_gene_counts_per_tier.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 8. Sub-panel module score (mean LFC per panel per contrast) ===
# =====================================================================
# Used for sub-panel-level summary heatmap (compact view)

build_module_score <- function(panel_name, deg_annotated_list) {
  panel_genes <- sub_panels[[panel_name]]
  scores <- sapply(names(deg_annotated_list), function(cn) {
    d <- deg_annotated_list[[cn]]
    lfc_vals <- d$log2FoldChange[d$SYMBOL %in% panel_genes & !is.na(d$log2FoldChange)]
    if (length(lfc_vals) == 0) return(NA)
    mean(lfc_vals, na.rm = TRUE)
  })
  scores
}

module_score_mat <- sapply(sub_panel_order, function(p) {
  build_module_score(p, deg_annotated)
})
module_score_mat <- t(module_score_mat)  # row = subpanel, col = contrast

# One-sample t-test against 0 (panel-level activation per contrast)
module_pval_mat <- sapply(names(deg_annotated), function(cn) {
  d <- deg_annotated[[cn]]
  sapply(sub_panel_order, function(p) {
    panel_genes <- sub_panels[[p]]
    lfc_vals <- d$log2FoldChange[d$SYMBOL %in% panel_genes & !is.na(d$log2FoldChange)]
    if (length(lfc_vals) < 3) return(NA)
    tryCatch(t.test(lfc_vals, mu = 0)$p.value, error = function(e) NA)
  })
})
rownames(module_pval_mat) <- sub_panel_order

# Save module scores
ms_long <- module_score_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SubPanel") %>%
  tidyr::pivot_longer(-SubPanel, names_to = "Contrast", values_to = "mean_LFC") %>%
  dplyr::left_join(
    module_pval_mat %>% as.data.frame() %>%
      tibble::rownames_to_column("SubPanel") %>%
      tidyr::pivot_longer(-SubPanel, names_to = "Contrast", values_to = "t_test_pval"),
    by = c("SubPanel","Contrast")
  ) %>%
  dplyr::mutate(SuperCat = super_cat_map[SubPanel])
write.csv(ms_long,
          file.path(out_root, "tables", "module_LFC_scores_long.csv"),
          row.names = FALSE)

# Per-tier sub-panel summary heatmap (compact)
draw_module_summary_heatmap <- function(tier_name) {
  cs <- contrast_tiers$contrast[contrast_tiers$tier == tier_name]
  if (length(cs) == 0) return(invisible(NULL))
  
  ms <- module_score_mat[, cs, drop = FALSE]
  pv <- module_pval_mat[, cs, drop = FALSE]
  
  sig_mat <- ifelse(is.na(pv), "",
                    ifelse(pv < 0.001, "***",
                           ifelse(pv < 0.01,  "**",
                                  ifelse(pv < 0.05,  "*", ""))))
  
  row_groups <- data.frame(
    SuperCat = factor(super_cat_map[rownames(ms)], levels = super_cat_order),
    row.names = rownames(ms)
  )
  ord <- order(row_groups$SuperCat)
  ms <- ms[ord, , drop = FALSE]
  pv <- pv[ord, , drop = FALSE]
  sig_mat <- sig_mat[ord, , drop = FALSE]
  row_groups <- row_groups[ord, , drop = FALSE]
  
  gaps_row <- cumsum(table(droplevels(row_groups$SuperCat)))
  gaps_row <- as.numeric(gaps_row[-length(gaps_row)])
  
  lfc_max <- max(abs(ms), na.rm = TRUE); lfc_max <- min(lfc_max, 3)
  breaks  <- seq(-lfc_max, lfc_max, length.out = 101)
  
  ctrl_group <- unique(contrast_tiers$control_group[contrast_tiers$tier == tier_name])
  ttl <- paste0("Sub-panel module score (mean LFC) — ", tier_name,
                "\n(Control = ", ctrl_group, ",  one-sample t-test vs 0)")
  
  pheatmap(ms,
           color = colorRampPalette(c("#0055FF","white","#C73E1D"))(100),
           breaks = breaks,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = row_groups,
           annotation_colors = list(SuperCat = super_cat_palette),
           gaps_row = gaps_row,
           display_numbers = sig_mat,
           number_color = "black", fontsize_number = 11,
           fontsize_row = 10, fontsize_col = 9,
           cellwidth = 50, cellheight = 22,
           angle_col = 45,
           main = ttl,
           filename = file.path(out_root, tier_name, "heatmap",
                                "Heatmap_module_LFC_summary.png"),
           width = max(7, ncol(ms) * 1.6 + 4), height = 6.5)
}

for (tn in tier_names) draw_module_summary_heatmap(tn)

# =====================================================================
# ===== 9. Volcano plot (sub-panel colored, per contrast) =============
# =====================================================================
make_inflammation_volcano <- function(df, contrast_name, tier_name,
                                      n_each = n_each_side) {
  
  x_min <- min(df$log2FoldChange, na.rm = TRUE)
  x_max <- max(df$log2FoldChange, na.rm = TRUE)
  y_max <- max(df$neg_log10_padj, na.rm = TRUE)
  
  x_left  <- x_min * 1.10; x_right <- x_max * 1.10
  y_top   <- y_max * 1.10
  x_bar_d <- x_min * 1.05; x_bar_u <- x_max * 1.05
  x_lab_d <- x_min * 0.93; x_lab_u <- x_max * 0.93
  y_lab_t <- y_max * 0.92; y_lab_b <- y_max * 0.10
  
  panel_sig <- df %>% dplyr::filter(SubPanel != "Other", sig)
  
  down_labels <- panel_sig %>%
    dplyr::filter(log2FoldChange < 0) %>%
    dplyr::arrange(SuperCat, SubPanel, dplyr::desc(neg_log10_padj)) %>%
    dplyr::slice_head(n = n_each)
  up_labels <- panel_sig %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::arrange(SuperCat, SubPanel, dplyr::desc(neg_log10_padj)) %>%
    dplyr::slice_head(n = n_each)
  
  spacing_max <- y_max * 0.045
  avail_h     <- y_lab_t - y_lab_b
  
  position_labels <- function(labs, x_pos) {
    if (nrow(labs) == 0) return(labs)
    n <- nrow(labs); needed <- (n - 1) * spacing_max
    labs$y_label <- if (n == 1) y_lab_t
    else if (needed < avail_h) seq(y_lab_t, y_lab_t - needed, length.out = n)
    else seq(y_lab_t, y_lab_b, length.out = n)
    labs$x_label <- x_pos
    labs
  }
  down_labels <- position_labels(down_labels, x_lab_d)
  up_labels   <- position_labels(up_labels,   x_lab_u)
  
  ti <- contrast_tiers[contrast_tiers$contrast == contrast_name, ]
  ttl <- paste0("[INFLAMMATION] Volcano — ", contrast_name,
                "\n", tier_name, "  |  Control = ", ti$control_group)
  
  n_sig_dn <- sum(df$sig & df$log2FoldChange < 0, na.rm = TRUE)
  n_sig_up <- sum(df$sig & df$log2FoldChange > 0, na.rm = TRUE)
  
  p <- ggplot() +
    # Background (Other / ns)
    geom_point(data = df %>% dplyr::filter(SubPanel == "Other", dot_class == "ns"),
               aes(x = log2FoldChange, y = neg_log10_padj),
               color = "grey85", alpha = 0.4, size = 0.7) +
    geom_point(data = df %>% dplyr::filter(SubPanel == "Other",
                                           dot_class %in% c("sub_up","sub_down")),
               aes(x = log2FoldChange, y = neg_log10_padj),
               color = "grey70", alpha = 0.5, size = 0.8) +
    geom_point(data = df %>% dplyr::filter(SubPanel == "Other",
                                           dot_class %in% c("sig_up","sig_down")),
               aes(x = log2FoldChange, y = neg_log10_padj),
               color = "grey50", alpha = 0.7, size = 1.0) +
    # Panel genes — color by sub-panel
    geom_point(data = df %>% dplyr::filter(SubPanel != "Other"),
               aes(x = log2FoldChange, y = neg_log10_padj, color = SubPanel),
               alpha = 0.75, size = 1.5) +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", color = "black", linewidth = 0.35) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed", color = "black", linewidth = 0.35)
  
  if (nrow(down_labels) > 0) {
    p <- p +
      geom_segment(data = down_labels,
                   aes(x = log2FoldChange, y = neg_log10_padj,
                       xend = x_label, yend = y_label,
                       color = SubPanel),
                   linetype = "dashed", linewidth = 0.2, alpha = 0.5,
                   show.legend = FALSE) +
      geom_point(data = down_labels,
                 aes(x = log2FoldChange, y = neg_log10_padj,
                     fill = SubPanel),
                 shape = 21, color = "black", stroke = 0.4, size = 2.4,
                 show.legend = FALSE) +
      geom_text(data = down_labels,
                aes(x = x_label, y = y_label, label = SYMBOL,
                    color = SubPanel),
                hjust = 0, vjust = 0.5,
                size = 2.9, fontface = "italic", family = plot_font,
                show.legend = FALSE)
  }
  if (nrow(up_labels) > 0) {
    p <- p +
      geom_segment(data = up_labels,
                   aes(x = log2FoldChange, y = neg_log10_padj,
                       xend = x_label, yend = y_label,
                       color = SubPanel),
                   linetype = "dashed", linewidth = 0.2, alpha = 0.5,
                   show.legend = FALSE) +
      geom_point(data = up_labels,
                 aes(x = log2FoldChange, y = neg_log10_padj,
                     fill = SubPanel),
                 shape = 21, color = "black", stroke = 0.4, size = 2.4,
                 show.legend = FALSE) +
      geom_text(data = up_labels,
                aes(x = x_label, y = y_label, label = SYMBOL,
                    color = SubPanel),
                hjust = 1, vjust = 0.5,
                size = 2.9, fontface = "italic", family = plot_font,
                show.legend = FALSE)
  }
  
  p <- p +
    annotate("text", x = x_left,  y = y_top, hjust = 0, vjust = 1,
             label = sprintf("Down\nSig.: %d", n_sig_dn),
             color = "#0055FF", size = 3.4, lineheight = 1.0,
             family = plot_font, fontface = "bold") +
    annotate("text", x = x_right, y = y_top, hjust = 1, vjust = 1,
             label = sprintf("Up\nSig.: %d", n_sig_up),
             color = "#C73E1D", size = 3.4, lineheight = 1.0,
             family = plot_font, fontface = "bold") +
    scale_color_manual(values = sub_panel_palette, name = "Sub-panel",
                       breaks = sub_panel_order, drop = TRUE) +
    scale_fill_manual(values = sub_panel_palette, name = "Sub-panel",
                      breaks = sub_panel_order, drop = TRUE) +
    coord_cartesian(xlim = c(x_left, x_right), ylim = c(0, y_top), clip = "off") +
    labs(x = expression(log[2]~"Fold Change"),
         y = expression(-log[10]~italic(P)["adj"]),
         title = ttl) +
    theme_classic(base_size = 11, base_family = plot_font) +
    theme(
      text             = element_text(family = plot_font),
      plot.title       = element_text(face = "bold", size = 12,
                                      hjust = 0.5, lineheight = 1.1),
      axis.title       = element_text(face = "bold", size = 11),
      axis.text        = element_text(size = 10, color = "black"),
      legend.position  = "right",
      legend.text      = element_text(size = 8),
      legend.title     = element_text(size = 9, face = "bold"),
      legend.key.size  = unit(0.4, "cm"),
      panel.grid       = element_blank(),
      plot.margin      = margin(15, 15, 15, 15)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
           fill  = "none")
  
  return(p)
}

cat("\n===== Generating volcanoes (sub-panel colored) =====\n")
for (nm in names(deg_annotated)) {
  tn <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  p  <- make_inflammation_volcano(deg_annotated[[nm]], nm, tn)
  ggsave(file.path(out_root, tn, "volcano",
                   paste0("Volcano_INFLAM_", nm, ".png")),
         p, width = 12, height = 9, dpi = 300, bg = "white")
  ggsave(file.path(out_root, tn, "volcano",
                   paste0("Volcano_INFLAM_", nm, ".pdf")),
         p, width = 12, height = 9, device = cairo_pdf)
  cat("  ", tn, " — ", nm, "\n", sep = "")
}

# =====================================================================
# ===== 10. KEGG enrichment (per contrast, inflammation-focused) ======
# =====================================================================
universe_entrez <- unique(na.omit(as.character(anno$ENTREZID)))

run_kegg_ora <- function(deg_df, padj_cut = padj_cutoff,
                         lfc_cut = lfc_cutoff, universe) {
  res_list <- list()
  gene_lists <- list(
    up   = deg_df %>% dplyr::filter(padj < padj_cut, log2FoldChange >  lfc_cut) %>% dplyr::pull(ENTREZID),
    down = deg_df %>% dplyr::filter(padj < padj_cut, log2FoldChange < -lfc_cut) %>% dplyr::pull(ENTREZID),
    all  = deg_df %>% dplyr::filter(padj < padj_cut, abs(log2FoldChange) > lfc_cut) %>% dplyr::pull(ENTREZID)
  )
  for (set_nm in names(gene_lists)) {
    g <- as.character(gene_lists[[set_nm]])
    if (length(g) < 10) { res_list[[set_nm]] <- NULL; next }
    res <- tryCatch(
      enrichKEGG(gene = g, organism = "mmu", universe = universe,
                 pvalueCutoff = 0.1, qvalueCutoff = 0.25, keyType = "kegg"),
      error = function(e) { message("ORA failed: ", e$message); NULL }
    )
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      res <- setReadable(res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    }
    res_list[[set_nm]] <- res
  }
  res_list
}

run_kegg_gsea <- function(deg_df) {
  rnk <- deg_df %>%
    dplyr::filter(!is.na(ENTREZID), !is.na(log2FoldChange)) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))
  geneList <- setNames(rnk$log2FoldChange, as.character(rnk$ENTREZID))
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  res <- tryCatch(
    gseKEGG(geneList = geneList, organism = "mmu",
            pvalueCutoff = 0.25, minGSSize = 10, maxGSSize = 500,
            seed = TRUE, verbose = FALSE),
    error = function(e) { message("GSEA failed: ", e$message); NULL }
  )
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    res <- setReadable(res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  }
  res
}

cat("\n===== Running KEGG ORA + GSEA per contrast =====\n")
kegg_ora  <- lapply(deg_gene, run_kegg_ora, universe = universe_entrez)
kegg_gsea <- lapply(deg_gene, run_kegg_gsea)
saveRDS(kegg_ora,  file.path(out_root, "kegg_ora_results.rds"))
saveRDS(kegg_gsea, file.path(out_root, "kegg_gsea_results.rds"))

# Export per-tier CSVs
for (nm in names(kegg_ora)) {
  tn <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  for (set_nm in names(kegg_ora[[nm]])) {
    r <- kegg_ora[[nm]][[set_nm]]
    if (!is.null(r) && nrow(as.data.frame(r)) > 0) {
      write.csv(as.data.frame(r),
                file.path(out_root, tn, "enrich",
                          paste0("ORA_", nm, "_", set_nm, ".csv")),
                row.names = FALSE)
    }
  }
  r <- kegg_gsea[[nm]]
  if (!is.null(r) && nrow(as.data.frame(r)) > 0) {
    write.csv(as.data.frame(r),
              file.path(out_root, tn, "enrich",
                        paste0("GSEA_", nm, ".csv")),
              row.names = FALSE)
  }
}

# Inflammation-focused KEGG categories
pathway_categories <- list(
  Inflammation = c(
    "mmu04060","mmu04061","mmu04062","mmu04064","mmu04066","mmu04612",
    "mmu04620","mmu04621","mmu04622","mmu04623","mmu04625","mmu04630",
    "mmu04640","mmu04650","mmu04657","mmu04658","mmu04660","mmu04666",
    "mmu04668","mmu04672","mmu05321"
  ),
  Apop_Resolution = c(
    "mmu04210","mmu04215","mmu04216","mmu04217","mmu04115",
    "mmu04141","mmu04142","mmu04145","mmu04140","mmu04137","mmu04218"
  )
)

make_inflam_lollipop <- function(ora_res, pathway_cats, contrast_name,
                                 tier_name, top_per_cat = 12) {
  if (is.null(ora_res) || nrow(as.data.frame(ora_res)) == 0) return(NULL)
  df <- as.data.frame(ora_res)
  df$Category <- NA
  for (cat_nm in names(pathway_cats)) {
    df$Category[df$ID %in% pathway_cats[[cat_nm]]] <- cat_nm
  }
  df <- df[!is.na(df$Category), ]
  if (nrow(df) == 0) return(NULL)
  df$logp <- -log10(df$pvalue)
  df <- df %>%
    dplyr::group_by(Category) %>%
    dplyr::slice_max(logp, n = top_per_cat) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Category, logp)
  df$Description <- factor(df$Description, levels = df$Description)
  df$Category    <- factor(df$Category,
                           levels = c("Inflammation","Apop_Resolution"))
  cat_colors <- c(Inflammation = "#C73E1D", Apop_Resolution = "#1B7837")
  
  ti <- contrast_tiers[contrast_tiers$contrast == contrast_name, ]
  
  ggplot(df, aes(x = logp, y = Description, color = Category)) +
    geom_segment(aes(x = 0, xend = logp, yend = Description), linewidth = 0.6) +
    geom_point(aes(size = Count)) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(values = cat_colors, guide = "none") +
    scale_size_continuous(range = c(3, 8), name = "Gene count") +
    labs(title = paste0("[INFLAMMATION] KEGG enrichment — ", contrast_name),
         subtitle = paste0(tier_name, "  |  Control = ", ti$control_group),
         x = expression(-log[10]~"(p-value)"), y = NULL) +
    theme_bw(base_size = 11, base_family = plot_font) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "grey92"),
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"))
}

cat("\n===== Generating KEGG lollipop plots =====\n")
for (nm in names(kegg_ora)) {
  if (is.null(kegg_ora[[nm]]$all)) next
  tn <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  p <- make_inflam_lollipop(kegg_ora[[nm]]$all, pathway_categories, nm, tn)
  if (!is.null(p)) {
    ggsave(file.path(out_root, tn, "enrich",
                     paste0("Lollipop_INFLAM_", nm, ".png")),
           p, width = 10, height = 9, dpi = 300)
  }
}

# Disease-pathway blacklist for dotplot
dotplot_blacklist <- c(
  "Herpes simplex virus 1 infection","Epstein-Barr virus infection",
  "Kaposi sarcoma-associated herpesvirus infection",
  "Human papillomavirus infection","Human cytomegalovirus infection",
  "Human immunodeficiency virus 1 infection",
  "Human T-cell leukemia virus 1 infection",
  "Hepatitis B","Hepatitis C","Measles","Influenza A",
  "Coronavirus disease - COVID-19","Viral myocarditis",
  "Viral carcinogenesis","Viral life cycle - HIV-1",
  "Staphylococcus aureus infection","Tuberculosis","Leishmaniasis",
  "Malaria","Chagas disease","African trypanosomiasis","Amoebiasis",
  "Toxoplasmosis","Pertussis","Legionellosis","Salmonella infection",
  "Yersinia infection","Vibrio cholerae infection",
  "Pathogenic Escherichia coli infection","Shigellosis",
  "Bacterial invasion of epithelial cells",
  "Epithelial cell signaling in Helicobacter pylori infection",
  "Pathways in cancer","Transcriptional misregulation in cancer",
  "MicroRNAs in cancer","Proteoglycans in cancer",
  "Chemical carcinogenesis - DNA adducts",
  "Chemical carcinogenesis - receptor activation",
  "Chemical carcinogenesis - reactive oxygen species",
  "Dilated cardiomyopathy","Hypertrophic cardiomyopathy",
  "Arrhythmogenic right ventricular cardiomyopathy",
  "Alzheimer disease","Parkinson disease","Huntington disease",
  "Amyotrophic lateral sclerosis","Prion disease",
  "Pathways of neurodegeneration - multiple diseases"
)
filter_kegg_result <- function(res, blacklist) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(res)
  res_new <- res
  res_new@result <- res@result[!res@result$Description %in% blacklist, ]
  res_new
}

cat("\n===== Generating filtered KEGG dotplots =====\n")
for (nm in names(kegg_ora)) {
  res_full <- kegg_ora[[nm]]$all
  if (is.null(res_full) || nrow(as.data.frame(res_full)) == 0) next
  res_filt <- filter_kegg_result(res_full, dotplot_blacklist)
  n_filt   <- nrow(as.data.frame(res_filt))
  if (n_filt == 0) next
  tn <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  
  p_filt <- dotplot(res_filt, showCategory = min(20, n_filt), font.size = 9) +
    ggtitle(paste0("[INFLAMMATION] KEGG ORA (filtered) — ", nm, "\n", tn)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          text = element_text(family = plot_font))
  ggsave(file.path(out_root, tn, "enrich",
                   paste0("Dotplot_KEGG_filtered_", nm, ".png")),
         p_filt, width = 9, height = 9, dpi = 300)
}

# GSEA enrichment plots for key inflammation pathways
target_pathways <- c(
  "mmu04060","mmu04062","mmu04064","mmu04066","mmu04612","mmu04620",
  "mmu04621","mmu04622","mmu04625","mmu04630","mmu04657","mmu04666",
  "mmu04668","mmu05321",
  "mmu04210","mmu04216","mmu04217","mmu04115","mmu04141","mmu04142","mmu04140"
)

cat("\n===== Generating GSEA plots for target pathways =====\n")
for (nm in names(kegg_gsea)) {
  gsea_res <- kegg_gsea[[nm]]
  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) next
  hits <- intersect(target_pathways, gsea_res@result$ID)
  if (length(hits) == 0) next
  tn <- contrast_tiers$tier[contrast_tiers$contrast == nm]
  for (pid in hits) {
    desc <- gsea_res@result$Description[gsea_res@result$ID == pid]
    p <- tryCatch(
      gseaplot2(gsea_res, geneSetID = pid,
                title = paste0("[INFLAMMATION] ", nm, "  |  ", tn,
                               "\n", desc),
                pvalue_table = TRUE),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      ggsave(file.path(out_root, tn, "enrich",
                       paste0("GSEA_", nm, "_", pid, ".png")),
             p, width = 8, height = 6.5, dpi = 300)
    }
  }
}

# =====================================================================
# ===== 11. GSVA supplementary (sample-level, expression-based) =======
# =====================================================================
cat("\n===== Computing GSVA scores (supplementary) =====\n")

# Build sample-wise expression matrix at gene level
# (collapse transcripts → max-expression isoform per gene)
vst_mat <- assay(vsd)
tx_ids_nov <- sub("\\..*$", "", rownames(vst_mat))
sym_map    <- anno$SYMBOL[match(tx_ids_nov, anno$TXID)]

# Collapse to gene level (best isoform = highest mean expression)
gene_expr_list <- split(seq_len(nrow(vst_mat)), sym_map)
gene_expr_list <- gene_expr_list[!is.na(names(gene_expr_list)) &
                                   names(gene_expr_list) != ""]
gene_expr_mat <- t(sapply(gene_expr_list, function(idx) {
  if (length(idx) == 1) return(vst_mat[idx, ])
  means <- rowMeans(vst_mat[idx, , drop = FALSE])
  vst_mat[idx[which.max(means)], ]
}))
rownames(gene_expr_mat) <- names(gene_expr_list)

# Subset to genes present in any sub-panel
panel_genes_all <- unique(unlist(sub_panels))
gene_expr_panel <- gene_expr_mat[rownames(gene_expr_mat) %in% panel_genes_all, ]
cat("  Gene-level expression matrix: ", nrow(gene_expr_panel), " genes × ",
    ncol(gene_expr_panel), " samples\n", sep = "")

# GSVA — sub-panel level (modern API: gsvaParam)
gsva_res <- tryCatch({
  gsva_param <- gsvaParam(exprData = gene_expr_panel,
                          geneSets = sub_panels,
                          minSize = 5, maxSize = 500,
                          kcdf = "Gaussian")
  gsva(gsva_param, verbose = FALSE)
}, error = function(e) {
  # Fallback to legacy API
  message("Using legacy GSVA API: ", e$message)
  gsva(gene_expr_panel, sub_panels, method = "gsva",
       kcdf = "Gaussian", verbose = FALSE)
})

# GSVA heatmap (sample-level)
group_factor <- factor(as.character(colData(dds)$group),
                       levels = names(group_colors_full))
col_ord <- order(group_factor)
gsva_sorted <- gsva_res[, col_ord]
ann_col_gsva <- data.frame(group = group_factor[col_ord],
                           row.names = colnames(gsva_sorted))
gaps_col <- cumsum(table(group_factor))
gaps_col <- as.numeric(gaps_col[-length(gaps_col)])

row_groups_gsva <- data.frame(
  SuperCat = factor(super_cat_map[rownames(gsva_sorted)],
                    levels = super_cat_order),
  row.names = rownames(gsva_sorted)
)
ord <- order(row_groups_gsva$SuperCat)
gsva_sorted <- gsva_sorted[ord, , drop = FALSE]
row_groups_gsva <- row_groups_gsva[ord, , drop = FALSE]
gaps_row_gsva <- cumsum(table(droplevels(row_groups_gsva$SuperCat)))
gaps_row_gsva <- as.numeric(gaps_row_gsva[-length(gaps_row_gsva)])

pheatmap(gsva_sorted,
         color = colorRampPalette(c("#0055FF","white","#C73E1D"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_groups_gsva,
         annotation_col = ann_col_gsva,
         annotation_colors = list(SuperCat = super_cat_palette,
                                  group = group_colors_full),
         gaps_row = gaps_row_gsva, gaps_col = gaps_col,
         fontsize_row = 10, fontsize_col = 7,
         cellwidth = 14, cellheight = 22, angle_col = 45,
         main = "[Supplementary] GSVA sub-panel score — All 24 samples (8 groups)",
         filename = file.path(out_root, "GSVA_supplementary", "heatmap",
                              "GSVA_subpanel_24samples.png"),
         width = 12, height = 6.5)

# Group-mean GSVA
group_vec_chr <- as.character(colData(dds)$group)
gsva_group_means <- sapply(names(group_colors_full), function(g) {
  rowMeans(gsva_res[, group_vec_chr == g, drop = FALSE], na.rm = TRUE)
})

row_groups_gsva_gm <- data.frame(
  SuperCat = factor(super_cat_map[rownames(gsva_group_means)],
                    levels = super_cat_order),
  row.names = rownames(gsva_group_means)
)
ord <- order(row_groups_gsva_gm$SuperCat)
gsva_group_means <- gsva_group_means[ord, , drop = FALSE]
row_groups_gsva_gm <- row_groups_gsva_gm[ord, , drop = FALSE]
gaps_row_gm <- cumsum(table(droplevels(row_groups_gsva_gm$SuperCat)))
gaps_row_gm <- as.numeric(gaps_row_gm[-length(gaps_row_gm)])

pheatmap(gsva_group_means,
         color = colorRampPalette(c("#0055FF","white","#C73E1D"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_groups_gsva_gm,
         annotation_colors = list(SuperCat = super_cat_palette),
         gaps_row = gaps_row_gm,
         fontsize_row = 10, fontsize_col = 9,
         cellwidth = 36, cellheight = 22, angle_col = 45,
         main = "[Supplementary] GSVA sub-panel — Group mean (8 groups)",
         filename = file.path(out_root, "GSVA_supplementary", "heatmap",
                              "GSVA_subpanel_groupmean.png"),
         width = 9, height = 6.5)

# Save GSVA matrix
gsva_long <- gsva_res %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SubPanel") %>%
  tidyr::pivot_longer(-SubPanel, names_to = "sample",
                      values_to = "GSVA_score") %>%
  dplyr::mutate(group = as.character(colData(dds)$group)[match(sample, rownames(colData(dds)))],
                SuperCat = super_cat_map[SubPanel])
write.csv(gsva_long,
          file.path(out_root, "GSVA_supplementary", "GSVA_scores_long.csv"),
          row.names = FALSE)

# =====================================================================
# ===== 12. Final summary =============================================
# =====================================================================
cat("\n===== Step 3 complete =====\n")
cat("Output root: ", out_root, "/\n\n", sep = "")
cat("Per-tier outputs (Tier1 ~ Tier6):\n")
cat("  heatmap/  : Heatmap_LFC_Inflammation.png            (7 sub-panels, LFC)\n")
cat("              Heatmap_LFC_Apoptosis_Resolution.png    (2 sub-panels, LFC)\n")
cat("              Heatmap_module_LFC_summary.png          (9-panel compact)\n")
cat("              LFC_matrix_*.csv\n")
cat("  volcano/  : Volcano_INFLAM_<contrast>.{png,pdf}     (9 sub-panel colored)\n")
cat("  enrich/   : Lollipop_INFLAM_<contrast>.png          (Inflam vs Apop/Resol)\n")
cat("              Dotplot_KEGG_filtered_<contrast>.png    (disease pathways removed)\n")
cat("              GSEA_<contrast>_<pathway>.png\n")
cat("              ORA_<contrast>_{up,down,all}.csv\n")
cat("              GSEA_<contrast>.csv\n\n")
cat("Tables/:\n")
cat("  panel_definitions.csv                  (9 sub-panel full gene list)\n")
cat("  module_LFC_scores_long.csv             (contrast × subpanel × LFC mean)\n")
cat("  Heatmap_gene_counts_per_tier.csv\n\n")
cat("GSVA_supplementary/:\n")
cat("  heatmap/GSVA_subpanel_24samples.png    (sample-level)\n")
cat("  heatmap/GSVA_subpanel_groupmean.png    (8-group mean)\n")
cat("  GSVA_scores_long.csv\n")