# =====================================================================
# Step 8b: Final stem cell GO BP bubble plots (curated pathway list)
# Generates 2 figures: Layout A (basal, 4 contrasts) + Layout B (full, 14)
# =====================================================================

# ===== Curated pathway list with category assignment =================
curated_pathways <- tibble::tribble(
  ~Category,                       ~Description,
  # Category 1: Stem cell core
  "Stem cell core",                "stem cell development",
  "Stem cell core",                "stem cell proliferation",
  "Stem cell core",                "regulation of stem cell proliferation",
  "Stem cell core",                "regulation of hematopoietic stem cell proliferation",
  "Stem cell core",                "cell fate commitment",
  "Stem cell core",                "negative regulation of cell fate commitment",
  # Category 2: Intestinal epithelium
  "Intestinal epithelium",         "epithelial cell proliferation",
  "Intestinal epithelium",         "regulation of epithelial cell proliferation",
  "Intestinal epithelium",         "positive regulation of epithelial cell proliferation",
  "Intestinal epithelium",         "regulation of epithelial cell differentiation",
  "Intestinal epithelium",         "intestinal epithelial cell differentiation",
  "Intestinal epithelium",         "maintenance of gastrointestinal epithelium",
  "Intestinal epithelium",         "intestinal absorption",
  "Intestinal epithelium",         "regulation of intestinal absorption",
  "Intestinal epithelium",         "regulation of intestinal lipid absorption",
  # Category 3: Signaling
  #"Signaling",                     #"positive regulation of Notch signaling pathway",
  #"Signaling",                     #"non-canonical Wnt signaling pathway",
  #"Signaling",                     #"negative regulation of canonical Wnt signaling pathway",
  # Category 4: Epithelial architecture / Polarity
  "Epithelial architecture",       "cell-cell junction organization",
  "Epithelial architecture",       "cell junction assembly",
  "Epithelial architecture",       "regulation of cell junction assembly",
  "Epithelial architecture",       "cell-cell junction maintenance",
  "Epithelial architecture",       "cell junction disassembly",
  "Epithelial architecture",       "tight junction organization",
  "Epithelial architecture",       "establishment of apical/basal cell polarity"
)

category_order <- c("Stem cell core","Intestinal epithelium",
                    "Signaling","Epithelial architecture")
curated_pathways$Category <- factor(curated_pathways$Category,
                                    levels = category_order)

# Category color palette (consistent with manuscript narrative)
category_palette <- c(
  "Stem cell core"          = "#2E8B57",   # sea green (ISC color)
  "Intestinal epithelium"   = "#8C510A",   # brown (Enterocyte color)
  "Signaling"               = "#4A90E2",   # blue (Wnt/Notch color)
  "Epithelial architecture" = "#7BAFD4"    # light blue
)

cat("Curated pathways: ", nrow(curated_pathways), " total across ",
    length(category_order), " categories\n", sep = "")

build_curated_long <- function(enrich_list, baseline_map) {
  c2b <- list()
  for (bl in names(baseline_map)) {
    for (cn in baseline_map[[bl]]) c2b[[cn]] <- bl
  }
  
  # ===== Build matched data =====
  rows <- list()
  for (cn in names(c2b)) {
    r <- enrich_list[[cn]]
    if (is.null(r) || nrow(as.data.frame(r)) == 0) next
    df <- as.data.frame(r)
    sub <- df[df$Description %in% curated_pathways$Description, , drop = FALSE]
    if (nrow(sub) == 0) next
    sub$Contrast <- cn
    sub$Baseline <- c2b[[cn]]
    rows[[cn]] <- sub[, c("ID","Description","Count","pvalue","p.adjust",
                          "Contrast","Baseline")]
  }
  out <- if (length(rows) > 0) do.call(rbind, rows) else NULL
  
  # ===== ★ Add placeholder rows for ALL contrasts × ALL pathways =====
  # This ensures empty contrasts (no enrichment) are still shown as
  # empty columns, and empty pathways still appear as blank rows.
  placeholder <- expand.grid(
    Contrast    = names(c2b),
    Description = curated_pathways$Description,
    stringsAsFactors = FALSE
  )
  placeholder$Baseline <- unlist(c2b[placeholder$Contrast])
  placeholder$ID       <- NA_character_
  placeholder$Count    <- NA_real_
  placeholder$pvalue   <- NA_real_
  placeholder$p.adjust <- NA_real_
  
  # Order columns same as `out`
  placeholder <- placeholder[, c("ID","Description","Count","pvalue",
                                 "p.adjust","Contrast","Baseline")]
  
  # Combine: real data first, then placeholders for missing combos
  if (!is.null(out)) {
    # Mark which combos already have real data
    out$key <- paste(out$Contrast, out$Description, sep = "||")
    placeholder$key <- paste(placeholder$Contrast,
                             placeholder$Description, sep = "||")
    placeholder <- placeholder[!placeholder$key %in% out$key, , drop = FALSE]
    out$key <- NULL
    placeholder$key <- NULL
    out <- rbind(out, placeholder)
  } else {
    out <- placeholder
  }
  
  out$logp <- -log10(out$pvalue)    # NA stays NA
  out <- dplyr::left_join(out, curated_pathways, by = "Description")
  out
}

# ===== Build curated bubble plot =====================================
make_curated_bubble <- function(enrich_list,
                                baseline_map,
                                baseline_levels,
                                baseline_palette,
                                title_main,
                                title_sub) {
  long <- build_curated_long(enrich_list, baseline_map)
  if (is.null(long) || nrow(long) == 0) {
    message("No curated pathways found in enrichment")
    return(NULL)
  }
  long$Baseline <- factor(long$Baseline, levels = baseline_levels)
  
  # Pathway ordering: by category (preserve curated_pathways order)
  pathway_order_factor <- factor(curated_pathways$Description,
                                 levels = curated_pathways$Description)
  long$Description <- factor(long$Description,
                             levels = rev(curated_pathways$Description))
  
  # Contrast order
  contrast_order <- unlist(baseline_map, use.names = FALSE)
  long$Contrast <- factor(long$Contrast, levels = contrast_order)
  
  # Build plot
  p <- ggplot(long, aes(x = Contrast, y = Description)) +
    geom_point(aes(size = Count, fill = logp),
               shape = 21, color = "black", stroke = 0.4, alpha = 0.95,
               na.rm = TRUE)            # ★ NA 값은 그냥 안 그림 (placeholder)
  
  if (has_ggh4x) {
    # Both row (Category) and column (Baseline) colored strips
    strip_fills_x <- unname(baseline_palette[baseline_levels])
    # For Y: use reversed category order to match reversed pathways
    strip_fills_y <- unname(category_palette[rev(category_order)])
    
    p <- p + ggh4x::facet_grid2(
      Category ~ Baseline,
      scales = "free", space = "free",
      labeller = labeller(Baseline = baseline_strip_labels),
      strip = ggh4x::strip_themed(
        # Column (Baseline) strips
        background_x = lapply(strip_fills_x, function(col)
          element_rect(fill = col, color = "black", linewidth = 0.8)),
        text_x = lapply(strip_fills_x, function(col) {
          text_col <- ifelse(col %in% c("#1F4E79","#A61C00"), "white", "black")
          element_text(color = text_col, face = "bold", size = 11,
                       family = plot_font, lineheight = 0.95)
        }),
        # Row (Category) strips
        background_y = lapply(unname(category_palette[category_order]),
                              function(col) {
                                element_rect(fill = col, color = "black",
                                             linewidth = 0.8)
                              }),
        text_y = lapply(unname(category_palette[category_order]),
                        function(col) {
                          text_col <- ifelse(col %in% c("#2E8B57","#8C510A","#4A90E2"),
                                             "white", "black")
                          element_text(color = text_col, face = "bold",
                                       size = 10, family = plot_font,
                                       angle = 90, lineheight = 0.95)
                        })
      )
    )
  } else {
    p <- p + facet_grid(Category ~ Baseline,
                        scales = "free", space = "free",
                        labeller = labeller(Baseline = baseline_strip_labels))
  }
  
  p <- p +
    scale_fill_gradient2(
      low = "#5865f2", mid = "#faffff", high = "#ff3333",
      midpoint = 3,                          # ★ GO BP는 -log10p가 KEGG보다 낮으므로 midpoint=3
      name = expression(-log[10]~italic(P)),
      limits = c(0, NA),
      guide = guide_colorbar(barwidth = 1, barheight = 8)
    ) +
    scale_size_continuous(range = c(2, 8), name = "Gene count",
                          breaks = c(10, 30, 60)) +
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
      panel.spacing.y     = unit(0.4, "lines"),
      strip.placement     = "outside",
      legend.position     = "right",
      legend.title        = element_text(face = "bold", size = 10),
      legend.text         = element_text(size = 9),
      legend.box          = "vertical",
      legend.spacing.y    = unit(0.3, "cm"),
      plot.margin         = margin(15, 15, 15, 15)
    )
  
  list(plot = p, data = long, n_pathways = length(unique(long$Description)))
}

# ===== Generate plots ================================================
generate_curated_bubble <- function(enrich_list, layout_label) {
  if (layout_label == "basal") {
    bmap <- baseline_map_basal
    bl   <- baseline_levels_basal
    pal  <- baseline_palette_basal
    sub_label <- "basal (no TNFa) - 4 contrasts"
    fig_prefix <- "Fig_5"
    fig_id <- "f_GOBP"
    fig_w <- 11
  } else {
    bmap <- baseline_map_full
    bl   <- baseline_levels_full
    pal  <- baseline_palette_full
    sub_label <- "full 14 contrasts"
    fig_prefix <- "Fig_S5"
    fig_id <- "f_GOBP"
    fig_w <- 15
  }
  
  title_main <- paste0("GO BP - Stem cell + lineage curated pathways (", sub_label, ")")
  title_sub  <- paste0(nrow(curated_pathways),
                       " curated pathways across 4 categories; bubble size = gene count, color = -log10(p)")
  
  result <- make_curated_bubble(enrich_list, bmap, bl, pal,
                                title_main, title_sub)
  
  if (is.null(result)) {
    cat("  [", layout_label, "] No curated pathways found in enrichment\n", sep = "")
    return(invisible(NULL))
  }
  
  fname <- paste0(fig_prefix, fig_id, "_curated_", layout_label)
  fig_h <- max(10, result$n_pathways * 0.35 + 5)
  
  ggsave(file.path(out_root, paste0(fname, ".png")),
         result$plot, width = fig_w, height = fig_h,
         dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(out_root, paste0(fname, ".pdf")),
         result$plot, width = fig_w, height = fig_h,
         device = cairo_pdf, limitsize = FALSE)
  
  write.csv(result$data,
            file.path(out_root, "tables",
                      paste0("GOBP_curated_bubble_", layout_label, ".csv")),
            row.names = FALSE)
  
  cat("  [", layout_label, "] ", result$n_pathways,
      " pathways shown -> ", fname, ".png\n", sep = "")
}

cat("\n===== Generating curated GO BP bubble plots =====\n")
cat("\nLayout A (basal, no TNFa):\n")
generate_curated_bubble(go_bp_list, "basal")

cat("\nLayout B (full 14 contrasts):\n")
generate_curated_bubble(go_bp_list, "full")

cat("\n===== Curated bubble plots complete =====\n")


# =====================================================================
# Bubble plot — print(p) 가능 + 폰트 중앙 관리
# =====================================================================

# =====================================================================
# Bubble plot — print(p) 가능 + 폰트 중앙 관리 + Fig 4f 스타일 strip
# =====================================================================

# ===== 폰트 / 크기 — 튜닝 포인트 (여기만 만져서 조정) ================
bubble_params <- list(
  # 텍스트
  title        = 22,    # ★ 메인 타이틀
  subtitle     = 16,    # ★ 서브타이틀
  axis_x       = 16,    # ★ 하단 contrast 라벨
  axis_y       = 18,    # ★ 좌측 pathway 이름
  strip_x      = 13,    # ★ 상단 베이스라인 색 막대 두께 (글자는 투명)
  strip_y      = 14,    # ★ 우측 카테고리 라벨 (Stem cell core 등)
  legend_ttl   = 15,    # ★ 범례 제목
  legend_txt   = 14,    # ★ 범례 값
  # 버블
  size_range   = c(3, 12),   # ★ 버블 최소/최대 크기
  stroke       = 0.5,        # ★ 버블 외곽선
  # 색
  midpoint     = 3,          # ★ -log10p 색 변환점 (GO BP=3, KEGG=5 권장)
  # 디바이스
  panel_spacing = 0.5        # ★ 패널 간 간격 (lines)
)

# ===== Contrast label map — 모두 "Org+BMDM ..." 통일 =================
# (basal: 4 contrasts, full: 14 contrasts 모두 포함)
contrast_label_map <- c(
  # Basal layout (4 contrasts)
  T1_OrgBMDM1K_vs_Org         = "Org+BMDM 1K",
  T1_OrgBMDM5K_vs_Org         = "Org+BMDM 5K",
  T2_OrgBMDM1K_vs_Mac         = "Org+BMDM 1K",
  T2_OrgBMDM5K_vs_Mac         = "Org+BMDM 5K",
  # Full layout 추가 contrasts
  T3_OrgTNFa_vs_Org           = "Org+TNFα",
  T3_OrgTNFaBMDM1K_vs_Org     = "Org+BMDM 1K+TNFα",
  T3_OrgTNFaBMDM5K_vs_Org     = "Org+BMDM 5K+TNFα",
  T4_OrgTNFaBMDM1K_vs_OrgTNFa = "Org+BMDM 1K+TNFα",
  T4_OrgTNFaBMDM5K_vs_OrgTNFa = "Org+BMDM 5K+TNFα",
  T5_MacTNFa_vs_Mac           = "Org+TNFα",
  T5_OrgTNFaBMDM1K_vs_Mac     = "Org+BMDM 1K+TNFα",
  T5_OrgTNFaBMDM5K_vs_Mac     = "Org+BMDM 5K+TNFα",
  T6_OrgTNFaBMDM1K_vs_MacTNFa = "Org+BMDM 1K+TNFα",
  T6_OrgTNFaBMDM5K_vs_MacTNFa = "Org+BMDM 5K+TNFα"
)

# ===== make_curated_bubble — strip 색 막대만 + 우측 Baseline 범례 ====
make_curated_bubble <- function(enrich_list,
                                baseline_map,
                                baseline_levels,
                                baseline_palette,
                                title_main,
                                title_sub,
                                fs = bubble_params) {
  long <- build_curated_long(enrich_list, baseline_map)
  if (is.null(long) || nrow(long) == 0) {
    message("No curated pathways found in enrichment")
    return(NULL)
  }
  long$Baseline <- factor(long$Baseline, levels = baseline_levels)
  long$Description <- factor(long$Description,
                             levels = rev(curated_pathways$Description))
  
  contrast_order <- unlist(baseline_map, use.names = FALSE)
  long$Contrast <- factor(long$Contrast, levels = contrast_order)
  
  # ★ Baseline 색 범례용 dummy 데이터
  baseline_legend_data <- data.frame(
    Baseline = factor(baseline_levels, levels = baseline_levels)
  )
  
  p <- ggplot(long, aes(x = Contrast, y = Description)) +
    geom_point(aes(size = Count, fill = logp),
               shape = 21, color = "black",
               stroke = fs$stroke, alpha = 0.95,
               na.rm = TRUE) +
    # ★ Baseline 색 범례 생성용 invisible dummy
    geom_point(data = baseline_legend_data,
               aes(color = Baseline),
               x = NA, y = NA, size = 0, na.rm = TRUE,
               inherit.aes = FALSE, show.legend = TRUE)
  
  if (has_ggh4x) {
    strip_fills_x <- unname(baseline_palette[baseline_levels])
    
    p <- p + ggh4x::facet_grid2(
      Category ~ Baseline,
      scales = "free", space = "free",
      strip = ggh4x::strip_themed(
        # ★ Column (Baseline) strips — 색 막대만, 글자 투명
        background_x = lapply(strip_fills_x, function(col)
          element_rect(fill = col, color = "black", linewidth = 0.6)),
        text_x = lapply(strip_fills_x, function(col)
          element_text(color  = NA,                  # ★ 글자 투명
                       size   = fs$strip_x,
                       margin = margin(t = 5, b = 5))),
        # Row (Category) strips — 글자 보임 (그대로 유지)
        background_y = lapply(unname(category_palette[category_order]),
                              function(col)
                                element_rect(fill = col, color = "black",
                                             linewidth = 0.8)),
        text_y = lapply(unname(category_palette[category_order]),
                        function(col) {
                          text_col <- ifelse(col %in% c("#2E8B57","#8C510A","#4A90E2"),
                                             "white", "black")
                          element_text(color = text_col, face = "bold",
                                       size = fs$strip_y,
                                       family = plot_font,
                                       angle = 90, lineheight = 0.95)
                        })
      )
    )
  } else {
    p <- p + facet_grid(Category ~ Baseline,
                        scales = "free", space = "free")
  }
  
  p <- p +
    scale_fill_gradient2(
      low = "#5865f2", mid = "#faffff", high = "#ff3333",
      midpoint = fs$midpoint,
      name   = expression(-log[10]~italic(P)),
      limits = c(0, NA),
      guide  = guide_colorbar(barwidth = 1.4, barheight = 12, order = 1)
    ) +
    scale_size_continuous(
      range  = fs$size_range,
      name   = "Gene count",
      breaks = c(10, 30, 60),
      guide  = guide_legend(order = 2)
    ) +
    # ★ 오른쪽 Baseline 색 범례
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
         title = title_main, subtitle = title_sub) +
    theme_bw(base_size = 11, base_family = plot_font) +
    theme(
      text             = element_text(family = plot_font),
      plot.title       = element_text(face = "bold", hjust = 0.5,
                                      size = fs$title,
                                      margin = margin(b = 6)),
      plot.subtitle    = element_text(hjust = 0.5,
                                      size = fs$subtitle,
                                      color = "grey30",
                                      margin = margin(b = 14)),
      axis.text.x      = element_text(size = fs$axis_x,
                                      face = "bold", color = "black",
                                      angle = 45, hjust = 1, vjust = 1),
      axis.text.y      = element_text(size = fs$axis_y, color = "black"),
      axis.ticks       = element_line(color = "black", linewidth = 0.5),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA,
                                      linewidth = 0.7),
      panel.spacing.x  = unit(fs$panel_spacing, "lines"),
      panel.spacing.y  = unit(fs$panel_spacing, "lines"),
      strip.placement  = "outside",
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = fs$legend_ttl),
      legend.text      = element_text(size = fs$legend_txt),
      legend.box       = "vertical",
      legend.spacing.y = unit(0.4, "cm"),
      plot.margin      = margin(15, 15, 15, 15)
    )
  
  list(plot = p, data = long, n_pathways = length(unique(long$Description)))
}

# ===== Layout A (basal, 4 contrasts) =================================
cat("\n===== Building Layout A (basal, no TNFα) =====\n")
result_basal <- make_curated_bubble(
  enrich_list      = go_bp_list,
  baseline_map     = baseline_map_basal,
  baseline_levels  = baseline_levels_basal,
  baseline_palette = baseline_palette_basal,
  title_main       = "GO BP - Stem cell + lineage curated pathways (basal, no TNFα)",
  title_sub        = paste0(nrow(curated_pathways),
                            " curated pathways across 4 categories; bubble size = gene count, color = -log10(p)")
)
p_basal <- result_basal$plot
cat("  Layout A: ", result_basal$n_pathways, " pathways\n", sep = "")

# ===== Layout B (full, 14 contrasts) =================================
cat("\n===== Building Layout B (full 14 contrasts) =====\n")
result_full <- make_curated_bubble(
  enrich_list      = go_bp_list,
  baseline_map     = baseline_map_full,
  baseline_levels  = baseline_levels_full,
  baseline_palette = baseline_palette_full,
  title_main       = "GO BP - Stem cell + lineage curated pathways (full 14 contrasts)",
  title_sub        = paste0(nrow(curated_pathways),
                            " curated pathways across 4 categories; bubble size = gene count, color = -log10(p)")
)
p_full <- result_full$plot
cat("  Layout B: ", result_full$n_pathways, " pathways\n", sep = "")

# =====================================================================
# 화면 출력 (한 번에 하나씩)
# =====================================================================
# ★ Layout A 보기
while (!is.null(dev.list())) dev.off()
print(p_basal)

# ★ Layout B 보기 (Layout A 확인 후 실행)
# while (!is.null(dev.list())) dev.off()
# print(p_full)