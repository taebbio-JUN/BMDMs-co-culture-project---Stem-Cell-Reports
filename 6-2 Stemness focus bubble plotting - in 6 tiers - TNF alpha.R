
# ===== 폰트 / 크기 — 튜닝 포인트 (여기만 만져서 조정) ================
bubble_params <- list(
  # 텍스트
  title        = 22,    # ★ 메인 타이틀
  subtitle     = 16,    # ★ 서브타이틀
  axis_x       = 16,    # ★ 하단 contrast 라벨
  axis_y       = 20,    # ★ 좌측 pathway 이름
  strip_x      = 13,    # ★ 상단 베이스라인 색 막대 두께 (글자는 투명)
  strip_y      = 14,    # ★ 우측 카테고리 라벨 (Stem cell core 등)
  legend_ttl   = 16,    # ★ 범례 제목
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



# =====================================================================
# Layout B (full, 14 contrasts) — 화면 출력 + 저장
# =====================================================================

# ===== 크기 지정 =====================================================
fig_w_full <- 16     # ★ 가로 (inches) — Full은 contrast 많아서 길게
fig_h_full <- 16     # ★ 세로 (inches)

# ===== 화면 출력 =====================================================
while (!is.null(dev.list())) dev.off()
print(p_full)

# ▼▼▼ 화면 확인 후 마음에 들면 아래 저장 코드 실행 ▼▼▼

# ===== 저장 (PNG / PDF / TIFF) =======================================
fname_full <- "Fig_S5f_GOBP_curated_full"

# PNG
ggsave(file.path(out_root, paste0(fname_full, ".png")),
       p_full, width = fig_w_full, height = fig_h_full,
       dpi = 300, bg = "white", limitsize = FALSE)

# PDF (vector)
ggsave(file.path(out_root, paste0(fname_full, ".pdf")),
       p_full, width = fig_w_full, height = fig_h_full,
       device = cairo_pdf, limitsize = FALSE)

# TIFF (출판용, 300 dpi + LZW)
ggsave(file.path(out_root, paste0(fname_full, ".tiff")),
       p_full, width = fig_w_full, height = fig_h_full,
       device = "tiff", dpi = 300, bg = "white",
       compression = "lzw", limitsize = FALSE)

cat("Layout B saved: ", fname_full, ".{png,pdf,tiff} (",
    fig_w_full, " × ", fig_h_full, " in)\n", sep = "")

# ===== Data table 저장 ===============================================
write.csv(result_full$data,
          file.path(out_root, "tables", "GOBP_curated_bubble_full.csv"),
          row.names = FALSE)