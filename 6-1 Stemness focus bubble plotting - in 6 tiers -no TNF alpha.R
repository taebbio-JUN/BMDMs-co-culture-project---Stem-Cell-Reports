# =====================================================================
# 화면 출력 + Zoom 창 크기로 저장 (직접 크기 지정)
# =====================================================================

# ===== Zoom 창에서 보이는 크기 직접 입력 =============================
# Layout B (Full 14 contrasts) - 가로가 매우 김
fig_w_full <- 22    # ★ 직접 입력 (inches) - 풀스크린 Zoom 비율
fig_h_full <- 14    # ★ 직접 입력

# Layout A (Basal 4 contrasts) - 가로가 짧음
fig_w_basal <- 13
fig_h_basal <- 14

# ===== Layout A 화면 + 저장 ==========================================
while (!is.null(dev.list())) dev.off()
print(p_basal)

fname_basal <- "Fig_5f_GOBP_curated_basal"

ggsave(file.path(out_root, paste0(fname_basal, ".png")),
       p_basal, width = fig_w_basal, height = fig_h_basal,
       dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_root, paste0(fname_basal, ".pdf")),
       p_basal, width = fig_w_basal, height = fig_h_basal,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(out_root, paste0(fname_basal, ".tiff")),
       p_basal, width = fig_w_basal, height = fig_h_basal,
       device = "tiff", dpi = 300, bg = "white",
       compression = "lzw", limitsize = FALSE)

cat("Layout A saved (", fig_w_basal, " × ", fig_h_basal, " in)\n", sep = "")

# ===== Layout B 화면 + 저장 ==========================================
while (!is.null(dev.list())) dev.off()
print(p_full)

fname_full <- "Fig_S5f_GOBP_curated_full"

ggsave(file.path(out_root, paste0(fname_full, ".png")),
       p_full, width = fig_w_full, height = fig_h_full,
       dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_root, paste0(fname_full, ".pdf")),
       p_full, width = fig_w_full, height = fig_h_full,
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(out_root, paste0(fname_full, ".tiff")),
       p_full, width = fig_w_full, height = fig_h_full,
       device = "tiff", dpi = 300, bg = "white",
       compression = "lzw", limitsize = FALSE)

cat("Layout B saved (", fig_w_full, " × ", fig_h_full, " in)\n", sep = "")

# ===== Data tables ===================================================
write.csv(result_basal$data,
          file.path(out_root, "tables", "GOBP_curated_bubble_basal.csv"),
          row.names = FALSE)
write.csv(result_full$data,
          file.path(out_root, "tables", "GOBP_curated_bubble_full.csv"),
          row.names = FALSE)