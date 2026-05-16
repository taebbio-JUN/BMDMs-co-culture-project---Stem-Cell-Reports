# =====================================================================
# Project: HN00273522 — Bulk RNA-seq DE analysis
# Step 0: Inspect transcript count matrix structure
# =====================================================================

# ---- 0. Working directory & file path ----
# setwd("C:/your/project/path")   # 필요 시 수정
file_path <- "Transcript_count_matrix.csv"

# ---- 1. Load (base R first, no namespace conflicts) ----
counts_raw <- read.csv(
  file_path,
  header        = TRUE,
  row.names     = NULL,        # 일단 row.names 지정하지 않고 1열 확인
  check.names   = FALSE,       # 샘플명 자동 변형 방지 (e.g. "Sample-1" → "Sample.1")
  stringsAsFactors = FALSE
)

# ---- 2. Basic dimensions & column overview ----
cat("===== Dimensions =====\n")
cat("Rows (transcripts):", nrow(counts_raw), "\n")
cat("Cols (incl. ID col):", ncol(counts_raw), "\n\n")

cat("===== Column names =====\n")
print(colnames(counts_raw))
cat("\n")

cat("===== Column classes =====\n")
print(sapply(counts_raw, class))
cat("\n")

# ---- 3. Head / Tail / Random rows ----
cat("===== head(6) =====\n");  print(head(counts_raw, 6))
cat("\n===== tail(3) =====\n"); print(tail(counts_raw, 3))
cat("\n===== random 3 rows =====\n")
set.seed(1); print(counts_raw[sample(nrow(counts_raw), 3), ])

# ---- 4. ID column inspection (1st column 가정) ----
id_col <- counts_raw[[1]]
cat("\n===== ID column (col 1) =====\n")
cat("Class           :", class(id_col), "\n")
cat("Unique IDs      :", length(unique(id_col)), "/", length(id_col), "\n")
cat("Any duplicated  :", any(duplicated(id_col)), "\n")
cat("Any NA          :", any(is.na(id_col)), "\n")
cat("Example IDs     :\n"); print(head(id_col, 5))
# Ensembl(ENSMUST/ENSMUSG/ENST) vs RefSeq(NM_) vs Symbol 판별용

# ---- 5. Count matrix integrity check (assume cols 2:N are samples) ----
mat_part   <- counts_raw[, -1, drop = FALSE]
is_numeric <- sapply(mat_part, is.numeric)

cat("\n===== Sample columns =====\n")
cat("N sample cols   :", ncol(mat_part), "\n")
cat("All numeric?    :", all(is_numeric), "\n")
if (!all(is_numeric)) {
  cat("Non-numeric cols:\n"); print(names(mat_part)[!is_numeric])
}

# NA / negative / non-integer 검사 (DESeq2는 정수 카운트 요구)
cat("Any NA in counts:", any(is.na(mat_part)), "\n")
if (all(is_numeric)) {
  cat("Any negative    :", any(mat_part < 0, na.rm = TRUE), "\n")
  cat("All integer-like:", all(mat_part == floor(mat_part), na.rm = TRUE), "\n")
}

# ---- 6. Per-sample summary (library size & zero-inflation) ----
if (all(is_numeric)) {
  lib_size  <- colSums(mat_part, na.rm = TRUE)
  zero_frac <- colMeans(mat_part == 0, na.rm = TRUE)
  qc <- data.frame(
    sample          = names(mat_part),
    lib_size        = lib_size,
    mean_count      = round(colMeans(mat_part), 2),
    median_count    = apply(mat_part, 2, median),
    pct_zero        = round(zero_frac * 100, 2)
  )
  cat("\n===== Per-sample QC =====\n")
  print(qc, row.names = FALSE)
}

# ---- 7. Per-transcript summary (low-count rows) ----
if (all(is_numeric)) {
  rs <- rowSums(mat_part, na.rm = TRUE)
  cat("\n===== Per-transcript row sums =====\n")
  print(summary(rs))
  cat("Rows with sum = 0 :", sum(rs == 0), "\n")
  cat("Rows with sum < 10:", sum(rs < 10), "\n")
}

# ---- 8. Save snapshot for reference ----
str(counts_raw, list.len = 5)