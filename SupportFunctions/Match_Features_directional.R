match_features_bidirectional <- function(df1,
                                         df2,
                                         mz1_col = "mz",
                                         rt1_col = "RT",
                                         mz2_col = "mz",
                                         rt2_col = "RT..min.",
                                         mz_tol_ppm = 5,
                                         rt_tol_min = 0.1) {
  # Helper: convert with decimal commas to numeric
  as_num <- function(x) {
    x <- gsub(",", ".", as.character(x))
    as.numeric(x)
  }
  
  # Column checks
  for (nm in c(mz1_col, rt1_col)) {
    if (!nm %in% names(df1)) {
      stop("df1 must contain column: ", nm)
    }
  }
  for (nm in c(mz2_col, rt2_col)) {
    if (!nm %in% names(df2)) {
      stop("df2 must contain column: ", nm)
    }
  }
  
  # Extract mz/RT and convert
  mz1 <- as_num(df1[[mz1_col]])
  rt1 <- as_num(df1[[rt1_col]])
  mz2 <- as_num(df2[[mz2_col]])
  rt2 <- as_num(df2[[rt2_col]])
  
  # Keep only rows with non-NA mz & RT
  keep1 <- stats::complete.cases(mz1, rt1)
  keep2 <- stats::complete.cases(mz2, rt2)
  
  mz1v <- mz1[keep1]
  rt1v <- rt1[keep1]
  mz2v <- mz2[keep2]
  rt2v <- rt2[keep2]
  
  if (!length(mz1v) || !length(mz2v)) {
    message("No non-NA mz/RT values to compare.")
    return(data.frame())
  }
  
  n1 <- length(mz1v)
  n2 <- length(mz2v)
  
  # Precompute pairwise differences
  mz_diff  <- abs(outer(mz1v, mz2v, `-`))    # n1 x n2
  rt_diff  <- abs(outer(rt1v, rt2v, `-`))    # n1 x n2
  
  # ppm relative to df1 (rows)
  mz_den1   <- matrix(mz1v, nrow = n1, ncol = n2)
  ppm_diff1 <- mz_diff / mz_den1 * 1e6
  
  # ppm relative to df2 (columns)
  mz_den2   <- matrix(mz2v, nrow = n1, ncol = n2, byrow = TRUE)
  ppm_diff2 <- mz_diff / mz_den2 * 1e6
  
  orig1_idx <- which(keep1)   # original row numbers in df1
  orig2_idx <- which(keep2)   # original row numbers in df2
  
  ## ---- STEP 1: df1 -> df2 (df1 as reference) ----
  dir1_list <- lapply(seq_len(n1), function(i) {
    cand <- which(ppm_diff1[i, ] <= mz_tol_ppm & rt_diff[i, ] <= rt_tol_min)
    if (!length(cand)) return(NULL)
    # best match = minimal ppm error (relative to df1)
    j <- cand[which.min(ppm_diff1[i, cand])]
    data.frame(
      df1_row   = orig1_idx[i],
      df2_row   = orig2_idx[j],
      mz1       = mz1v[i],
      RT1       = rt1v[i],
      mz2       = mz2v[j],
      RT2       = rt2v[j],
      ppm_error = ppm_diff1[i, j],
      RT_diff   = rt_diff[i, j],
      direction = "1->2",
      stringsAsFactors = FALSE
    )
  })
  dir1 <- do.call(rbind, Filter(Negate(is.null), dir1_list))
  
  ## ---- STEP 2: df2 -> df1 (df2 as reference) ----
  dir2_list <- lapply(seq_len(n2), function(j) {
    cand <- which(ppm_diff2[, j] <= mz_tol_ppm & rt_diff[, j] <= rt_tol_min)
    if (!length(cand)) return(NULL)
    # best match = minimal ppm error (relative to df2)
    i <- cand[which.min(ppm_diff2[cand, j])]
    data.frame(
      df1_row   = orig1_idx[i],
      df2_row   = orig2_idx[j],
      mz1       = mz1v[i],
      RT1       = rt1v[i],
      mz2       = mz2v[j],
      RT2       = rt2v[j],
      ppm_error = ppm_diff2[i, j],
      RT_diff   = rt_diff[i, j],
      direction = "2->1",
      stringsAsFactors = FALSE
    )
  })
  dir2 <- do.call(rbind, Filter(Negate(is.null), dir2_list))
  
  ## ---- Combine and make unique pairs ----
  if (is.null(dir1) && is.null(dir2)) {
    message("No matches found within the specified tolerances.")
    return(data.frame())
  }
  
  all_matches <- rbind(
    if (!is.null(dir1)) dir1 else NULL,
    if (!is.null(dir2)) dir2 else NULL
  )
  
  # uniqueness based on (df1_row, df2_row)
  key <- paste(all_matches$df1_row, all_matches$df2_row, sep = "_")
  unique_idx <- !duplicated(key)
  result <- all_matches[unique_idx, ]
  
  # Sort by df1_row, then df2_row
  result <- result[order(result$df1_row, result$df2_row), ]
  
  print(result)
  invisible(result)
}
