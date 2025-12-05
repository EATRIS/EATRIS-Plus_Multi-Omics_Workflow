DNAmAge_own <- function (x, clocks = "all", toBetas = FALSE, fastImp = FALSE, 
          normalize = FALSE, age, cell.count = TRUE, cell.count.reference = "blood gse35069 complete", cell.count.data=NULL,
          min.perc = 0.8, ...) 
{
  available.clocks <- c("Horvath", "Hannum", "Levine", 
                        "BNN", "skinHorvath", "PedBE", "Wu", 
                        "TL", "BLUP", "EN", "all")
  method <- match(clocks, available.clocks)
  if (any(is.na(method))) {
    stop("You wrote the name of an unavailable clock. Available clocks are:\n            Horvath, Hannum, Levine, BNN, skinHorvath, PedBE, Wu, TL, \n            BLUP and EN")
  }
  if (length(available.clocks) %in% method) {
    method <- seq_len(length(available.clocks) - 1)
  }
  cpgs <- methylclock:::getInputCpgValues(x, toBetas)
  #Quickfix
  #if (!all(available.clocks %in% ls(.GlobalEnv))) {
  #  load_DNAm_Clocks_data()
  #}
  if (inherits(x, "data.frame")) {
    cpgs.names <- as.character(x[, 1, drop = TRUE])
    if (length(grep("cg", cpgs.names)) == 0) {
      stop("First column should contain CpG names")
    }
    cpgs <- t(as.matrix(x[, -1]))
    colnames(cpgs) <- cpgs.names
  }
  else if (inherits(x, "matrix")) {
    cpgs <- t(x)
  }
  else if (inherits(x, "ExpressionSet")) {
    cpgs <- t(Biobase::exprs(x))
  }
  else if (inherits(x, "GenomicRatioSet")) {
    cpgs <- t(minfi::getBeta(x))
  }
  else {
    stop("x must be a data.frame or a 'GenomicRatioSet' or a \n             'ExpressionSet' object")
  }
  if (toBetas) {
    toBeta <- function(m) {
      2^m/(2^m + 1)
    }
    cpgs <- toBeta(cpgs)
  }
  if (any(cpgs < -0.1 | cpgs > 1.1, na.rm = TRUE)) {
    stop("Data seems to do not be beta values. Check your data or set 'toBetas=TRUE'")
  }
  cpgs.all <- c(coefHorvath$CpGmarker, coefHannum$CpGmarker, 
                coefLevine$CpGmarker, coefSkin$CpGmarker, coefPedBE$CpGmarker, 
                coefWu$CpGmarker, coefTL$CpGmarker, coefBLUP$CpGmarker, 
                coefEN$CpGmarker)
  cpgs.in <- intersect(cpgs.all, colnames(cpgs))
  miss <- apply(cpgs[, cpgs.in], 2, function(x) any(is.na(x)))
  if (any(miss)) {
    cpgs.imp <- cpgs_imputation(miss, cpgs, fastImp, cpgs.in=cpgs.in, ...)
  }
  else {
    cpgs.imp <- cpgs
  }
  if (1 %in% method) {
    DNAmAge <- methylclock:::predAge(cpgs.imp, coefHorvath, intercept = TRUE, 
                       min.perc)
    horvath <-  methylclock:::anti.trafo(DNAmAge)
    Horvath <- data.frame(id = rownames(cpgs.imp), Horvath = horvath)
  }
  if (2 %in% method) {
    hannum <- methylclock:::predAge(cpgs.imp, coefHannum, intercept = FALSE, 
                      min.perc)
    Hannum <- data.frame(id = rownames(cpgs.imp), Hannum = hannum)
  }
  if (3 %in% method) {
    levine <- methylclock:::predAge(cpgs.imp, coefLevine, intercept = TRUE, 
                      min.perc)
    Levine <- data.frame(id = rownames(cpgs.imp), Levine = levine)
  }
  if (4 %in% method) {
    if (any(!coefHorvath$CpGmarker[-1] %in% colnames(cpgs.imp))) {
      warning("Bayesian method cannot be estimated")
      bn <- rep(NA, nrow(cpgs.imp))
    }
    else {
      cpgs.bn <- t(cpgs.imp[, coefHorvath$CpGmarker[-1]])
      bn <- try(methylclock:::main_NewModel1Clean(cpgs.bn), TRUE)
      if (inherits(bn, "try-error")) {
        warning("Bayesian method produced an error")
        bn <- rep(NA, nrow(cpgs.imp))
      }
    }
    BNN <- data.frame(id = rownames(cpgs.imp), BNN = bn)
  }
  if (5 %in% method) {
    skinHorvath <- methylclock:::predAge(cpgs.imp, coefSkin, intercept = TRUE, 
                           min.perc)
    skinHorvath <- methylclock:::anti.trafo(skinHorvath)
    skinHorvath <- data.frame(id = rownames(cpgs.imp), skinHorvath = skinHorvath)
  }
  if (6 %in% method) {
    pedBE <- methylclock:::predAge(cpgs.imp, coefPedBE, intercept = TRUE, 
                     min.perc)
    pedBE <- methylclock:::anti.trafo(pedBE)
    PedBE <- data.frame(id = rownames(cpgs.imp), PedBE = pedBE)
  }
  if (7 %in% method) {
    wu <- methylclock:::predAge(cpgs.imp, coefWu, intercept = TRUE, min.perc)
    wu <- methylclock:::anti.trafo(wu)/12
    Wu <- data.frame(id = rownames(cpgs.imp), Wu = wu)
  }
  if (8 %in% method) {
    tl <- methylclock:::predAge(cpgs.imp, coefTL, intercept = TRUE, min.perc)
    TL <- data.frame(id = rownames(cpgs.imp), TL = tl)
  }
  if (9 %in% method | 10 %in% method) {
    colCpGs <- colnames(cpgs.imp)
    cpgs.imp <- t(apply(cpgs.imp, 1, scale))
    colnames(cpgs.imp) <- colCpGs
    if (9 %in% method) {
      blup <- methylclock:::predAge(cpgs.imp, coefBLUP, intercept = TRUE, 
                      min.perc)
      BLUP <- data.frame(id = rownames(cpgs.imp), BLUP = blup)
    }
    if (10 %in% method) {
      en <- methylclock:::predAge(cpgs.imp, coefEN, intercept = TRUE, 
                    min.perc)
      EN <- data.frame(id = rownames(cpgs.imp), EN = en)
    }
  }
  if (!missing(age)) {
    if (!cell.count) {
      if (1 %in% method) {
        Horvath <- methylclock:::ageAcc1(Horvath, age, lab = "Horvath")
      }
      if (2 %in% method) {
        Hannum <- methylclock:::ageAcc1(Hannum, age, lab = "Hannum")
      }
      if (3 %in% method) {
        Levine <- methylclock:::ageAcc1(Levine, age, lab = "Levine")
      }
      if (4 %in% method) {
        BNN <- methylclock:::ageAcc1(BNN, age, lab = "BNN")
      }
      if (5 %in% method) {
        skinHorvath <- methylclock:::ageAcc1(skinHorvath, age, lab = "skinHorvath")
      }
      if (6 %in% method) {
        PedBE <- methylclock:::ageAcc1(PedBE, age, lab = "PedBE")
      }
      if (7 %in% method) {
        Wu <- methylclock:::ageAcc1(Wu, age, lab = "Wu")
      }
      if (8 %in% method) {
        TL <- methylclock:::ageAcc1(TL, age, lab = "TL")
      }
      if (9 %in% method) {
        BLUP <- methylclock:::ageAcc1(BLUP, age, lab = "BLUP")
      }
      if (10 %in% method) {
        EN <- methylclock:::ageAcc1(EN, age, lab = "EN")
      }
    }
    else {
      if (is.null(cell.count.data)){
      cell.counts <- try(meffilEstimateCellCountsFromBetas(t(cpgs), 
                                                           cell.count.reference), TRUE)
      } else {
        cell.counts <- cell.count.data
      }
      if (inherits(cell.counts, "try-error")) {
        stop("cell counts cannot be estimated since\n                meffilEstimateCellCountsFromBetas function is giving an error.\n                Probably your data do not have any of the required CpGs for that\n                reference panel.")
      }
      else {
        ok <- which(apply(cell.counts, 2, IQR) > 1e-05)
        cell.counts <- cell.counts[, ok]
        df <- data.frame(age = age, cell.counts)
        if (1 %in% method) {
          Horvath <- methylclock:::ageAcc2(Horvath, df, lab = "Horvath")
        }
        if (2 %in% method) {
          Hannum <- methylclock:::ageAcc2(Hannum, df, lab = "Hannum")
        }
        if (3 %in% method) {
          Levine <- methylclock:::ageAcc2(Levine, df, lab = "Levine")
        }
        if (4 %in% method) {
          BNN <- methylclock:::ageAcc2(BNN, df, lab = "BNN")
        }
        if (5 %in% method) {
          skinHorvath <- methylclock:::ageAcc2(skinHorvath, df, lab = "Hovarth2")
        }
        if (6 %in% method) {
          PedBE <- methylclock:::ageAcc2(PedBE, df, lab = "PedBE")
        }
        if (7 %in% method) {
          Wu <- methylclock:::ageAcc2(Wu, df, lab = "Wu")
        }
        if (8 %in% method) {
          TL <- methylclock:::ageAcc2(TL, df, lab = "TL")
        }
        if (9 %in% method) {
          BLUP <- methylclock:::ageAcc2(BLUP, df, lab = "BLUP")
        }
        if (10 %in% method) {
          EN <- methylclock:::ageAcc2(EN, df, lab = "EN")
        }
      }
    }
  }
  else {
    cell.count <- FALSE
  }
  out <- NULL
  if (1 %in% method) {
    out <- Horvath
  }
  if (2 %in% method) {
    if (is.null(out)) {
      out <- Hannum
    }
    else {
      out <- out %>% full_join(Hannum, by = "id")
    }
  }
  if (3 %in% method) {
    if (is.null(out)) {
      out <- Levine
    }
    else {
      out <- out %>% full_join(Levine, by = "id")
    }
  }
  if (4 %in% method) {
    if (is.null(out)) {
      out <- BNN
    }
    else {
      out <- out %>% full_join(BNN, by = "id")
    }
  }
  if (5 %in% method) {
    if (is.null(out)) {
      out <- skinHorvath
    }
    else {
      out <- out %>% full_join(skinHorvath, by = "id")
    }
  }
  if (6 %in% method) {
    if (is.null(out)) {
      out <- PedBE
    }
    else {
      out <- out %>% full_join(PedBE, by = "id")
    }
  }
  if (7 %in% method) {
    if (is.null(out)) {
      out <- Wu
    }
    else {
      out <- out %>% full_join(Wu, by = "id")
    }
  }
  if (8 %in% method) {
    if (is.null(out)) {
      out <- TL
    }
    else {
      out <- out %>% full_join(TL, by = "id")
    }
  }
  if (9 %in% method) {
    if (is.null(out)) {
      out <- BLUP
    }
    else {
      out <- out %>% full_join(BLUP, by = "id")
    }
  }
  if (10 %in% method) {
    if (is.null(out)) {
      out <- EN
    }
    else {
      out <- out %>% full_join(EN, by = "id")
    }
  }
  out <- tibble::as_tibble(out)
  if (!missing(age)) {
    out <- add_column(out, age = age)
  }
  if (cell.count) {
    attr(out, "cell_proportion") <- cell.counts
  }
  out
}

cpgs_imputation <- function (imp, cpgs, fastImp, cpgs.in, ...) 
{
  cpgs.imp <- NULL
  if (any(imp)) {
    if (fastImp) {
      cat(paste("Imputing missing data of", sum(imp), 
                "CpGs .... \n"))
      mm <- apply(cpgs[, cpgs.in], 2, median, na.rm = TRUE)
      cpgs.imp <- sweep(cpgs[, cpgs.in], 2, STATS = mm, 
                        FUN = function(x, s) ifelse(is.na(x), s, x))
    }
    else {
      quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
      }
      cat(paste("Imputing missing data of the entire matrix .... \n"))
      cpgs.imp <- quiet(t(impute.knn(t(cpgs), ...)$data))
    }
    cat("Data imputed. Starting DNAm clock estimation ... \n")
  }
  else {
    cpgs.imp <- cpgs
  }
  return(cpgs.imp)
}
