# Cardinal NEST - analysis module
# Univariate (meansTest/contrastTest), multivariate (PCA/NMF), supervised
# (PLS/OPLS-DA, SSC), clustering (spatialKMeans), colocalization. Class-based
# methods share a grouping mechanism (pixelData column or coord/run split);
# results export to CSV, PNG/PDF, or back onto the dataset.

library(shiny)
library(ggplot2)

# Muffle the benign matter coercion warning when file-backed spectra are
# pulled into a local matrix; real warnings still surface.
with_muffled_coercion_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("coercing input to a local matrix", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Estimated in-memory spectra size: n_features * n_pixels * 8 bytes.
estimate_spectra_bytes <- function(data) {
  as.numeric(nrow(data)) * as.numeric(ncol(data)) * 8
}

# Materialize matter spectra to a base matrix in column chunks (a single
# as.matrix() overflows matter's 32-bit size computation).
ensure_in_memory <- function(data, na_action = "zero", impute_stat = "mean") {
  s <- tryCatch(Cardinal::spectra(data), error = function(e) NULL)
  if (is.null(s)) {
    return(data)
  }
  if (is.matrix(s) && is.numeric(s) && !anyNA(s)) {
    # Already a clean in-memory numeric matrix; leave the object alone.
    return(data)
  }
  # Preserve NAs; the chosen missing-value policy applies below.
  Cardinal::spectra(data) <- materialize_matrix(s, na_action = "leave")
  apply_na_handling(data, na_action, impute_stat)
}

# Column-block materialization staying under matter's ~268M-element
# 32-bit overflow point; small in-memory matrices pass through.
materialize_matrix <- function(x, max_elems = 1e8, na_action = "zero") {
  if (is.matrix(x) && is.numeric(x) && !anyNA(x)) {
    # Already a clean in-memory numeric matrix; nothing to do.
    return(x)
  }
  if (is.matrix(x)) {
    out <- as.matrix(x)
  } else {
    nr <- nrow(x)
    nc <- ncol(x)
    if (is.null(nr) || is.null(nc)) {
      out <- as.matrix(x)
    } else {
      out <- matrix(0, nrow = nr, ncol = nc)
      block_ncol <- max(1L, as.integer(floor(max_elems / max(nr, 1))))
      start <- 1L
      while (start <= nc) {
        end <- min(start + block_ncol - 1L, nc)
        cols <- start:end
        out[, cols] <- as.matrix(x[, cols, drop = FALSE])
        start <- end + 1L
      }
    }
  }
  if (!is.numeric(out)) {
    stop("Failed to materialize spectra into a numeric matrix.", call. = FALSE)
  }
  if (anyNA(out) && identical(na_action, "zero")) {
    # Sparse gaps materialize as NA; in MSI missing intensity means zero
    # signal, and Cardinal's engines fail on NAs anyway.
    warning("Spectra matrix contained NAs; imputing missing values with 0.",
            call. = FALSE)
    out[is.na(out)] <- 0
  }
  # na_action == "leave": caller applies its own policy via apply_na_handling().
  out
}

# NA policy on materialized spectra: "leave", "zero", "remove_features",
# "remove_pixels", or "impute" (feature mean/median, fallback 0).
apply_na_handling <- function(data,
                              method = c("zero", "leave", "remove_features",
                                         "remove_pixels", "impute"),
                              impute_stat = c("mean", "median")) {
  method <- match.arg(method)
  impute_stat <- match.arg(impute_stat)
  if (identical(method, "leave")) {
    return(data)
  }
  s <- tryCatch(Cardinal::spectra(data), error = function(e) NULL)
  if (is.null(s) || !anyNA(s)) {
    return(data)
  }
  s <- as.matrix(s)

  if (identical(method, "zero")) {
    s[is.na(s)] <- 0
    Cardinal::spectra(data) <- s
    return(data)
  }

  if (identical(method, "impute")) {
    if (identical(impute_stat, "mean")) {
      imputed <- rowMeans(s, na.rm = TRUE)
    } else {
      imputed <- apply(s, 1, stats::median, na.rm = TRUE)
    }
    imputed[is.nan(imputed)] <- 0  # features with no observed values -> 0
    na_idx <- which(is.na(s), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
      s[na_idx] <- imputed[na_idx[, 1]]
    }
    Cardinal::spectra(data) <- s
    return(data)
  }

  if (identical(method, "remove_features")) {
    keep <- !apply(s, 1, anyNA)
    if (!any(keep)) {
      stop("Every feature contains missing values; none remain after removal.",
           call. = FALSE)
    }
    return(data[keep, ])
  }

  keep <- !apply(s, 2, anyNA)
  if (!any(keep)) {
    stop("Every pixel contains missing values; none remain after removal.",
         call. = FALSE)
  }
  data[, keep]
}

# pixelData columns usable as a class variable (factor/char, or numeric
# with <= max_levels distinct values; needs >= 2 levels).
group_candidate_columns <- function(data, max_levels = 20) {
  pd <- as.data.frame(Cardinal::pixelData(data))
  keep <- character()
  for (nm in names(pd)) {
    v <- pd[[nm]]
    nu <- length(unique(v[!is.na(v)]))
    is_categorical <- is.factor(v) || is.character(v) || is.logical(v)
    if ((is_categorical || nu <= max_levels) && nu >= 2) {
      keep <- c(keep, nm)
    }
  }
  keep
}

# Grouping factor: coord splits `axis` at `threshold` (default median),
# run gives one group per acquisition run.
make_group_labels <- function(data, method = c("coord", "run"),
                              axis = "x", threshold = NULL) {
  method <- match.arg(method)
  if (method == "run") {
    return(factor(as.character(Cardinal::run(data))))
  }
  cd <- as.data.frame(Cardinal::coord(data))
  if (!axis %in% names(cd)) {
    stop(sprintf("Coordinate axis '%s' not found.", axis), call. = FALSE)
  }
  vals <- cd[[axis]]
  if (is.null(threshold) || is.na(threshold)) {
    threshold <- stats::median(vals, na.rm = TRUE)
  }
  lo <- paste0(axis, " <= ", signif(threshold, 4))
  hi <- paste0(axis, " > ", signif(threshold, 4))
  factor(ifelse(vals <= threshold, lo, hi), levels = c(lo, hi))
}

# Cardinal::meansTest across `class_labels`; `sample_labels` aggregate
# pixels to sample means first (avoids pseudo-replication), NULL = each
# pixel its own sample. FDR recomputed when a non-BH method is requested.
run_means_test <- function(data, class_labels, sample_labels = NULL,
                           fdr_method = "BH") {
  class_labels <- as.factor(class_labels)
  if (nlevels(class_labels) < 2) {
    stop("The class variable must have at least 2 groups.", call. = FALSE)
  }
  if (length(class_labels) != ncol(data)) {
    stop("Class labels must have one value per pixel.", call. = FALSE)
  }

  d <- data
  Cardinal::pData(d)$.nest_class <- class_labels

  if (is.null(sample_labels)) {
    samples <- factor(seq_len(ncol(d)))
  } else {
    if (length(sample_labels) != ncol(data)) {
      stop("Sample labels must have one value per pixel.", call. = FALSE)
    }
    samples <- as.factor(sample_labels)
  }

  mt <- Cardinal::meansTest(d, fixed = ~ .nest_class, samples = samples)
  res <- as.data.frame(Cardinal::topFeatures(mt, n = Inf, sort.by = "statistic"))

  # topFeatures uses p.adjust("fdr"); recompute for other methods
  # ("fdr" and "BH" are synonymous in p.adjust).
  if (is.null(fdr_method) || identical(fdr_method, "none")) {
    res$fdr <- res$pvalue
  } else if (!identical(fdr_method, "fdr") && !identical(fdr_method, "BH")) {
    res$fdr <- stats::p.adjust(res$pvalue, method = fdr_method)
  }

  res <- res[order(res$pvalue), , drop = FALSE]
  # Keep the model so a contrast test can reuse it.
  attr(res, "means_model") <- mt
  res
}

# Cardinal::contrastTest on a meansTest model; `method` = contrast type,
# `emm_adjust` = emmeans-internal adjustment, FDR recomputed per contrast
# for non-BH methods. A compatible `model` is reused, else fit here and
# returned as attr(res, "means_model").
run_contrast_test <- function(data, class_labels, sample_labels = NULL,
                              method = "pairwise", emm_adjust = "none",
                              fdr_method = "BH", model = NULL) {
  class_labels <- as.factor(class_labels)
  if (nlevels(class_labels) < 2) {
    stop("The class variable must have at least 2 groups.", call. = FALSE)
  }
  if (length(class_labels) != ncol(data)) {
    stop("Class labels must have one value per pixel.", call. = FALSE)
  }

  if (is.null(model)) {
    d <- data
    Cardinal::pData(d)$.nest_class <- class_labels

    if (is.null(sample_labels)) {
      samples <- factor(seq_len(ncol(d)))
    } else {
      if (length(sample_labels) != ncol(data)) {
        stop("Sample labels must have one value per pixel.", call. = FALSE)
      }
      samples <- as.factor(sample_labels)
    }

    mt <- Cardinal::meansTest(d, fixed = ~ .nest_class, samples = samples)
  } else {
    mt <- model
  }

  ct <- Cardinal::contrastTest(mt, specs = ~ .nest_class,
                               method = method, emm_adjust = emm_adjust)
  res <- as.data.frame(Cardinal::topFeatures(ct, n = Inf))

  # topFeatures uses p.adjust("fdr"); recompute for other methods.
  pval_cols <- grep("\\.pvalue$", names(res), value = TRUE)
  if (is.null(fdr_method) || identical(fdr_method, "none")) {
    for (pv in pval_cols) {
      fdr_col <- sub("\\.pvalue$", ".fdr", pv)
      res[[fdr_col]] <- res[[pv]]
    }
  } else if (!identical(fdr_method, "fdr") && !identical(fdr_method, "BH")) {
    for (pv in pval_cols) {
      fdr_col <- sub("\\.pvalue$", ".fdr", pv)
      res[[fdr_col]] <- stats::p.adjust(res[[pv]], method = fdr_method)
    }
  }

  # Sort by min p-value across contrasts so any-significant features
  # surface first.
  if (length(pval_cols) > 0) {
    min_p <- do.call(pmin, c(res[pval_cols], list(na.rm = TRUE)))
    res <- res[order(min_p), , drop = FALSE]
  }
  # Keep the model so later contrast runs can reuse it.
  attr(res, "means_model") <- mt
  res
}

# Cardinal::PCA -> list(model, scores, loadings, var_explained, coord).
run_pca_analysis <- function(data, ncomp = 3) {
  ncomp <- as.integer(ncomp)
  if (is.na(ncomp) || ncomp < 2) ncomp <- 2
  # irlba requires ncomp < min(n_features, n_pixels).
  max_comp <- min(nrow(data), ncol(data)) - 1L
  if (max_comp < 2L) {
    stop("Dataset is too small for PCA (need at least 3 features and 3 pixels).",
         call. = FALSE)
  }
  ncomp <- min(ncomp, max_comp)

  pc <- Cardinal::PCA(data, ncomp = ncomp)

  scores <- as.data.frame(pc$x)
  loadings <- as.data.frame(pc$rotation)
  loadings$mz <- Cardinal::mz(data)

  ve <- pc$sdev^2 / sum(pc$sdev^2)
  names(ve) <- colnames(pc$x)

  list(
    model = pc,
    scores = scores,
    loadings = loadings,
    var_explained = ve,
    coord = as.data.frame(Cardinal::coord(data))
  )
}

# Top features by absolute loading for a given principal component.
top_loadings <- function(loadings, pc = "PC1", n = 10) {
  if (!pc %in% names(loadings)) {
    stop(sprintf("Component '%s' not found in loadings.", pc), call. = FALSE)
  }
  ord <- order(abs(loadings[[pc]]), decreasing = TRUE)
  out <- loadings[ord, c("mz", pc), drop = FALSE]
  utils::head(out, n)
}

# Dataset record with PCA scores as pixelData columns (PCA_*).
make_pca_record <- function(data, scores, id, filename = NULL, ext = NULL) {
  cols <- stats::setNames(as.list(scores), paste0("PCA_", colnames(scores)))
  make_result_record(data, cols, id, paste0("PCA (", ncol(scores), " PCs)"),
                     filename = filename, ext = ext)
}

# Cardinal::NMF; per-component contribution approximated from score-column
# sums of squares (no eigenvalue-based variance measure).
run_nmf_analysis <- function(data, ncomp = 3) {
  ncomp <- as.integer(ncomp)
  if (is.na(ncomp) || ncomp < 2) ncomp <- 2
  # Rank must be < min(features, pixels), as with PCA.
  max_comp <- min(nrow(data), ncol(data)) - 1L
  if (max_comp < 2L) {
    stop("Dataset is too small for NMF (need at least 3 features and 3 pixels).",
         call. = FALSE)
  }
  ncomp <- min(ncomp, max_comp)

  nmf <- Cardinal::NMF(data, ncomp = ncomp)

  scores <- as.data.frame(nmf$x)
  loadings <- as.data.frame(nmf$activation)
  loadings$mz <- Cardinal::mz(data)

  ss <- colSums(as.matrix(nmf$x)^2)
  ve <- if (sum(ss) > 0) ss / sum(ss) else rep(1 / length(ss), length(ss))
  names(ve) <- colnames(nmf$x)

  list(
    model = nmf,
    scores = scores,
    loadings = loadings,
    var_explained = ve,
    coord = as.data.frame(Cardinal::coord(data))
  )
}

# Dataset record with NMF scores as pixelData columns (NMF_*).
make_nmf_record <- function(data, scores, id, filename = NULL, ext = NULL) {
  cols <- stats::setNames(as.list(scores), paste0("NMF_", colnames(scores)))
  make_result_record(data, cols, id, paste0("NMF (", ncol(scores), " comps)"),
                     filename = filename, ext = ext)
}

# PLS-DA or OPLS-DA via Cardinal::PLS/OPLS; class factor as response,
# predicted class = max fitted response. var_explained only for PLS.
run_discriminant_analysis <- function(data, class_labels, ncomp = 3,
                                       method = c("pls", "opls")) {
  method <- match.arg(method)
  class_labels <- as.factor(class_labels)
  if (nlevels(class_labels) < 2) {
    stop("The class variable must have at least 2 groups.", call. = FALSE)
  }
  if (length(class_labels) != ncol(data)) {
    stop("Class labels must have one value per pixel.", call. = FALSE)
  }

  ncomp <- as.integer(ncomp)
  if (is.na(ncomp) || ncomp < 1) ncomp <- 1
  # Latent components cannot exceed min(n_features, n_pixels) - 1.
  max_comp <- min(nrow(data), ncol(data)) - 1L
  if (max_comp < 1L) {
    stop("Dataset is too small for discriminant analysis (need at least 2 features and 2 pixels).",
         call. = FALSE)
  }
  ncomp <- min(ncomp, max_comp)

  fit_fn <- if (method == "pls") Cardinal::PLS else Cardinal::OPLS
  # Muffle only the benign matter "local matrix" coercion warning.
  model <- with_muffled_coercion_warning(
    fit_fn(data, class_labels, ncomp = ncomp)
  )

  comp_names <- colnames(model[["scores"]])
  scores <- as.data.frame(model[["scores"]])
  names(scores) <- comp_names

  loadings <- as.data.frame(model[["loadings"]])
  names(loadings) <- comp_names
  loadings$mz <- Cardinal::mz(data)

  # Per-component captured variance (PLS only).
  ve <- NULL
  if ("cvar" %in% names(model) && !is.null(model[["cvar"]])) {
    cv <- as.numeric(model[["cvar"]])
    if (isTRUE(sum(cv) > 0)) {
      ve <- cv / sum(cv)
      names(ve) <- comp_names
    }
  }

  fitted_vals <- as.matrix(fitted(model))
  predicted <- factor(colnames(fitted_vals)[max.col(fitted_vals, ties.method = "first")],
                      levels = levels(class_labels))
  confusion <- table(Predicted = predicted, Actual = class_labels)
  accuracy <- sum(predicted == class_labels) / length(class_labels)

  list(
    model = model,
    method = method,
    scores = scores,
    loadings = loadings,
    var_explained = ve,
    fitted = fitted_vals,
    predicted = predicted,
    actual = class_labels,
    accuracy = accuracy,
    confusion = confusion,
    coord = as.data.frame(Cardinal::coord(data))
  )
}

# Dataset record with DA scores + predicted class (PLS_/OPLS_*).
make_da_record <- function(data, da, id, filename = NULL, ext = NULL) {
  prefix <- toupper(da$method)
  cols <- stats::setNames(as.list(da$scores), paste0(prefix, "_", names(da$scores)))
  cols[[paste0(prefix, "_predicted")]] <- da$predicted
  make_result_record(data, cols, id,
                     paste0(prefix, "-DA (", ncol(da$scores), " comps)"),
                     filename = filename, ext = ext)
}

# Supervised spatialShrunkenCentroids; class factor as `y`.
# r = radius, s = shrinkage, weights = "gaussian"|"adaptive".
run_ssc_classification <- function(data, class_labels, r = 2, s = 9,
                                   weights = c("gaussian", "adaptive")) {
  weights <- match.arg(weights)
  class_labels <- as.factor(class_labels)
  if (nlevels(class_labels) < 2) {
    stop("The class variable must have at least 2 groups.", call. = FALSE)
  }
  if (length(class_labels) != ncol(data)) {
    stop("Class labels must have one value per pixel.", call. = FALSE)
  }

  if (is.null(r) || is.na(r) || r < 0) r <- 1
  if (is.null(s) || is.na(s) || s < 0) s <- 0
  r <- as.numeric(r)
  s <- as.numeric(s)

  # Muffle only the benign matter "local matrix" coercion warning.
  model <- with_muffled_coercion_warning(
    Cardinal::spatialShrunkenCentroids(data, y = class_labels, r = r, s = s,
                                       weights = weights)
  )
  # A vector s returns a list of models; use the first.
  if (is.null(model$class) && is.list(model) && length(model) >= 1L) {
    model <- model[[1L]]
  }

  predicted <- factor(as.character(model$class), levels = levels(class_labels))
  confusion <- table(Predicted = predicted, Actual = class_labels)
  accuracy <- sum(predicted == class_labels) / length(class_labels)

  list(
    model = model,
    r = r,
    s = s,
    weights = weights,
    predicted = predicted,
    actual = class_labels,
    accuracy = accuracy,
    confusion = confusion,
    statistic = model$statistic,
    mz = Cardinal::mz(data),
    coord = as.data.frame(Cardinal::coord(data))
  )
}

# Top features per class via topFeatures().
ssc_top_features <- function(ssc, n = 20) {
  if (is.null(ssc$model)) {
    return(data.frame(i = integer(0), mz = numeric(0),
                      class = character(0), statistic = numeric(0),
                      centers = numeric(0), sd = numeric(0)))
  }
  as.data.frame(Cardinal::topFeatures(ssc$model, n = n, sort.by = "statistic"))
}

# Canonicalize a CV metric matrix to metrics-as-rows x classes-as-columns;
# PLS and SSC return transposed orientations, so locate "Recall" to decide.
.cv_canonical_metric_matrix <- function(mat) {
  if (is.null(mat) || !is.matrix(mat)) return(mat)
  if ("Recall" %in% rownames(mat)) {
    mat
  } else if ("Recall" %in% colnames(mat)) {
    t(mat)
  } else {
    mat
  }
}

# Cardinal::crossValidate for PLS/SSC; `seed` makes the stratified folds
# reproducible. Spectra are fit as-is, so scores may be mildly optimistic
# vs a fully nested pipeline.
run_cross_validation <- function(data, class_labels, model_type = c("pls", "ssc"),
                                 ncomp = 2, r = 2, s = 9,
                                 weights = "gaussian", nfolds = 5, seed = NULL) {
  model_type <- match.arg(model_type)
  class_labels <- as.factor(class_labels)
  if (nlevels(class_labels) < 2) {
    stop("The class variable must have at least 2 groups.", call. = FALSE)
  }
  if (length(class_labels) != ncol(data)) {
    stop("Class labels must have one value per pixel.", call. = FALSE)
  }

  nfolds <- as.integer(nfolds)
  if (is.na(nfolds) || nfolds < 2) nfolds <- 5L
  nfolds <- min(nfolds, ncol(data))

  if (!is.null(seed) && !is.na(seed)) set.seed(as.integer(seed))

  # Stratified folds: shuffle within class, deal round-robin so every fold
  # contains members of every class.
  folds <- integer(ncol(data))
  for (cl in levels(class_labels)) {
    idx <- which(class_labels == cl)
    idx <- sample(idx)
    folds[idx] <- rep(seq_len(nfolds), length.out = length(idx))
  }

  fit_fn <- if (model_type == "pls") Cardinal::PLS else Cardinal::spatialShrunkenCentroids

  # trainProcess/testProcess = NULL: use the data as-is between folds.
  extra_args <- if (model_type == "pls") {
    list(ncomp = ncomp)
  } else {
    list(r = r, s = s, weights = weights, k = nlevels(class_labels))
  }

  cv <- with_muffled_coercion_warning(
    do.call(Cardinal::crossValidate,
            c(list(fit_fn, data, y = class_labels, folds = folds,
                   trainProcess = NULL, testProcess = NULL),
              extra_args))
  )

  # cv's `$`/`[[` are overridden; unclass() for plain list access.
  m <- unclass(cv@model)
  if (is.null(m[["average"]]) || !is.matrix(m[["average"]])) {
    stop("Cross-validation returned an unexpected model structure (no average matrix).",
         call. = FALSE)
  }
  cls_levels <- levels(class_labels)

  avg <- .cv_canonical_metric_matrix(m[["average"]])

  avg_df <- as.data.frame(avg)
  avg_df$metric <- rownames(avg)
  rownames(avg_df) <- NULL

  # Per-fold score matrices, canonicalized to metrics x class.
  fold_scores <- m[["scores"]]
  if (is.list(fold_scores)) {
    fold_scores <- lapply(fold_scores, .cv_canonical_metric_matrix)
  }

  # Mean per-class recall across folds (not overall accuracy).
  mean_recall <- NA_real_
  if (is.list(fold_scores) && length(fold_scores) > 0) {
    recalls <- vapply(fold_scores, function(fs) {
      if (is.null(fs) || !is.matrix(fs) || all(is.na(fs))) return(NA_real_)
      if (!"Recall" %in% rownames(fs)) return(NA_real_)
      mean(fs["Recall", ], na.rm = TRUE)
    }, numeric(1))
    if (any(!is.na(recalls))) mean_recall <- mean(recalls[!is.na(recalls)])
  }

  list(
    cv = cv,
    model_type = model_type,
    nfolds = nfolds,
    average = avg,
    average_df = avg_df,
    fold_scores = fold_scores,
    folds = folds,
    mean_recall = mean_recall,
    class_levels = cls_levels
  )
}

# Per-class metric bars from the CV average matrix.
plot_cross_validation <- function(cv_result) {
  avg <- cv_result[["average"]]
  cls <- cv_result[["class_levels"]]
  metric_names <- rownames(avg)
  df <- data.frame(
    class = rep(cls, each = length(metric_names)),
    metric = rep(metric_names, length(cls)),
    value = as.numeric(avg)
  )
  ggplot(df, aes(x = class, y = value, fill = metric)) +
    geom_col(position = "dodge", na.rm = TRUE) +
    scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
    labs(title = sprintf("Cross-validation (%d folds, %s)",
                         cv_result[["nfolds"]], toupper(cv_result[["model_type"]])),
         x = "Class", y = "Score", fill = "Metric") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
}

# Dataset record with SSC predictions as pixelData column ssc_class.
make_ssc_record <- function(data, ssc, id, filename = NULL, ext = NULL) {
  make_result_record(data, list(ssc_class = ssc$predicted), id,
                     sprintf("SSC classification (r = %g, s = %g)",
                             ssc$r, ssc$s),
                     filename = filename, ext = ext)
}

# Cardinal::spatialKMeans clustering into k groups (radius r).
run_clustering <- function(data, k = 3, r = 1) {
  k <- as.numeric(k)
  r <- as.numeric(r)
  if (is.na(k) || k < 2) k <- 2
  if (ncol(data) < 2L) {
    stop("Dataset has too few pixels for clustering (need at least 2).",
         call. = FALSE)
  }
  k <- min(k, ncol(data))

  coord <- as.data.frame(Cardinal::coord(data))

  model <- Cardinal::spatialKMeans(data, r = r, k = k)
  cluster <- as.integer(model[["cluster"]])
  extra <- list(r = r,
                tot_withinss = model[["tot.withinss"]],
                betweenss = model[["betweenss"]])

  cluster <- factor(cluster, levels = sort(unique(cluster)))

  list(
    method = "kmeans",
    model = model,
    k = nlevels(cluster),
    cluster = cluster,
    sizes = table(cluster),
    coord = coord,
    extra = extra
  )
}

# Dataset record with clusters as pixelData column kmeans_cluster.
make_clustering_record <- function(data, cl, id, filename = NULL, ext = NULL) {
  cols <- stats::setNames(list(cl$cluster), paste0(cl$method, "_cluster"))
  make_result_record(data, cols, id,
                     sprintf("k-means clustering (k = %d)", cl$k),
                     filename = filename, ext = ext)
}

# Rank features by colocalization with a reference m/z
# (Cardinal::colocalized) or a 0/1 region mask (per-feature Pearson
# correlation). Exactly one of `mz`/`mask` required.
run_colocalization <- function(data, mz = NULL, mask = NULL, n = 20) {
  n_feat <- nrow(data)
  has_mz <- !is.null(mz) && !is.na(mz)
  has_mask <- !is.null(mask)
  if (has_mz && has_mask) {
    stop("Provide either a reference m/z or a region, not both.", call. = FALSE)
  }
  if (has_mz) {
    mz <- as.numeric(mz)
    # colocalized() errors when the m/z isn't a feature; snap to nearest.
    feat_mz <- Cardinal::mz(data)
    if (length(feat_mz) > 0) {
      mz <- feat_mz[which.min(abs(feat_mz - mz))]
    }
    coloc <- Cardinal::colocalized(data, mz = mz, n = n_feat)
    df <- as.data.frame(coloc)
    ref_label <- sprintf("m/z %s", signif(mz, 6))
  } else if (has_mask) {
    mask <- as.logical(mask)
    if (length(mask) != ncol(data)) {
      stop("The region mask must have one value per pixel.", call. = FALSE)
    }
    if (!any(mask, na.rm = TRUE)) {
      stop("The selected region contains no pixels.", call. = FALSE)
    }
    # Column-block materialization; direct as.matrix() overflows on big data.
    intensities <- materialize_matrix(Cardinal::spectra(data))  # features x pixels
    refv <- as.numeric(mask)
    cor_vals <- suppressWarnings(
      apply(intensities, 1, function(x) stats::cor(x, refv))
    )
    cor_vals[!is.finite(cor_vals)] <- 0
    df <- data.frame(mz = Cardinal::mz(data), correlation = cor_vals)
    df <- df[order(-df$correlation), , drop = FALSE]
    rownames(df) <- NULL
    ref_label <- "region"
    mz <- NA_real_
  } else {
    stop("Provide either a reference m/z or a region to colocalize against.",
         call. = FALSE)
  }

  n <- as.integer(n)
  if (is.na(n) || n < 1) n <- 1L
  top_mz <- if ("mz" %in% names(df)) utils::head(df$mz, n) else numeric(0)

  list(
    ref_label = ref_label,
    ref_mz = mz,
    table = df,
    top_mz = top_mz
  )
}

# -log10(p) vs m/z with an FDR <= alpha threshold line and highlights.
plot_univariate <- function(res, alpha = 0.05) {
  neglp <- -log10(res$pvalue)
  neglp[!is.finite(neglp)] <- max(neglp[is.finite(neglp)], 0, na.rm = TRUE)
  sig <- !is.na(res$fdr) & res$fdr <= alpha
  cols <- ifelse(sig, "#d7301f", "#4d4d4d")
  plot(res$mz, neglp, pch = 19, col = cols,
       xlab = "m/z", ylab = expression(-log[10](p)),
       main = "Univariate significance")
  thr <- suppressWarnings(min(neglp[sig], na.rm = TRUE))
  if (is.finite(thr)) {
    abline(h = thr, lty = 2, col = "#d7301f")
  }
  legend("topright", legend = c(sprintf("FDR <= %.3g", alpha), "n.s."),
         col = c("#d7301f", "#4d4d4d"), pch = 19, bty = "n")
  invisible(TRUE)
}

# Same significance plot on the contrast test's <contrast>.pvalue/.fdr
# column naming (first contrast).
plot_contrast <- function(res, alpha = 0.05) {
  pval_col <- grep("\\.pvalue$", names(res), value = TRUE)
  if (length(pval_col) > 0) {
    pval_col <- pval_col[1]
    fdr_col <- sub("\\.pvalue$", ".fdr", pval_col)
  } else if ("pvalue" %in% names(res)) {
    pval_col <- "pvalue"
    fdr_col <- "fdr"
  } else {
    plot.new(); text(0.5, 0.5, "No contrast p-values to display.")
    return(invisible(FALSE))
  }
  p <- res[[pval_col]]
  f <- res[[fdr_col]]
  neglp <- -log10(p)
  neglp[!is.finite(neglp)] <- max(neglp[is.finite(neglp)], 0, na.rm = TRUE)
  sig <- !is.na(f) & f <= alpha
  cols <- ifelse(sig, "#d7301f", "#4d4d4d")
  plot(res$mz, neglp, pch = 19, col = cols,
       xlab = "m/z", ylab = expression(-log[10](p)),
       main = "Contrast significance")
  thr <- suppressWarnings(min(neglp[sig], na.rm = TRUE))
  if (is.finite(thr)) {
    abline(h = thr, lty = 2, col = "#d7301f")
  }
  legend("topright", legend = c(sprintf("FDR <= %.3g", alpha), "n.s."),
         col = c("#d7301f", "#4d4d4d"), pch = 19, bty = "n")
  invisible(TRUE)
}

# Decomposition scores scatter (shared by PCA/NMF); `ve` annotates axes.
plot_decomp_scores <- function(scores, ve, comp_x, comp_y, groups = NULL,
                               main = "Scores") {
  axis_lab <- function(cn) {
    if (!is.null(ve) && cn %in% names(ve)) {
      sprintf("%s (%.1f%%)", cn, 100 * ve[[cn]])
    } else {
      cn
    }
  }
  if (!is.null(groups)) {
    groups <- as.factor(groups)
    palette <- grDevices::hcl.colors(max(nlevels(groups), 2), "Dark 3")
    cols <- palette[as.integer(groups)]
  } else {
    cols <- "#377eb8"
  }
  plot(scores[[comp_x]], scores[[comp_y]], pch = 19, col = cols,
       xlab = axis_lab(comp_x), ylab = axis_lab(comp_y), main = main)
  abline(h = 0, v = 0, lty = 3, col = "grey70")
  if (!is.null(groups)) {
    legend("topright", legend = levels(groups),
           col = grDevices::hcl.colors(max(nlevels(groups), 2), "Dark 3"),
           pch = 19, bty = "n")
  }
  invisible(TRUE)
}

# Generic per-component contribution bar chart (shared by PCA/NMF).
plot_decomp_variance <- function(ve, main, ylab) {
  bp <- barplot(100 * ve, names.arg = names(ve),
                ylab = ylab, main = main, col = "#377eb8")
  invisible(bp)
}

# Generic loadings/activation stick plot vs m/z for one component.
plot_decomp_loadings <- function(loadings, comp, main) {
  if (!comp %in% names(loadings)) comp <- setdiff(names(loadings), "mz")[1]
  plot(loadings$mz, loadings[[comp]], type = "h", col = "#377eb8",
       xlab = "m/z", ylab = sprintf("Loading (%s)", comp), main = main)
  abline(h = 0, lty = 3, col = "grey70")
  invisible(TRUE)
}

# Spatial score map via image(model, type = "x"); NULL comp shows all.
plot_decomp_image <- function(model, scores, comp = NULL) {
  idx <- if (is.null(comp)) {
    seq_len(ncol(scores))
  } else {
    m <- match(comp, colnames(scores))
    if (is.na(m)) 1L else m
  }
  print(Cardinal::image(model, type = "x", i = idx, scale = TRUE))
  invisible(TRUE)
}

# PCA scores scatter for two components, optionally colored by group labels.
plot_pca_scores <- function(pca, pc_x = "PC1", pc_y = "PC2", groups = NULL) {
  plot_decomp_scores(pca$scores, pca$var_explained, pc_x, pc_y,
                     groups = groups, main = "PCA scores")
}

# Variance-explained bar chart for the PCA model.
plot_pca_variance <- function(pca) {
  plot_decomp_variance(pca$var_explained, "PCA variance explained",
                       "Variance explained (%)")
}

# PCA loadings stick plot for one principal component.
plot_pca_loadings <- function(pca, comp = "PC1") {
  plot_decomp_loadings(pca$loadings, comp, sprintf("PCA loadings (%s)", comp))
}

# PCA spatial score map (image(pca, type = "x")).
plot_pca_image <- function(pca, comp = NULL) {
  plot_decomp_image(pca$model, pca$scores, comp)
}

# NMF scores scatter for two components, optionally colored by group labels.
plot_nmf_scores <- function(nmf, comp_x = "C1", comp_y = "C2", groups = NULL) {
  plot_decomp_scores(nmf$scores, nmf$var_explained, comp_x, comp_y,
                     groups = groups, main = "NMF scores")
}

# Component-contribution bar chart for the NMF model.
plot_nmf_variance <- function(nmf) {
  plot_decomp_variance(nmf$var_explained, "NMF component contribution",
                       "Relative contribution (%)")
}

# NMF activation (loadings) stick plot for one component.
plot_nmf_loadings <- function(nmf, comp = "C1") {
  plot_decomp_loadings(nmf$loadings, comp, sprintf("NMF activation (%s)", comp))
}

# NMF spatial score map (image(nmf, type = "x")).
plot_nmf_image <- function(nmf, comp = NULL) {
  plot_decomp_image(nmf$model, nmf$scores, comp)
}

# DA scores scatter colored by class; stripchart fallback for one
# component. Axes show captured variance % when known (PLS).
plot_da_scores <- function(da, comp_x = "C1", comp_y = "C2", groups = NULL) {
  scores <- da$scores
  comps <- names(scores)
  if (!comp_x %in% comps) comp_x <- comps[1]
  if (!comp_y %in% comps) comp_y <- comps[min(2, length(comps))]
  if (is.null(groups)) groups <- da$actual
  groups <- as.factor(groups)

  ve <- da$var_explained
  axis_lab <- function(cn) {
    if (!is.null(ve) && cn %in% names(ve)) {
      sprintf("%s (%.1f%%)", cn, 100 * ve[[cn]])
    } else {
      cn
    }
  }

  palette <- grDevices::hcl.colors(max(nlevels(groups), 2), "Dark 3")
  main <- sprintf("%s-DA scores", toupper(da$method))

  if (length(comps) < 2L) {
    # One component: show its distribution across classes.
    stripchart(scores[[comp_x]] ~ groups, vertical = TRUE, method = "jitter",
               pch = 19, col = palette,
               ylab = axis_lab(comp_x), xlab = "Class", main = main)
  } else {
    plot(scores[[comp_x]], scores[[comp_y]], pch = 19,
         col = palette[as.integer(groups)],
         xlab = axis_lab(comp_x), ylab = axis_lab(comp_y), main = main)
    abline(h = 0, v = 0, lty = 3, col = "grey70")
  }
  legend("topright", legend = levels(groups), col = palette,
         pch = 19, bty = "n")
  invisible(TRUE)
}

# Fitted class-membership response map via image(model, "response").
plot_da_response <- function(da) {
  print(Cardinal::image(da$model, type = "response", scale = TRUE))
  invisible(TRUE)
}

# Coefficients identifying the features driving class prediction.
plot_da_coefficients <- function(da) {
  print(Cardinal::plot(da$model, type = "coefficients", linewidth = 2))
  invisible(TRUE)
}

# Predicted-class map; 1-D strip fallback.
plot_ssc_class <- function(ssc) {
  coord <- ssc$coord
  classes <- as.factor(ssc$predicted)
  n <- nlevels(classes)
  palette <- grDevices::hcl.colors(max(n, 2), "Dark 3")
  main <- sprintf("SSC predicted class (r = %g, s = %g)", ssc$r, ssc$s)

  if (all(c("x", "y") %in% names(coord))) {
    df <- data.frame(coord, class = classes)
    ggplot(df, aes(x = x, y = y, fill = class)) +
      geom_tile() +
      coord_equal() +
      scale_fill_manual(values = palette, name = "Class") +
      labs(title = main, x = "x", y = "y") +
      theme_minimal()
  } else {
    df <- data.frame(pixel = seq_along(classes), class = classes)
    ggplot(df, aes(x = pixel, y = 1, fill = class)) +
      geom_tile() +
      scale_fill_manual(values = palette, name = "Class") +
      labs(title = main, x = "pixel", y = "") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# Per-class fitted membership probability map.
plot_ssc_probability <- function(ssc) {
  print(Cardinal::image(ssc$model, type = "probability"))
  invisible(TRUE)
}

# Shrunken t-statistics identifying class-distinguishing features.
plot_ssc_statistic <- function(ssc) {
  print(Cardinal::plot(ssc$model, type = "statistic", linewidth = 2))
  invisible(TRUE)
}

# Cluster map; 1-D strip fallback.
plot_clustering <- function(cl) {
  coord <- cl$coord
  clusters <- cl$cluster
  n <- nlevels(clusters)
  palette <- grDevices::hcl.colors(max(n, 2), "Dark 3")
  label <- "Spatial k-means"
  main <- sprintf("%s clustering (k = %d)", label, n)

  if (all(c("x", "y") %in% names(coord))) {
    df <- data.frame(coord, cluster = clusters)
    ggplot(df, aes(x = x, y = y, fill = cluster)) +
      geom_tile() +
      coord_equal() +
      scale_fill_manual(values = palette, name = "Cluster") +
      labs(title = main, x = "x", y = "y") +
      theme_minimal()
  } else {
    df <- data.frame(pixel = seq_along(clusters), cluster = clusters)
    ggplot(df, aes(x = pixel, y = 1, fill = cluster)) +
      geom_tile() +
      scale_fill_manual(values = palette, name = "Cluster") +
      labs(title = main, x = "pixel", y = "") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# Ion images of the top co-localized features.
plot_colocalization <- function(coloc, data, top = 3) {
  mzs <- utils::head(coloc$top_mz, top)
  if (length(mzs) == 0L) {
    plot.new(); text(0.5, 0.5, "No co-localized features to display.")
    return(invisible(FALSE))
  }
  print(Cardinal::image(data, mz = mzs, scale = TRUE))
  invisible(TRUE)
}

analysisUI <- function(id) {
  ns <- NS(id)

  tagList(
    shinyjs::useShinyjs(),
    fluidRow(
      column(12,
             h3("Statistical Analysis"),
             p("Run univariate (mean tests), multivariate (PCA/NMF), supervised, clustering, or colocalization analysis on the active dataset."),
             uiOutput(ns("active_dataset"))
      )
    ),

    fluidRow(
      column(4,
        bslib::card(
          bslib::card_header("Grouping"),
          bslib::card_body(
            radioButtons(ns("group_source"), "Class variable source",
                         choices = c("Existing metadata column" = "existing",
                                     "Create from coordinates/run" = "create"),
                         selected = "existing"),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'existing'", ns("group_source")),
              selectInput(ns("group_col"), "Class column", choices = character(0))
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'create'", ns("group_source")),
              selectInput(ns("create_method"), "Method",
                          choices = c("Coordinate threshold" = "coord",
                                      "Run name" = "run")),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'coord'", ns("create_method")),
                selectInput(ns("create_axis"), "Axis",
                            choices = c("x" = "x", "y" = "y")),
                numericInput(ns("create_threshold"),
                             "Threshold (blank = median)", value = NA, step = "any")
              )
            ),

            selectInput(ns("sample_col"),
                        "Sample/replicate column (optional)",
                        choices = c("None (per pixel)" = ""),
                        selectize = FALSE)
          )
        ),

        bslib::card(
          bslib::card_header("Analysis"),
          bslib::card_body(
            selectInput(ns("analysis_type"), "Analysis",
                        choices = c("Univariate (mean test)" = "univariate",
                                    "Contrast test (post-hoc)" = "contrast",
                                    "Multivariate (PCA)" = "pca",
                                    "Multivariate (NMF)" = "nmf",
                                    "PLS-DA (supervised)" = "plsda",
                                    "OPLS-DA (supervised)" = "oplsda",
                                    "SSC classification (supervised)" = "sscda",
                                    "k-means clustering" = "kmeans",
                                    "Colocalization" = "coloc")),

            selectInput(ns("handle_na"), "Handle NAs",
                        choices = c("Leave as missing" = "leave",
                                    "Replace with zero" = "zero",
                                    "Remove affected features" = "remove_features",
                                    "Remove affected pixels" = "remove_pixels",
                                    "Impute values" = "impute"),
                        selected = "zero"),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'impute'", ns("handle_na")),
              selectInput(ns("impute_stat"), "Impute statistic",
                          choices = c("Mean" = "mean",
                                      "Median" = "median"),
                          selected = "mean")
            ),

            div(class = "text-muted small",
                "How missing intensities are treated before analysis. 'Impute values' uses each feature's selected statistic (mean or median) of observed intensities. 'Leave as missing' keeps NAs, which some analyses may not support."),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'univariate'", ns("analysis_type")),
              selectInput(ns("fdr_method"), "Multiple-testing correction",
                          choices = c("Benjamini-Hochberg (FDR)" = "BH",
                                      "Bonferroni" = "bonferroni",
                                      "None" = "none")),
              numericInput(ns("alpha"), "Significance threshold (FDR)",
                           value = 0.05, min = 0, max = 1, step = 0.01)
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'contrast'", ns("analysis_type")),
              selectInput(ns("contrast_method"), "Contrast method",
                          choices = c("Pairwise" = "pairwise",
                                      "Reverse pairwise" = "revpairwise",
                                      "Dunnett (vs first)" = "dunnett",
                                      "Tukey" = "tukey")),
              selectInput(ns("contrast_emm_adjust"), "emmeans p-value adjustment",
                          choices = c("None" = "none",
                                      "Tukey" = "tukey",
                                      "Bonferroni" = "bonferroni",
                                      "FDR (BH)" = "fdr")),
              selectInput(ns("contrast_fdr_method"), "Feature-level FDR correction",
                          choices = c("Benjamini-Hochberg (FDR)" = "BH",
                                      "Bonferroni" = "bonferroni",
                                      "None" = "none")),
              numericInput(ns("contrast_alpha"), "Significance threshold (FDR)",
                           value = 0.05, min = 0, max = 1, step = 0.01),
              div(class = "text-muted small",
                  "Post-hoc contrasts on the mean-test models via Cardinal::contrastTest (emmeans).")
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'pca'", ns("analysis_type")),
              numericInput(ns("ncomp"), "Number of components",
                           value = 3, min = 2, step = 1)
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'nmf'", ns("analysis_type")),
              numericInput(ns("nmf_ncomp"), "Number of components",
                           value = 3, min = 2, step = 1),
              div(class = "text-muted small",
                  "Non-negative matrix factorization (Cardinal::NMF).")
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'plsda' || input['%s'] == 'oplsda'",
                                  ns("analysis_type"), ns("analysis_type")),
              numericInput(ns("da_ncomp"), "Number of components",
                           value = 3, min = 1, step = 1),
              checkboxInput(ns("da_cv"), "Run cross-validation", value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s']", ns("da_cv")),
                numericInput(ns("da_nfolds"), "Number of folds",
                             value = 5, min = 2, step = 1)
              ),
              div(class = "text-muted small",
                  "Supervised: requires the class variable defined above.")
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'sscda'", ns("analysis_type")),
              numericInput(ns("ssc_r"), "Spatial neighborhood radius (r)",
                           value = 2, min = 0, step = 1),
              numericInput(ns("ssc_s"), "Shrinkage / sparsity (s)",
                           value = 9, min = 0, step = 1),
              selectInput(ns("ssc_weights"), "Spatial weights",
                          choices = c("Gaussian" = "gaussian",
                                      "Adaptive" = "adaptive")),
              checkboxInput(ns("ssc_cv"), "Run cross-validation", value = FALSE),
              conditionalPanel(
                condition = sprintf("input['%s']", ns("ssc_cv")),
                numericInput(ns("ssc_nfolds"), "Number of folds",
                             value = 5, min = 2, step = 1)
              ),
              div(class = "text-muted small",
                  "Supervised spatial shrunken centroids (Cardinal::spatialShrunkenCentroids with y = class); requires the class variable defined above.")
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'kmeans'", ns("analysis_type")),
              numericInput(ns("clust_k"), "Number of clusters (k)",
                           value = 3, min = 2, step = 1),
              numericInput(ns("clust_r"), "Spatial neighborhood radius (r)",
                           value = 1, min = 0, step = 1),
              div(class = "text-muted small",
                  "Spatially-aware k-means (Cardinal::spatialKMeans).")
            ),

            conditionalPanel(
              condition = sprintf("input['%s'] == 'coloc'", ns("analysis_type")),
              radioButtons(ns("coloc_ref"), "Colocalize with",
                           choices = c("A reference m/z" = "mz",
                                       "A region (class level)" = "region"),
                           selected = "mz"),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'mz'", ns("coloc_ref")),
                numericInput(ns("coloc_mz"), "Reference m/z",
                             value = NA, step = "any")
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'region'", ns("coloc_ref")),
                selectInput(ns("coloc_level"), "Region (class level)",
                            choices = character(0))
              ),
              numericInput(ns("coloc_top"), "Ion images to display",
                           value = 3, min = 1, max = 9, step = 1),
              div(class = "text-muted small",
                  "Finds features co-localized with the chosen ion or region (Cardinal::colocalized).")
            ),

            actionButton(ns("run_analysis"), "Run analysis",
                         icon = icon("play"), class = "btn-success w-100")
          )
        )
      ),

      column(8,
        bslib::card(
          bslib::card_header("Status"),
          bslib::card_body(uiOutput(ns("status")))
        ),
        bslib::card(
          bslib::card_header(
            div(class = "d-flex justify-content-between align-items-center",
              span("Results"),
              uiOutput(ns("result_controls"))
            )
          ),
          bslib::card_body(
            bslib::navset_tab(
              bslib::nav_panel("Plot",
                shinycssloaders::withSpinner(plotOutput(ns("result_plot"), height = "60vh"))
              ),
              bslib::nav_panel("Table",
                DT::DTOutput(ns("result_table"))
              )
            )
          )
        ),
        bslib::card(
          bslib::card_header("Export"),
          bslib::card_body(
            fluidRow(
              column(4, downloadButton(ns("export_csv"), "Download table (.csv)",
                                       class = "btn-outline-primary w-100")),
              column(4,
                selectInput(ns("plot_format"), NULL,
                            choices = c("PNG" = "png", "PDF" = "pdf"), width = "100%"),
                downloadButton(ns("export_plot"), "Download plot",
                               class = "btn-outline-primary w-100")
              ),
              column(4, uiOutput(ns("save_scores_ui")))
            )
          )
        )
      )
    )
  )
}

analysisServer <- function(id, datasets, active_id,
                           perf_log = NULL, perf_state = NULL,
                           analysis_cache = NULL) {
  moduleServer(id, function(input, output, session) {
    library(Cardinal)
    ns <- session$ns

    # Inputs that never affect the result (run/export/save/view pickers).
    volatile_inputs <- c("run_analysis", "plot_format", "export_csv",
                         "export_plot", "save_scores", "save_scores_confirm",
                         "save_filename", "save_ext")
    volatile_suffix <- "(_view|_x|_y)$"

    # Cache key: dataset fingerprint + every non-volatile input.
    analysis_signature <- function(data) {
      keys <- setdiff(names(reactiveValuesToList(input)), volatile_inputs)
      keys <- sort(keys[!grepl(volatile_suffix, keys)])
      vals <- vapply(keys, function(k) {
        v <- input[[k]]
        if (is.null(v)) "" else paste(as.character(v), collapse = ",")
      }, character(1))
      paste(dataset_fingerprint(data), active_id(),
            paste(keys, vals, sep = "=", collapse = "|"), sep = "||")
    }

    # Signature of dataset + grouping + sample inputs; decides whether a
    # persisted meansTest model can be reused by a contrast test.
    model_signature <- function(data) {
      keys <- c("group_source", "group_col", "create_method", "create_axis",
                "create_threshold", "sample_col")
      vals <- vapply(keys, function(k) {
        v <- input[[k]]
        if (is.null(v)) "" else paste(as.character(v), collapse = ",")
      }, character(1))
      paste(dataset_fingerprint(data), active_id(),
            paste(vals, collapse = "|"), sep = "||")
    }

    # Record an elapsed-time row into the shared perf log when available.
    record_perf <- function(operation, seconds, note = "") {
      if (is.null(perf_log)) return(invisible())
      perf_log(perf_log_append(perf_log(), "analysis", operation, seconds, note))
    }

    mse_data <- active_dataset(datasets, active_id)

    status <- reactiveVal(list(type = "info", message = "No analysis run yet."))
    result <- reactiveVal(NULL)  # list(kind, table, plot_fn, ...)
    # Persisted meansTest fit + signature for contrast-test reuse.
    means_model <- reactiveVal(NULL)

    output$active_dataset <- active_dataset_ui(datasets, active_id)

    observe({
      data <- mse_data()
      if (is.null(data)) {
        updateSelectInput(session, "group_col", choices = character(0))
        updateSelectInput(session, "sample_col", choices = c("None (per pixel)" = ""))
        return()
      }
      cols <- group_candidate_columns(data)
      updateSelectInput(session, "group_col", choices = cols)

      pd_cols <- names(as.data.frame(Cardinal::pixelData(data)))
      updateSelectInput(session, "sample_col",
                        choices = c("None (per pixel)" = "", stats::setNames(pd_cols, pd_cols)))
    })

    current_group_labels <- reactive({
      data <- mse_data()
      if (is.null(data)) return(NULL)

      if (identical(input$group_source, "existing")) {
        req(input$group_col)
        pd <- as.data.frame(Cardinal::pixelData(data))
        if (!input$group_col %in% names(pd)) return(NULL)
        as.factor(pd[[input$group_col]])
      } else {
        method <- input$create_method
        if (identical(method, "run")) {
          make_group_labels(data, method = "run")
        } else {
          make_group_labels(data, method = "coord",
                            axis = input$create_axis,
                            threshold = input$create_threshold)
        }
      }
    })

    current_sample_labels <- reactive({
      data <- mse_data()
      if (is.null(data)) return(NULL)
      if (is.null(input$sample_col) || !nzchar(input$sample_col)) return(NULL)
      pd <- as.data.frame(Cardinal::pixelData(data))
      if (!input$sample_col %in% names(pd)) return(NULL)
      as.factor(pd[[input$sample_col]])
    })

    # Coloc region-level choices follow the current grouping.
    observe({
      labs <- tryCatch(current_group_labels(), error = function(e) NULL)
      if (is.null(labs)) {
        updateSelectInput(session, "coloc_level", choices = character(0))
        return()
      }
      lv <- levels(as.factor(labs))
      updateSelectInput(session, "coloc_level", choices = lv)
    })

    observe({
      if (is.null(mse_data())) {
        shinyjs::disable("run_analysis")
      } else {
        shinyjs::enable("run_analysis")
      }
    })

    observeEvent(input$run_analysis, {
      data <- mse_data()
      if (is.null(data)) {
        status(list(type = "danger", message = "No active dataset selected."))
        return()
      }

      # Engines fail on file-backed matter storage; materialize once here.
      na_action <- if (is.null(input$handle_na)) "zero" else input$handle_na
      impute_stat <- if (is.null(input$impute_stat)) "mean" else input$impute_stat
      data <- tryCatch(
        with_progress_stages("Loading spectra into memory for analysis",
                     ensure_in_memory(data, na_action = na_action,
                                      impute_stat = impute_stat)),
        error = function(e) {
          status(list(type = "danger",
                      message = sprintf("Could not load spectra into memory (%.1f GB estimated): %s",
                                        estimate_spectra_bytes(data) / 1024^3,
                                        conditionMessage(e))))
          NULL
        }
      )
      if (is.null(data)) return()

      # Reuse a cached result for an identical signature when caching is on.
      cache_on <- !is.null(analysis_cache) &&
        (is.null(perf_state) || isTRUE(perf_state$cache_enabled))
      sig <- if (cache_on) analysis_signature(data) else NULL
      if (cache_on) {
        hit <- analysis_cache$get(sig)
        if (!is.null(hit)) {
          result(hit$result)
          status(hit$status)
          record_perf(input$analysis_type, 0, note = "cached")
          return()
        }
      }
      run_t0 <- proc.time()[["elapsed"]]

      if (identical(input$analysis_type, "univariate")) {
        labs <- tryCatch(current_group_labels(), error = function(e) NULL)
        if (is.null(labs)) {
          status(list(type = "danger",
                      message = "Could not resolve a class variable for the test."))
          return()
        }
        alpha <- input$alpha
        res <- with_progress_stages("Running mean test", {
          tryCatch(
            run_means_test(data, labs,
                           sample_labels = current_sample_labels(),
                           fdr_method = input$fdr_method),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(res)) return()

        # Persist the model so a matching contrast test can reuse it.
        means_model(list(model = attr(res, "means_model"),
                         sig = model_signature(data)))

        n_sig <- sum(!is.na(res$fdr) & res$fdr <= alpha)
        result(list(
          kind = "univariate",
          table = res,
          plot_fn = function() plot_univariate(res, alpha = alpha)
        ))
        status(list(type = "success",
                    message = sprintf("Mean test complete: %d of %d features significant (FDR <= %.3g).",
                                      n_sig, nrow(res), alpha)))

      } else if (identical(input$analysis_type, "contrast")) {
        labs <- tryCatch(current_group_labels(), error = function(e) NULL)
        if (is.null(labs)) {
          status(list(type = "danger",
                      message = "Could not resolve a class variable for the contrast test."))
          return()
        }
        alpha <- input$contrast_alpha
        # Reuse the persisted meansTest fit when the signature matches.
        mm <- means_model()
        reuse_model <- if (!is.null(mm) &&
                           identical(mm$sig, model_signature(data))) {
          mm$model
        } else NULL
        res <- with_progress_stages("Running contrast test", {
          tryCatch(
            run_contrast_test(data, labs,
                              sample_labels = current_sample_labels(),
                              method = input$contrast_method,
                              emm_adjust = input$contrast_emm_adjust,
                              fdr_method = input$contrast_fdr_method,
                              model = reuse_model),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(res)) return()

        means_model(list(model = attr(res, "means_model"),
                         sig = model_signature(data)))

        fdr_cols <- grep("\\.fdr$", names(res), value = TRUE)
        n_sig <- if (length(fdr_cols) > 0) sum(!is.na(res[[fdr_cols[1]]]) & res[[fdr_cols[1]]] <= alpha) else 0
        result(list(
          kind = "contrast",
          table = res,
          plot_fn = function() plot_contrast(res, alpha = alpha)
        ))
        status(list(type = "success",
                    message = sprintf("Contrast test complete: %d of %d features significant (FDR <= %.3g).",
                                      n_sig, nrow(res), alpha)))

      } else if (identical(input$analysis_type, "pca")) {
        pca <- with_progress_stages("Running PCA", {
          tryCatch(
            run_pca_analysis(data, ncomp = input$ncomp),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(pca)) return()

        labs <- tryCatch(current_group_labels(), error = function(e) NULL)
        comps <- colnames(pca$scores)
        tbl <- cbind(pca$coord, pca$scores)
        result(list(
          kind = "pca",
          table = tbl,
          pca = pca,
          groups = labs,
          plot_fn = function() {
            plot_pca_scores(pca, pc_x = comps[1], pc_y = comps[min(2, length(comps))],
                            groups = labs)
          }
        ))
        ve <- pca$var_explained
        status(list(type = "success",
                    message = sprintf("PCA complete: %d components, %.1f%% variance explained.",
                                      length(ve), 100 * sum(ve))))

      } else if (identical(input$analysis_type, "nmf")) {
        nmf <- with_progress_stages("Running NMF", {
          tryCatch(
            run_nmf_analysis(data, ncomp = input$nmf_ncomp),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(nmf)) return()

        labs <- tryCatch(current_group_labels(), error = function(e) NULL)
        comps <- colnames(nmf$scores)
        tbl <- cbind(nmf$coord, nmf$scores)
        result(list(
          kind = "nmf",
          table = tbl,
          nmf = nmf,
          groups = labs,
          plot_fn = function() {
            plot_nmf_scores(nmf, comp_x = comps[1], comp_y = comps[min(2, length(comps))],
                            groups = labs)
          }
        ))
        status(list(type = "success",
                    message = sprintf("NMF complete: %d components extracted over %d features.",
                                      length(comps), nrow(nmf$loadings))))

      } else if (input$analysis_type %in% c("plsda", "oplsda")) {
        method <- if (identical(input$analysis_type, "oplsda")) "opls" else "pls"
        labs <- tryCatch(current_group_labels(), error = function(e) NULL)
        if (is.null(labs)) {
          status(list(type = "danger",
                      message = "Could not resolve a class variable for the discriminant analysis."))
          return()
        }
        da <- with_progress_stages(sprintf("Running %s-DA", toupper(method)), {
          tryCatch(
            run_discriminant_analysis(data, labs, ncomp = input$da_ncomp,
                                      method = method),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(da)) return()

        comps <- names(da$scores)
        # Coords + scores + actual/predicted class.
        tbl <- cbind(da$coord, da$scores,
                     Actual = as.character(da$actual),
                     Predicted = as.character(da$predicted))
        result(list(
          kind = "da",
          method = method,
          table = tbl,
          da = da,
          groups = da$actual,
          plot_fn = function() {
            plot_da_scores(da, comp_x = comps[1],
                           comp_y = comps[min(2, length(comps))],
                           groups = da$actual)
          }
        ))
        ve_msg <- if (!is.null(da$var_explained)) {
          sprintf(", %.1f%% variance captured", 100 * sum(da$var_explained))
        } else {
          ""
        }
        status(list(type = "success",
                    message = sprintf("%s-DA complete: %d components, %.1f%% training accuracy%s.",
                                      toupper(method), length(comps),
                                      100 * da$accuracy, ve_msg)))

        # Optional cross-validation.
        if (isTRUE(input$da_cv)) {
          cv <- with_progress_stages("Running cross-validation", {
            tryCatch(
              run_cross_validation(data, labs, model_type = "pls",
                                   ncomp = input$da_ncomp,
                                   nfolds = input$da_nfolds),
              error = function(e) {
                status(list(type = "danger", message = conditionMessage(e)))
                NULL
              }
            )
          })
          if (!is.null(cv)) {
            cv_msg <- sprintf(" Cross-validation: %.1f%% mean recall (%d folds).",
                              100 * cv[["mean_recall"]], cv[["nfolds"]])
            res <- result()
            res[["cv"]] <- cv
            res$plot_fn <- function() {
              p1 <- plot_da_scores(da, comp_x = comps[1],
                                   comp_y = comps[min(2, length(comps))],
                                   groups = da$actual)
              print(p1)
              grid::grid.newpage()
              grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
              print(p1, vp = grid::viewport(layout.pos.row = 1))
              print(plot_cross_validation(cv),
                    vp = grid::viewport(layout.pos.row = 2))
            }
            result(res)
            status(list(type = "success",
                        message = sprintf("%s-DA complete: %d components, %.1f%% training accuracy%s.%s",
                                          toupper(method), length(comps),
                                          100 * da$accuracy, ve_msg, cv_msg)))
          }
        }

      } else if (identical(input$analysis_type, "sscda")) {
        labs <- tryCatch(current_group_labels(), error = function(e) NULL)
        if (is.null(labs)) {
          status(list(type = "danger",
                      message = "Could not resolve a class variable for the SSC classification."))
          return()
        }
        ssc <- with_progress_stages("Running SSC classification", {
          tryCatch(
            run_ssc_classification(
              data, labs,
              r = if (is.null(input$ssc_r)) 2 else input$ssc_r,
              s = if (is.null(input$ssc_s)) 9 else input$ssc_s,
              weights = if (is.null(input$ssc_weights)) "gaussian" else input$ssc_weights),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(ssc)) return()

        # Coords + actual/predicted class.
        tbl <- cbind(ssc$coord,
                     Actual = as.character(ssc$actual),
                     Predicted = as.character(ssc$predicted))
        result(list(
          kind = "sscda",
          table = tbl,
          ssc = ssc,
          groups = ssc$actual,
          plot_fn = function() plot_ssc_class(ssc)
        ))
        status(list(type = "success",
                    message = sprintf("SSC classification complete: r = %g, s = %g, %.1f%% training accuracy.",
                                      ssc$r, ssc$s, 100 * ssc$accuracy)))

        # Optional cross-validation.
        if (isTRUE(input$ssc_cv)) {
          cv <- with_progress_stages("Running cross-validation", {
            tryCatch(
              run_cross_validation(data, labs, model_type = "ssc",
                                   r = ssc$r, s = ssc$s,
                                   weights = ssc$weights,
                                   nfolds = input$ssc_nfolds),
              error = function(e) {
                status(list(type = "danger", message = conditionMessage(e)))
                NULL
              }
            )
          })
          if (!is.null(cv)) {
            cv_msg <- sprintf(" Cross-validation: %.1f%% mean recall (%d folds).",
                              100 * cv[["mean_recall"]], cv[["nfolds"]])
            res <- result()
            res[["cv"]] <- cv
            res$plot_fn <- function() {
              print(plot_ssc_class(ssc))
              grid::grid.newpage()
              grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
              print(plot_ssc_class(ssc),
                    vp = grid::viewport(layout.pos.row = 1))
              print(plot_cross_validation(cv),
                    vp = grid::viewport(layout.pos.row = 2))
            }
            result(res)
            status(list(type = "success",
                        message = sprintf("SSC classification complete: r = %g, s = %g, %.1f%% training accuracy.%s",
                                          ssc$r, ssc$s, 100 * ssc$accuracy, cv_msg)))
          }
        }

      } else if (identical(input$analysis_type, "coloc")) {
        mz_ref <- NULL
        mask <- NULL
        if (identical(input$coloc_ref, "region")) {
          labs <- tryCatch(current_group_labels(), error = function(e) NULL)
          if (is.null(labs) || is.null(input$coloc_level) ||
              !nzchar(input$coloc_level)) {
            status(list(type = "danger",
                        message = "Select a class variable and region level to colocalize against."))
            return()
          }
          mask <- as.character(labs) == input$coloc_level
          ref_desc <- sprintf("region '%s'", input$coloc_level)
        } else {
          if (is.null(input$coloc_mz) || is.na(input$coloc_mz)) {
            status(list(type = "danger",
                        message = "Enter a reference m/z to colocalize against."))
            return()
          }
          mz_ref <- input$coloc_mz
          ref_desc <- sprintf("m/z %s", signif(input$coloc_mz, 6))
        }
        top <- if (is.null(input$coloc_top) || is.na(input$coloc_top)) 3L else as.integer(input$coloc_top)
        coloc <- with_progress_stages("Finding co-localized features", {
          tryCatch(
            run_colocalization(data, mz = mz_ref, mask = mask, n = 50),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(coloc)) return()

        result(list(
          kind = "coloc",
          table = coloc$table,
          coloc = coloc,
          plot_fn = function() plot_colocalization(coloc, data, top = top)
        ))
        status(list(type = "success",
                    message = sprintf("Colocalization complete: %d features ranked against %s.",
                                      nrow(coloc$table), ref_desc)))

      } else {
        cl <- with_progress_stages("Running k-means clustering", {
          tryCatch(
            run_clustering(data, k = input$clust_k,
                           r = if (is.null(input$clust_r)) 1 else input$clust_r),
            error = function(e) {
              status(list(type = "danger", message = conditionMessage(e)))
              NULL
            }
          )
        })
        if (is.null(cl)) return()

        # Coords + per-pixel cluster assignment.
        tbl <- cbind(cl$coord, Cluster = as.integer(cl$cluster))
        result(list(
          kind = "clustering",
          method = "kmeans",
          table = tbl,
          clustering = cl,
          plot_fn = function() plot_clustering(cl)
        ))
        status(list(type = "success",
                    message = sprintf("k-means clustering complete: %d clusters over %d pixels (sizes: %s).",
                                      cl$k, length(cl$cluster),
                                      paste(as.integer(cl$sizes), collapse = ", "))))
      }

      # On success, record timing and memoize the result.
      if (identical(status()$type, "success")) {
        record_perf(input$analysis_type, proc.time()[["elapsed"]] - run_t0)
        if (cache_on && !is.null(sig)) {
          analysis_cache$set(sig, list(result = result(), status = status()))
          # Bump cache_version so the Performance-tab counter re-renders.
          if (!is.null(perf_state)) {
            perf_state$cache_version <-
              (if (is.null(perf_state$cache_version)) 0L else perf_state$cache_version) + 1L
          }
        }
      }
    })

    output$result_plot <- renderPlot({
      res <- result()
      validate(need(!is.null(res), "Run an analysis to see results."))
      res$plot_fn()
    })

    output$result_table <- DT::renderDT({
      res <- result()
      validate(need(!is.null(res), "Run an analysis to see results."))
      tbl <- res$table
      num_cols <- vapply(tbl, is.numeric, logical(1))
      DT::datatable(tbl, rownames = FALSE,
                    options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatSignif(columns = names(tbl)[num_cols], digits = 5)
    })

    # Plot-view choices per result kind.
    decomp_view_choices <- c("Scores scatter" = "scores",
                             "Loadings" = "loadings",
                             "Component contribution" = "variance",
                             "Spatial score map" = "image")

    da_view_choices <- c("Scores scatter" = "scores",
                         "Response image" = "response",
                         "Coefficients" = "coefficients")

    ssc_view_choices <- c("Predicted class map" = "class",
                          "Class probability" = "probability",
                          "Statistic (t-stat)" = "statistic")

    output$result_controls <- renderUI({
      res <- result()
      if (is.null(res)) return(NULL)
      if (res$kind %in% c("pca", "nmf")) {
        decomp <- res[[res$kind]]
        comps <- colnames(decomp$scores)
        view_id <- paste0(res$kind, "_view")
        x_id <- paste0(res$kind, "_x")
        y_id <- paste0(res$kind, "_y")
        return(div(class = "d-flex gap-2 align-items-center flex-wrap",
          selectInput(ns(view_id), NULL, choices = decomp_view_choices,
                      selected = if (!is.null(input[[view_id]])) input[[view_id]] else "scores",
                      width = "170px"),
          selectInput(ns(x_id), NULL, choices = comps,
                      selected = if (!is.null(input[[x_id]]) && input[[x_id]] %in% comps) input[[x_id]] else comps[1],
                      width = "90px"),
          selectInput(ns(y_id), NULL, choices = comps,
                      selected = if (!is.null(input[[y_id]]) && input[[y_id]] %in% comps) input[[y_id]] else comps[min(2, length(comps))],
                      width = "90px")
        ))
      }
      if (identical(res$kind, "da")) {
        comps <- names(res$da$scores)
        # Add the CV plot as a selectable view when CV was run.
        view_choices <- da_view_choices
        if (!is.null(res$cv)) {
          view_choices <- c(view_choices, "Cross-validation" = "cv")
        }
        view <- if (!is.null(input$da_view) && input$da_view %in% view_choices) input$da_view else "scores"
        # Component pickers only apply to the scores scatter (2+ comps).
        show_comps <- identical(view, "scores") && length(comps) >= 2
        return(div(class = "d-flex gap-2 align-items-center flex-wrap",
          selectInput(ns("da_view"), NULL, choices = view_choices,
                      selected = view, width = "170px"),
          if (show_comps) selectInput(ns("da_x"), NULL, choices = comps,
                      selected = if (!is.null(input$da_x) && input$da_x %in% comps) input$da_x else comps[1],
                      width = "90px"),
          if (show_comps) selectInput(ns("da_y"), NULL, choices = comps,
                      selected = if (!is.null(input$da_y) && input$da_y %in% comps) input$da_y else comps[min(2, length(comps))],
                      width = "90px")
        ))
      }
      if (identical(res$kind, "sscda")) {
        # Add the CV plot as a selectable view when CV was run.
        view_choices <- ssc_view_choices
        if (!is.null(res$cv)) {
          view_choices <- c(view_choices, "Cross-validation" = "cv")
        }
        view <- if (!is.null(input$ssc_view) && input$ssc_view %in% view_choices) input$ssc_view else "class"
        return(div(class = "d-flex gap-2 align-items-center flex-wrap",
          selectInput(ns("ssc_view"), NULL, choices = view_choices,
                      selected = view, width = "190px")
        ))
      }
      NULL
    })

    # Plot closure for the chosen decomposition view; shared by PCA/NMF.
    decomp_plot_fn <- function(kind, decomp, groups, view, comp_x, comp_y) {
      if (identical(kind, "pca")) {
        switch(view,
          loadings = function() plot_pca_loadings(decomp, comp_x),
          variance = function() plot_pca_variance(decomp),
          image = function() plot_pca_image(decomp, comp_x),
          function() plot_pca_scores(decomp, comp_x, comp_y, groups))
      } else {
        switch(view,
          loadings = function() plot_nmf_loadings(decomp, comp_x),
          variance = function() plot_nmf_variance(decomp),
          image = function() plot_nmf_image(decomp, comp_x),
          function() plot_nmf_scores(decomp, comp_x, comp_y, groups))
      }
    }

    observeEvent(list(input$pca_view, input$pca_x, input$pca_y), {
      res <- result()
      if (is.null(res) || res$kind != "pca") return()
      req(input$pca_x)
      view <- if (is.null(input$pca_view)) "scores" else input$pca_view
      res$plot_fn <- decomp_plot_fn("pca", res$pca, res$groups, view,
                                    input$pca_x, input$pca_y)
      result(res)
    }, ignoreInit = TRUE)

    observeEvent(list(input$nmf_view, input$nmf_x, input$nmf_y), {
      res <- result()
      if (is.null(res) || res$kind != "nmf") return()
      req(input$nmf_x)
      view <- if (is.null(input$nmf_view)) "scores" else input$nmf_view
      res$plot_fn <- decomp_plot_fn("nmf", res$nmf, res$groups, view,
                                    input$nmf_x, input$nmf_y)
      result(res)
    }, ignoreInit = TRUE)

    # DA plot closure; `cv` enables the "cv" view.
    da_plot_fn <- function(da, groups, view, comp_x, comp_y, cv = NULL) {
      switch(view,
        response = function() plot_da_response(da),
        coefficients = function() plot_da_coefficients(da),
        cv = if (!is.null(cv)) function() plot_cross_validation(cv)
             else function() plot_da_scores(da, comp_x = comp_x, comp_y = comp_y, groups = groups),
        function() plot_da_scores(da, comp_x = comp_x, comp_y = comp_y, groups = groups))
    }

    observeEvent(list(input$da_view, input$da_x, input$da_y), {
      res <- result()
      if (is.null(res) || res$kind != "da") return()
      view <- if (is.null(input$da_view)) "scores" else input$da_view
      comps <- names(res$da$scores)
      dx <- if (!is.null(input$da_x) && input$da_x %in% comps) input$da_x else comps[1]
      dy <- if (!is.null(input$da_y) && input$da_y %in% comps) input$da_y else comps[min(2, length(comps))]
      res$plot_fn <- da_plot_fn(res$da, res$groups, view, dx, dy, cv = res$cv)
      result(res)
    }, ignoreInit = TRUE)

    # SSC plot closure; `cv` enables the "cv" view.
    ssc_plot_fn <- function(ssc, view, cv = NULL) {
      switch(view,
        probability = function() plot_ssc_probability(ssc),
        statistic = function() plot_ssc_statistic(ssc),
        cv = if (!is.null(cv)) function() plot_cross_validation(cv)
             else function() plot_ssc_class(ssc),
        function() plot_ssc_class(ssc))
    }

    observeEvent(input$ssc_view, {
      res <- result()
      if (is.null(res) || res$kind != "sscda") return()
      view <- if (is.null(input$ssc_view)) "class" else input$ssc_view
      res$plot_fn <- ssc_plot_fn(res$ssc, view, cv = res$cv)
      result(res)
    }, ignoreInit = TRUE)

    # Save button shown for score/class/cluster results.
    output$save_scores_ui <- renderUI({
      res <- result()
      if (is.null(res) || !res$kind %in% c("pca", "nmf", "da", "sscda", "clustering")) {
        return(div(class = "text-muted small",
                   "Run PCA, NMF, PLS-DA/OPLS-DA, SSC classification, or clustering to enable saving results."))
      }
      label <- switch(res$kind,
        clustering = "Save clusters as dataset",
        sscda = "Save predicted class as dataset",
        "Save scores as dataset")
      actionButton(ns("save_scores"), label,
                   icon = icon("floppy-disk"), class = "btn-outline-success w-100")
    })

    observeEvent(input$save_scores, {
      res <- result()
      req(!is.null(res), res$kind %in% c("pca", "nmf", "da", "sscda", "clustering"))
      default_name <- switch(res$kind,
        da = paste0(res$method, "da_scores"),
        sscda = "ssc_class",
        clustering = paste0(res$method, "_clusters"),
        nmf = "nmf_scores",
        "pca_scores")
      title <- switch(res$kind,
        da = sprintf("Save %s-DA scores as a new dataset", toupper(res$method)),
        sscda = "Save SSC predicted class as a new dataset",
        clustering = sprintf("Save %s clusters as a new dataset", "k-means"),
        nmf = "Save NMF scores as a new dataset",
        "Save PCA scores as a new dataset")
      save_dataset_modal(ns, title, default_name, "save_scores_confirm")
    })

    observeEvent(input$save_scores_confirm, {
      removeModal()
      res <- result()
      data <- mse_data()
      req(!is.null(res), res$kind %in% c("pca", "nmf", "da", "sscda", "clustering"), !is.null(data))
      ds <- datasets()
      new_id <- next_dataset_id(names(ds))
      record <- tryCatch(
        switch(res$kind,
          da = make_da_record(data, res$da, new_id,
                              filename = input$save_filename, ext = input$save_ext),
          sscda = make_ssc_record(data, res$ssc, new_id,
                              filename = input$save_filename, ext = input$save_ext),
          clustering = make_clustering_record(data, res$clustering, new_id,
                              filename = input$save_filename, ext = input$save_ext),
          nmf = make_nmf_record(data, res$nmf$scores, new_id,
                              filename = input$save_filename, ext = input$save_ext),
          make_pca_record(data, res$pca$scores, new_id,
                          filename = input$save_filename, ext = input$save_ext)
        ),
        error = function(e) {
          status(list(type = "danger", message = conditionMessage(e)))
          NULL
        }
      )
      if (is.null(record)) return()
      ds[[new_id]] <- record
      datasets(ds)
      active_id(new_id)
      status(list(type = "success",
                  message = sprintf("Scores saved as '%s'.", record$name)))
    })

    output$export_csv <- downloadHandler(
      filename = function() {
        res <- result()
        kind <- if (is.null(res)) "results" else res$kind
        paste0("analysis_", kind, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
      },
      content = function(file) {
        res <- result()
        if (is.null(res)) {
          utils::write.csv(data.frame(), file, row.names = FALSE)
        } else {
          utils::write.csv(res$table, file, row.names = FALSE)
        }
      }
    )

    output$export_plot <- export_plot_handler(
      result, function() input$plot_format,
      basename = function() {
        res <- result()
        paste0("analysis_", if (is.null(res)) "plot" else res$kind)
      },
      empty_msg = "No analysis result to export.")

    output$status <- render_status(status)
  })
}
