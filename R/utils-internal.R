# =============================================================================
# Internal helpers for HaploCatcher
# -----------------------------------------------------------------------------
# Every function in this file is hidden (the leading "." keeps it out of the
# package namespace / autocomplete) and undocumented for end users. These
# helpers centralize the logic that used to be copy-pasted across locus_cv(),
# locus_train(), locus_perm_cv() and auto_locus(), so the exported functions
# stay short, fast, and consistent. They also do the heavy input checking in
# one place, which makes error messages uniform and the public code lean.
# =============================================================================

# ---- argument validators ---------------------------------------------------

# Stop unless `x` is a single TRUE/FALSE value.
.hc_assert_logical <- function(x, name) {
  # reject anything that is not a length-one logical
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop(sprintf("Argument '%s' must be a single logical (TRUE or FALSE)!", name), call. = FALSE)
  # return invisibly so the helper can be used purely for its side effect
  invisible(TRUE)
}

# Stop unless `x` is a single proportion strictly inside (0, 1).
.hc_assert_prop <- function(x, name) {
  # reject non-numeric, non-scalar, or out-of-range values
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1)
    stop(sprintf("Argument '%s' must be a single number strictly between 0 and 1!", name), call. = FALSE)
  # invisible TRUE keeps call sites clean
  invisible(TRUE)
}

# Stop unless `x` is a single positive whole number (count-like argument).
.hc_assert_count <- function(x, name) {
  # reject non-numeric, non-scalar, or non-positive values
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1)
    stop(sprintf("Argument '%s' must be a single positive number!", name), call. = FALSE)
  # invisible TRUE for side-effect use
  invisible(TRUE)
}

# Stop unless `x` is NULL or a non-empty character vector (case-label argument).
.hc_assert_labels <- function(x, name) {
  # NULL means "use the default gene/het_gene/non_gene naming convention"
  if (is.null(x)) return(invisible(TRUE))
  # otherwise it must be a non-empty character vector of call labels
  if (!is.character(x) || length(x) == 0L)
    stop(sprintf("Argument '%s' must be NULL or a character vector of call labels!", name), call. = FALSE)
  # ok
  invisible(TRUE)
}

# ---- call-case identification ----------------------------------------------
# These two helpers decide which calls are heterozygous vs negative. By default
# they follow the package convention (positive = "<gene>", het = "het_<gene>",
# negative = "non_<gene>"); supplying explicit labels overrides the convention
# so users are not forced into that naming scheme.

# TRUE for calls that should be treated as heterozygous.
.hc_is_het <- function(calls, het_label = NULL) {
  # exact-match the supplied labels, or fall back to the "het_" prefix
  if (is.null(het_label)) grepl("het_", calls, ignore.case = TRUE)
  else as.character(calls) %in% het_label
}

# TRUE for calls that should be treated as the negative / wild-type case.
.hc_is_neg <- function(calls, neg_label = NULL) {
  # exact-match the supplied labels, or fall back to the "non_" prefix
  if (is.null(neg_label)) grepl("non_", calls, ignore.case = TRUE)
  else as.character(calls) %in% neg_label
}

# ---- core data-object validation -------------------------------------------

# Validate the three primary data objects plus the gene/chromosome identifiers.
# Run once at the top of every public entry point so downstream helpers can
# assume well-formed inputs.
.hc_check_inputs <- function(geno_mat, gene_file, marker_info, gene_name, chromosome) {

  # --- genotype matrix ---
  # the genotype matrix must be a numeric matrix
  if (!is.matrix(geno_mat) || !is.numeric(geno_mat))
    stop("'geno_mat' must be a numeric marker matrix!", call. = FALSE)
  # it must carry both individual (row) and marker (column) names
  if (is.null(rownames(geno_mat)) || is.null(colnames(geno_mat)))
    stop("'geno_mat' must have row names (individuals) and column names (markers)!", call. = FALSE)
  # individuals (rows) must be uniquely named
  if (anyDuplicated(rownames(geno_mat)))
    stop("There are duplicated individuals (row names) in 'geno_mat'!", call. = FALSE)
  # missing genotype values are not allowed
  if (anyNA(geno_mat))
    stop("'geno_mat' contains missing values; impute before use!", call. = FALSE)

  # --- gene file ---
  # the required columns must all be present
  need_gene <- c("Gene", "FullSampleName", "Call")
  if (!is.data.frame(gene_file) || !all(need_gene %in% names(gene_file)))
    stop("'gene_file' must contain the columns 'Gene', 'FullSampleName', and 'Call'!", call. = FALSE)
  # the requested gene must actually appear in the gene file
  if (!gene_name %in% gene_file$Gene)
    stop("'gene_name' was not found in the 'gene_file' 'Gene' column!", call. = FALSE)

  # --- marker information ---
  # the required columns must all be present
  need_marker <- c("Marker", "Chromosome", "BP_Position")
  if (!is.data.frame(marker_info) || !all(need_marker %in% names(marker_info)))
    stop("'marker_info' must contain the columns 'Marker', 'Chromosome', and 'BP_Position'!", call. = FALSE)
  # the requested chromosome must be present in the marker table
  if (!chromosome %in% marker_info$Chromosome)
    stop("'chromosome' was not found in the 'marker_info' 'Chromosome' column!", call. = FALSE)

  # all good
  invisible(TRUE)
}

# ---- classification builder ------------------------------------------------

# Subset the gene file to one gene, enforce a biallelic structure, and
# optionally drop heterozygous calls. `het_label` overrides the default
# "het_" naming convention. Returns the cleaned classification frame.
.hc_classification <- function(gene_file, gene_name, include_hets, verbose, het_label = NULL) {

  # keep only the records belonging to the requested gene
  classification <- gene_file[gene_file$Gene == gene_name, , drop = FALSE]

  # an individual may only appear once per gene
  if (anyDuplicated(classification$FullSampleName))
    stop("There are duplicated individuals in the 'gene_file' for this gene!", call. = FALSE)

  # flag which rows are heterozygous under the chosen convention/labels
  het_rows <- .hc_is_het(classification$Call, het_label)
  # categories present once hets are set aside (the homozygous classes)
  n_homo <- length(unique(classification$Call[!het_rows]))
  # a single homozygous category cannot train a classifier (no negative cases)
  if (n_homo < 2L)
    stop("Fewer than two homozygous call categories are present; need positive and negative cases!", call. = FALSE)
  # only biallelic loci are supported (at most two homozygous classes)
  if (n_homo > 2L)
    stop("More than two homozygous call categories were found; only biallelic loci are supported!", call. = FALSE)
  # asking to keep hets but supplying none signals a formatting/label problem
  if (include_hets && !any(het_rows))
    stop("'include_hets = TRUE' but no heterozygous calls were found; check your data or 'het_label'!", call. = FALSE)

  # drop heterozygous calls unless the user keeps them
  if (!include_hets) {
    # tell the user we are removing hets when verbose
    if (verbose) message("Note: Removing heterozygous calls from the data.")
    # keep only the non-heterozygous rows
    classification <- classification[!het_rows, , drop = FALSE]
  } else if (verbose) {
    # otherwise note that hets are retained
    message("Note: Heterozygous calls retained in the data.")
  }

  # hand back the cleaned classification frame
  classification
}

# ---- chromosome genotype subsetting ----------------------------------------

# Pull the chromosome-specific marker block of the genotype matrix for the
# individuals present in `classification`, and verify marker agreement.
.hc_geno_chr <- function(geno_mat, marker_info, chromosome, classification) {

  # markers that sit on the requested chromosome
  chr_markers <- marker_info[marker_info$Chromosome == chromosome, , drop = FALSE]
  # rows = individuals in the classification; columns = chromosome markers
  geno_chr <- geno_mat[rownames(geno_mat) %in% classification$FullSampleName,
                       colnames(geno_mat) %in% chr_markers$Marker, drop = FALSE]
  # every chromosome marker must exist in the genotype matrix
  if (ncol(geno_chr) != nrow(chr_markers))
    stop("Some markers in 'marker_info' are absent from 'geno_mat'; the files disagree!", call. = FALSE)
  # at least some individuals must overlap between files
  if (nrow(geno_chr) == 0L)
    stop("None of the gene-file individuals were found in 'geno_mat'!", call. = FALSE)

  # return both the genotype block and the marker metadata
  list(geno = geno_chr, markers = chr_markers)
}

# ---- fast correlation-based marker selection -------------------------------

# Select the top `ncor_markers` markers most correlated with the call.
# Speed note: the original code built a full m x m correlation matrix and then
# discarded all but one column. Here we correlate the encoded call vector
# against the marker block in a single pass -> O(m) instead of O(m^2).
.hc_top_markers <- function(geno_sub, call_vec, chr_markers, ncor_markers,
                            graph = FALSE, chromosome = NULL, gene_name = NULL) {

  # never request more markers than are available
  if (ncor_markers > ncol(geno_sub))
    stop("'ncor_markers' exceeds the number of markers on this chromosome; lower it!", call. = FALSE)

  # numerically encode the categorical call so it can be correlated
  call_num <- match(call_vec, unique(call_vec))
  # correlate every marker column against the encoded call in one shot
  r <- suppressWarnings(stats::cor(call_num, geno_sub))
  # collapse the 1 x m result to a rounded absolute correlation vector
  r <- abs(round(as.numeric(r), 10))

  # assemble a tidy table of marker, |r|, and genomic position
  corr <- data.frame(marker = colnames(geno_sub),
                     `|r|`  = r,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
  # attach basepair position by matching marker names (robust to ordering)
  corr$BP_Position  <- chr_markers$BP_Position[match(corr$marker, chr_markers$Marker)]
  # derive megabasepair position for plotting
  corr$MBP_Position <- round(as.numeric(corr$BP_Position) / 1e6, 2)

  # optional diagnostic scatter of correlation strength across the chromosome
  if (graph) {
    # plot |r| against physical position
    graphics::plot(x = corr$MBP_Position, y = corr$`|r|`,
                   main = paste("Marker correlations on", chromosome, "for", gene_name),
                   xlab = "Megabasepair (Mbp) Position",
                   ylab = "Absolute Correlation (|r|)")
  }

  # order markers from strongest to weakest correlation (NAs sink to the end)
  corr <- corr[order(corr$`|r|`, decreasing = TRUE), , drop = FALSE]
  # keep only the requested number of top markers
  corr <- corr[seq_len(ncor_markers), , drop = FALSE]

  # draw the selection threshold once markers are chosen
  if (graph) {
    # horizontal line at the weakest retained correlation
    graphics::abline(h = min(corr$`|r|`), col = "red", lty = 2)
    # legend describing the threshold
    graphics::legend("bottomleft",
                     legend = paste("Top", ncor_markers, "correlated markers"),
                     col = "red", lty = 2)
  }

  # return the selected-marker table
  corr
}

# ---- model-ready frame assembly --------------------------------------------

# Build a [FullSampleName, Call, <markers>] data frame aligned by individual,
# with identifiers/call as factors and markers numeric.
.hc_assemble <- function(part_cls, geno_chr, marker_names, call_levels) {

  # order the classification rows by individual name
  part <- part_cls[order(part_cls$FullSampleName), , drop = FALSE]
  # pull the selected markers for these individuals, ordered to match
  g <- geno_chr[rownames(geno_chr) %in% part$FullSampleName, marker_names, drop = FALSE]
  g <- g[order(rownames(g)), , drop = FALSE]

  # combine identifiers, call, and marker genotypes into one frame
  out <- data.frame(FullSampleName = part$FullSampleName,
                    Call           = part$Call,
                    g,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  # identifiers as a factor
  out$FullSampleName <- as.factor(out$FullSampleName)
  # call as a factor with a consistent, caller-supplied level set
  out$Call <- factor(out$Call, levels = call_levels)
  # markers coerced to numeric (defensive; they should already be numeric)
  out[, -(1:2)] <- lapply(out[, -(1:2), drop = FALSE], as.numeric)
  # hand back the assembled frame
  out
}

# ---- model fitting ---------------------------------------------------------

# Build a lean 5-fold CV control object. The original used
# method = "repeatedcv" with repeats = 1000 while also passing an explicit
# `index`, which makes caret ignore the repeats entirely -- wasted bookkeeping.
# A plain "cv" with the same folds is identical in effect and clearer.
.hc_train_control <- function(y) {
  # create five outcome-stratified folds
  folds <- caret::createFolds(y, k = 5)
  # configure simple k-fold CV over exactly those folds
  caret::trainControl(method = "cv", number = 5, index = folds,
                      classProbs = FALSE, savePredictions = "final",
                      verboseIter = FALSE)
}

# KNN tuning grid: a single k when neighbors are scarce, else odd k's.
.hc_knn_grid <- function(n_neighbors) {
  # with two or fewer neighbors there is nothing to tune
  if (n_neighbors <= 2) expand.grid(k = 1)
  # otherwise sweep odd k values up to the requested maximum
  else expand.grid(k = seq(1, n_neighbors, by = 2))
}

# RF tuning grid: a single mtry when predictors are scarce, else steps of 5.
.hc_rf_grid <- function(n_pred) {
  # too few predictors -> only mtry = 1 is meaningful
  if (n_pred < 5) expand.grid(mtry = 1)
  # otherwise test mtry = 1, 5, 10, ... up to the predictor count
  else expand.grid(mtry = c(1, seq(0, n_pred, by = 5)[-1]))
}

# Fit one caret model, falling back to leave-one-out CV if k-fold fails.
.hc_fit_one <- function(data_xy, method, grid, control) {
  # attempt the primary k-fold fit
  fit <- try(caret::train(Call ~ ., data = data_xy, method = method,
                          tuneGrid = grid, trControl = control), silent = TRUE)
  # on failure, warn and retry with leave-one-out CV
  if (inherits(fit, "try-error")) {
    warning(sprintf("k-fold CV failed for '%s'; retrying with leave-one-out CV.", method))
    fit <- try(caret::train(Call ~ ., data = data_xy, method = method,
                            tuneGrid = grid,
                            trControl = caret::trainControl(method = "LOOCV")), silent = TRUE)
    # if it still fails the data structure is the problem
    if (inherits(fit, "try-error"))
      stop(sprintf("Hyper-parameter tuning of '%s' failed; check the data structure!", method), call. = FALSE)
  }
  # return the fitted caret model
  fit
}

# Fit the requested model set on an assembled training frame.
# `models_request` is one of "all", "knn", "rf".
.hc_fit_models <- function(training, n_neighbors, models_request, verbose = FALSE) {
  # number of marker predictors (drop FullSampleName + Call)
  n_pred <- ncol(training) - 2
  # shared CV control built from the training outcome
  control <- .hc_train_control(training$Call)
  # predictors + outcome only (drop the FullSampleName identifier column)
  data_xy <- training[, -1, drop = FALSE]
  # announce model fitting when verbose
  if (verbose) message("Note: Running models...")

  # container for the fitted models
  out <- list()
  # fit KNN when requested
  if (models_request %in% c("all", "knn"))
    out$knn <- .hc_fit_one(data_xy, "knn", .hc_knn_grid(n_neighbors), control)
  # fit RF when requested
  if (models_request %in% c("all", "rf"))
    out$rf <- .hc_fit_one(data_xy, "rf", .hc_rf_grid(n_pred), control)

  # report completion when verbose
  if (verbose) message("Note: Done fitting models.")
  # return the named list of fits
  out
}

# ---- evaluation (predictions + confusion matrices) -------------------------

# Wrap a prediction vector into the tidy per-model prediction frame used
# throughout the package.
.hc_pred_frame <- function(ids, model_label, gene_name, observed, predicted, call_levels) {
  data.frame(FullSampleName = ids,
             Model          = model_label,
             Gene           = gene_name,
             Observed_Call  = factor(observed,  levels = call_levels),
             Predicted_Call = factor(predicted, levels = call_levels),
             stringsAsFactors = FALSE)
}

# Predict the test partition with each fitted model and build confusion
# matrices. KNN can fail when too many neighbor ties occur; that case yields
# NA predictions/confusion for KNN while RF proceeds normally.
.hc_eval_models <- function(fits, testing, gene_name, call_levels, verbose = FALSE) {

  # predictor + call columns only (predict ignores the extra Call column)
  newx <- testing[, -1, drop = FALSE]
  # observed calls and individual ids for the test set
  observed <- testing$Call
  ids      <- as.character(testing$FullSampleName)

  # --- random forest (treated as robust) ---
  rf_hat  <- stats::predict(fits$rf, newx)
  pred_rf <- .hc_pred_frame(ids, "Random Forest", gene_name, observed, rf_hat, call_levels)
  confu_rf <- caret::confusionMatrix(pred_rf$Observed_Call, pred_rf$Predicted_Call)

  # --- k-nearest neighbors (guard against tie failures) ---
  knn_hat <- try(stats::predict(fits$knn, newx), silent = TRUE)
  if (inherits(knn_hat, "try-error")) {
    # too many ties -> emit a warning and record NA results for KNN
    warning("Too many ties in the K-Nearest Neighbors predictions; returning NA for KNN.")
    pred_knn  <- .hc_pred_frame(ids, "K-Nearest Neighbors", gene_name, observed, NA, call_levels)
    confu_knn <- NA
  } else {
    # normal path: build the KNN prediction frame and confusion matrix
    pred_knn  <- .hc_pred_frame(ids, "K-Nearest Neighbors", gene_name, observed, knn_hat, call_levels)
    confu_knn <- caret::confusionMatrix(pred_knn$Observed_Call, pred_knn$Predicted_Call)
  }

  # optionally print the confusion tables
  if (verbose) {
    if (!identical(confu_knn, NA))
      print(knitr::kable(confu_knn$table, caption = "Confusion matrix: K-Nearest Neighbors"))
    print(knitr::kable(confu_rf$table, caption = "Confusion matrix: Random Forest"))
  }

  # return predictions and confusion matrices keyed by model
  list(preds = list(knn = pred_knn, rf = pred_rf),
       confu = list(knn = confu_knn, rf = confu_rf))
}

# ---- call-frequency reporting ----------------------------------------------

# Pretty-print the relative frequency of each call (verbose mode only).
.hc_freq_table <- function(call, caption) {
  # relative frequency of each call level
  tab <- data.frame(table(call) / length(call))
  # tidy column names
  names(tab) <- c("Call", "Frequency")
  # render as a knitr table
  print(knitr::kable(tab, caption = caption, digits = 2))
}

# ---- permutation summarization ---------------------------------------------

# Pull the overall accuracy/kappa row for one model from one permutation.
.hc_overall_row <- function(cm, model_label, perm) {
  # NA confusion (KNN ties) -> NA metrics
  if (identical(cm, NA))
    return(data.frame(Permutation = perm, Model = model_label,
                      Accuracy = NA_real_, Kappa = NA_real_))
  # otherwise extract accuracy and kappa from the confusion object
  data.frame(Permutation = perm, Model = model_label,
             Accuracy = unname(cm$overall["Accuracy"]),
             Kappa    = unname(cm$overall["Kappa"]))
}

# Pull the by-class sensitivity/specificity/precision/recall/balanced-accuracy
# rows for one model from one permutation. Handles both the multiclass (het)
# matrix form and the two-class vector form of confusionMatrix$byClass.
.hc_byclass_rows <- function(cm, model_label, perm) {
  # the five by-class metrics we keep
  keep <- c("Sensitivity", "Specificity", "Precision", "Recall", "Balanced Accuracy")
  # NA confusion -> a single NA row
  if (identical(cm, NA))
    return(data.frame(Permutation = perm, Model = model_label, Class = NA_character_,
                      Sensitivity = NA_real_, Specificity = NA_real_, Precision = NA_real_,
                      Recall = NA_real_, Balanced_Accuracy = NA_real_))
  bc <- cm$byClass
  if (is.matrix(bc)) {
    # multiclass: one row per class, class label from the row names
    df <- as.data.frame(bc[, keep, drop = FALSE], stringsAsFactors = FALSE)
    cls <- gsub("Class: ", "", rownames(bc))
  } else {
    # two-class: a named vector; the positive class is the single class
    df  <- as.data.frame(t(bc[keep]), stringsAsFactors = FALSE)
    cls <- if (!is.null(cm$positive)) cm$positive else "Overall"
  }
  # force consistent underscore column names (as.data.frame munges spaces to
  # dots, which would clash with the NA-row branch above and break rbind)
  names(df) <- c("Sensitivity", "Specificity", "Precision", "Recall", "Balanced_Accuracy")
  # prepend identifiers
  data.frame(Permutation = perm, Model = model_label, Class = cls, df,
             row.names = NULL, check.names = FALSE)
}

# Aggregate numeric columns of `x` by the grouping columns, returning
# Mean_/Min_/Max_/SD_ statistics with NAs treated as zero (matches original).
.hc_summarize_by <- function(x, group_cols) {
  # metric columns are the numeric columns that are neither grouping keys
  # nor the permutation index (Class/Model are character and excluded here)
  numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
  metric_cols <- setdiff(numeric_cols, c(group_cols, "Permutation"))
  # replace NA metrics with zero so summaries are defined even with KNN ties
  x[metric_cols] <- lapply(x[metric_cols], function(v) ifelse(is.na(v), 0, v))
  # build the grouping formula once
  fml <- stats::as.formula(paste(". ~", paste(group_cols, collapse = " + ")))
  # data restricted to grouping + metric columns for aggregate()
  sub <- x[, c(group_cols, metric_cols), drop = FALSE]
  # compute the four summary statistics
  agg <- function(f) stats::aggregate(fml, data = sub, FUN = f)
  m <- agg(mean); mn <- agg(min); mx <- agg(max); sdv <- agg(stats::sd)
  # start from the grouping columns of the mean table
  out <- m[, group_cols, drop = FALSE]
  # interleave Mean/Min/Max/SD for each metric
  for (mc in metric_cols)
    out <- cbind(out,
                 stats::setNames(data.frame(m[[mc]], mn[[mc]], mx[[mc]], sdv[[mc]]),
                                 paste0(c("Mean_", "Min_", "Max_", "SD_"), mc)))
  # return the assembled summary
  out
}

# Turn the raw list of locus_cv() results into the five-element summary object
# returned by locus_perm_cv(). Replaces two near-identical het/no-het blocks.
.hc_summarize_perm <- function(results, gene_name) {

  # ---- overall parameters (one row per model per permutation) ----
  overall <- do.call(rbind, lapply(seq_along(results), function(i) {
    cm <- results[[i]]$confusion_matrices
    rbind(.hc_overall_row(cm$knn, "K-Nearest Neighbors", i),
          .hc_overall_row(cm$rf,  "Random Forest",       i))
  }))

  # ---- by-class parameters (one row per class per model per permutation) ----
  byclass <- do.call(rbind, lapply(seq_along(results), function(i) {
    cm <- results[[i]]$confusion_matrices
    rbind(.hc_byclass_rows(cm$knn, "K-Nearest Neighbors", i),
          .hc_byclass_rows(cm$rf,  "Random Forest",       i))
  }))

  # ---- summaries across permutations ----
  # overall summarized by model
  overall_summary <- .hc_summarize_by(overall, "Model")
  # by-class summarized by model (and class when more than one class exists)
  group_cols <- if (length(unique(stats::na.omit(byclass$Class))) > 1) c("Model", "Class") else "Model"
  byclass_summary <- .hc_summarize_by(byclass, group_cols)

  # warn if KNN collapsed to all-NA (only one model effectively survived)
  if (all(is.na(overall$Accuracy[overall$Model == "K-Nearest Neighbors"])))
    warning(sprintf("Too many ties in the K-Nearest Neighbors predictions of %s; check results for NAs!", gene_name))

  # assemble the five-element return object
  list(Overall_Parameters   = overall,
       By_Class_Parameters  = byclass,
       Overall_Summary      = overall_summary,
       By_Class_Summary     = byclass_summary,
       Raw_Permutation_Info = stats::setNames(results, paste0("Permutation_", seq_along(results))))
}

# ---- consensus voting ------------------------------------------------------

# Tabulate per-model votes across runs into a vote matrix and a majority-rule
# consensus (ties -> "No_Call"). `pred_list` is a list of locus_pred() frames,
# each with FullSampleName in column 1 and a single prediction in column 2.
.hc_tabulate_votes <- function(pred_list) {

  # individuals are taken from the first run's prediction frame
  ids <- as.character(pred_list[[1]]$FullSampleName)
  # the full set of possible call levels (union across runs, factor levels)
  lvls <- levels(factor(unlist(lapply(pred_list, function(p) as.character(p[[2]])))))

  # wide vote matrix: one column per vote
  votes <- data.frame(FullSampleName = ids, stringsAsFactors = FALSE)
  # fill each vote column from the corresponding run
  for (i in seq_along(pred_list))
    votes[[paste0("Vote_", i)]] <- as.character(pred_list[[i]][[2]])

  # per-individual class counts across all votes
  counts <- t(apply(votes[, -1, drop = FALSE], 1, function(row)
    table(factor(row, levels = lvls))))
  # ensure a matrix even when there is a single individual
  counts <- matrix(counts, nrow = nrow(votes), dimnames = list(NULL, lvls))

  # majority-rule consensus call; ties resolve to "No_Call"
  consensus_call <- apply(counts, 1, function(cnt) {
    top <- names(cnt)[cnt == max(cnt)]
    if (length(top) > 1) "No_Call" else top
  })

  # assemble the consensus frame: id, consensus call, and per-class counts
  consensus <- data.frame(FullSampleName = ids,
                          Consensus_Call = consensus_call,
                          counts,
                          check.names = FALSE,
                          stringsAsFactors = FALSE)

  # return both the raw vote matrix and the consensus summary
  list(votes = votes, consensus = consensus)
}

# ---- best-model selection --------------------------------------------------

# Choose "knn" or "rf" from a locus_perm_cv() result based on the selection
# metric, breaking ties at random.
.hc_select_model <- function(perm_results, metric) {
  # map the user metric onto the summary column name
  col <- if (metric == "accuracy") "Mean_Accuracy" else "Mean_Kappa"
  # the overall summary table holds one row per model
  s <- perm_results$Overall_Summary
  # models achieving the maximum mean metric
  best <- s$Model[which(s[[col]] == max(s[[col]], na.rm = TRUE))]
  # break ties (or empty) by random choice between the two model types
  if (length(best) != 1L) {
    message("Note: model tie on ", metric, "; selecting at random.")
    return(sample(c("knn", "rf"), 1))
  }
  # translate the human-readable model name to its short code
  if (best == "Random Forest") "rf" else "knn"
}

# ---- robust forward prediction ---------------------------------------------

# Human-readable name for a model short code.
.hc_model_name <- function(m) if (m == "knn") "k-nearest neighbors" else "random forest"

# Train the chosen model and forward-predict, falling back to the other model
# if the chosen one fails (e.g. KNN throws "too many ties in knn" at predict
# time). Returns list(fit, pred, model_used). Errors only if BOTH models fail.
.hc_forward_predict <- function(train_args, geno_mat, genotypes_to_predict, primary, verbose = TRUE) {

  # the alternate model to fall back to when the primary one fails
  alt <- if (primary == "knn") "rf" else "knn"

  # train with one model and forward-predict; the caller wraps this in try()
  attempt <- function(model) {
    # reuse the training arguments but force the requested model
    a <- train_args
    a$models_request <- model
    # train then predict the requested genotypes
    fit <- do.call(locus_train, a)
    pred <- locus_pred(fit, geno_mat, genotypes_to_predict)
    # bundle the fit, predictions, and which model actually ran
    list(fit = fit, pred = pred, model_used = model)
  }

  # try the CV-selected (primary) model first
  res <- try(attempt(primary), silent = TRUE)
  # success -> return immediately
  if (!inherits(res, "try-error")) return(res)

  # primary failed: warn (surfacing the underlying error) and try the alternate
  warning(sprintf("Forward prediction with the %s model failed: %s. Falling back to the %s model.",
                  .hc_model_name(primary),
                  conditionMessage(attr(res, "condition")),
                  .hc_model_name(alt)), call. = FALSE)
  res <- try(attempt(alt), silent = TRUE)
  # the fallback worked -> note it and return
  if (!inherits(res, "try-error")) {
    if (verbose) message("Forward prediction succeeded with the ", .hc_model_name(alt), " model.")
    return(res)
  }

  # both models failed: stop with the underlying error message
  stop(sprintf("Forward prediction failed for both the k-nearest neighbors and random forest models: %s. Check your data and parameters.",
               conditionMessage(attr(res, "condition"))), call. = FALSE)
}
