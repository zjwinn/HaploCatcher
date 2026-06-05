#' Auto Locus: An Automated Pipeline for Locus Prediction
#'
#' Weaves the HaploCatcher functions into a single pipeline (Figure 1b of the
#' package paper): permutation cross-validation, best-model selection by kappa
#' or accuracy, then forward prediction either with one seeded model or by
#' majority-rule voting over many random models.
#'
#' @param geno_mat An imputed, number-coded genotypic matrix with n rows of individuals and m columns of markers. Row names are genotype IDs; column names are marker IDs. Missing data are not allowed. Numeric coding may vary as long as it is consistent across markers.
#' @param gene_file A data frame with at least the columns 'Gene', 'FullSampleName', and 'Call'. 'Gene' is the gene each observation belongs to, 'FullSampleName' matches a column name in the genotypic matrix, and 'Call' is the marker call for that genotype.
#' @param gene_name A character string matching a value in the 'Gene' column of `gene_file`.
#' @param marker_info A data frame with the columns 'Marker', 'Chromosome', and 'BP_Position'. Every marker in the genotypic matrix must be listed. If positions are unavailable a numeric dummy (1..m) may be used.
#' @param chromosome A character string matching a value in the 'Chromosome' column of `marker_info`.
#' @param training_genotypes Character vector of FullSampleNames used for cross-validation and to train the prediction model.
#' @param testing_genotypes Character vector of FullSampleNames to predict.
#' @param ncor_markers Number of top correlated markers to retain for training. Default 50.
#' @param n_neighbors Number of neighbors to consider in KNN. Default 50.
#' @param cv_percent_testing Proportion reserved for validation during CV, strictly between 0 and 1. Default 0.20.
#' @param cv_percent_training Proportion used for training during CV, strictly between 0 and 1. Default 0.80.
#' @param n_perms Number of cross-validation permutations. Default 30.
#' @param model_selection_parameter Metric for selecting the best model: "kappa" or "accuracy". Default "kappa".
#' @param n_votes Number of models to train and predict with when voting. Default 30.
#' @param set_seed Numeric seed for a single reproducible prediction (required when `predict_by_vote = FALSE`). Default NULL.
#' @param predict_by_vote Logical; predict by majority rule over many random models. Default FALSE.
#' @param include_hets Logical; keep heterozygous calls. Default FALSE.
#' @param include_models Logical; keep the trained models in the CV result (large). Default FALSE.
#' @param verbose Logical; print progress and plots. Default TRUE.
#' @param parallel Logical; run CV and voting in parallel. Default FALSE.
#' @param n_cores Number of cores for parallel processing. If NULL and `parallel = TRUE`, uses all available cores minus one.
#' @param plot_cv_results Logical; draw the cross-validation summary plot. Default TRUE.
#' @param het_label Optional character vector of `Call` values to treat as heterozygous. When NULL (default), the package convention (calls containing "het_") is used.
#' @param neg_label Optional character vector of `Call` values to treat as the negative/wild-type case (used only for plot facet labels). When NULL (default), the "non_" prefix is used.
#'
#' @return A list. When `predict_by_vote = FALSE`: `method`,
#'   `cross_validation_results`, `prediction_model`, and `predictions`. When
#'   `predict_by_vote = TRUE`: `method`, `cross_validation_results`,
#'   `predictions` (per-vote calls), and `consensus_predictions` (majority rule).
#'
#' @export
#'
#' @examples
#'
#' #refer to vignette for an in depth look at the auto_locus function
#' vignette("An_Intro_to_HaploCatcher", package = "HaploCatcher")
#'
#' @importFrom randomForest randomForest
#' @importFrom ggplot2 ggplot
#' @importFrom caret train knn3 createFolds trainControl confusionMatrix
#' @importFrom foreach %dopar% foreach

auto_locus <- function(geno_mat, gene_file, gene_name, marker_info, chromosome,
                       training_genotypes, testing_genotypes,
                       ncor_markers = 50, n_neighbors = 50,
                       cv_percent_testing = 0.2, cv_percent_training = 0.8,
                       n_perms = 30, model_selection_parameter = "kappa",
                       n_votes = 30, set_seed = NULL, predict_by_vote = FALSE,
                       include_hets = FALSE, include_models = FALSE,
                       verbose = TRUE, parallel = FALSE, n_cores = NULL,
                       plot_cv_results = TRUE, het_label = NULL, neg_label = NULL) {

  # validate any custom case labels
  .hc_assert_labels(het_label, "het_label")
  .hc_assert_labels(neg_label, "neg_label")
  # a non-voting run needs a seed to be reproducible
  if (!predict_by_vote && is.null(set_seed))
    stop("Set 'set_seed' for a reproducible single-model run, or set 'predict_by_vote = TRUE'!", call. = FALSE)
  # normalize the selection metric, defaulting to kappa on bad input
  if (!model_selection_parameter %in% c("kappa", "accuracy")) {
    warning("'model_selection_parameter' must be 'kappa' or 'accuracy'; defaulting to 'kappa'.")
    model_selection_parameter <- "kappa"
  }

  # restrict the training data to the supplied training genotypes
  train_geno_mat  <- geno_mat[rownames(geno_mat) %in% training_genotypes, , drop = FALSE]
  train_gene_file <- gene_file[gene_file$FullSampleName %in% training_genotypes, , drop = FALSE]

  # ---- Step I: permutation cross-validation on the training data ----
  if (verbose) message("Cross validating ", gene_name, "...")
  fit_cv <- locus_perm_cv(n_perms = n_perms, geno_mat = train_geno_mat,
                          gene_file = train_gene_file, gene_name = gene_name,
                          marker_info = marker_info, chromosome = chromosome,
                          ncor_markers = ncor_markers, n_neighbors = n_neighbors,
                          percent_testing = cv_percent_testing,
                          percent_training = cv_percent_training,
                          include_hets = include_hets, include_models = include_models,
                          verbose = verbose, parallel = parallel, n_cores = n_cores,
                          het_label = het_label)
  # optionally visualize the cross-validation (honoring custom case labels)
  if (plot_cv_results)
    suppressWarnings(plot_locus_perm_cv(fit_cv, het_label = het_label, neg_label = neg_label))

  # pick the best model by the chosen metric
  model_selected <- .hc_select_model(fit_cv, model_selection_parameter)
  if (verbose)
    message("Best model by ", model_selection_parameter, ": ",
            if (model_selected == "knn") "k-nearest neighbors" else "random forest",
            ". Moving to forward prediction.")

  # a reusable training-argument list for the forward-prediction model(s)
  train_args <- list(geno_mat = train_geno_mat, gene_file = train_gene_file,
                     gene_name = gene_name, marker_info = marker_info,
                     chromosome = chromosome, ncor_markers = ncor_markers,
                     n_neighbors = n_neighbors, include_hets = include_hets,
                     verbose = verbose, set_seed = set_seed,
                     models_request = model_selected, het_label = het_label)

  # ---- Step IIA: single seeded model ----
  if (!predict_by_vote) {
    # train + predict with the best model, falling back to the other on failure
    fp <- .hc_forward_predict(train_args, geno_mat, testing_genotypes, model_selected, verbose)
    # return the single-model result (prediction_model reflects the model used)
    return(list(method = "single model - single prediction",
                cross_validation_results = fit_cv,
                prediction_model = fp$fit,
                predictions = fp$pred))
  }

  # ---- Step IIB: consensus by voting over many random models ----
  # closure that trains+predicts one vote, falling back to the other model if
  # the best model fails (e.g. KNN ties) for that random draw
  one_vote <- function() {
    .hc_forward_predict(train_args, geno_mat, testing_genotypes, model_selected, verbose = FALSE)$pred
  }

  if (parallel) {
    # default to all cores but one when unset
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    # start the cluster and guarantee teardown
    cluster <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cluster)
    on.exit({ parallel::stopCluster(cluster); foreach::registerDoSEQ() }, add = TRUE)
    if (verbose) message("Running votes in parallel; text output suppressed.")
    # export the public functions and internal helpers to the workers
    helper_env <- environment(.hc_check_inputs)
    exports <- c("one_vote", "locus_train", "locus_pred", "train_args", "geno_mat",
                 "testing_genotypes", grep("^\\.hc_", ls(helper_env, all.names = TRUE), value = TRUE))
    # collect one prediction frame per vote
    pred_list <- foreach::foreach(i = 1:n_votes,
                                  .packages = c("caret", "randomForest"),
                                  .export = exports) %dopar% {
      one_vote()
    }
  } else {
    # sequential voting
    pred_list <- lapply(seq_len(n_votes), function(i) one_vote())
  }

  # tabulate the votes into a vote matrix and a majority-rule consensus
  tab <- .hc_tabulate_votes(pred_list)
  # return the voting result
  list(method = "multiple models - majority rule",
       cross_validation_results = fit_cv,
       predictions = tab$votes,
       consensus_predictions = tab$consensus)
}
