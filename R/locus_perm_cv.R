#' Haplotype Prediction: Permutation Cross Validation of KNN and RF Models
#'
#' Repeats [locus_cv()] over many random partitions (permutations) and
#' summarizes the overall and by-class performance of the KNN and RF models.
#' Can run sequentially or in parallel.
#'
#' @param n_perms Number of permutations to perform. Default 30.
#' @param geno_mat An imputed, number-coded genotypic matrix with n rows of individuals and m columns of markers. Row names are genotype IDs; column names are marker IDs. Missing data are not allowed. Numeric coding may vary as long as it is consistent across markers.
#' @param gene_file A data frame with at least the columns 'Gene', 'FullSampleName', and 'Call'. 'Gene' is the gene each observation belongs to, 'FullSampleName' matches a column name in the genotypic matrix, and 'Call' is the marker call for that genotype.
#' @param gene_name A character string matching a value in the 'Gene' column of `gene_file`.
#' @param marker_info A data frame with the columns 'Marker', 'Chromosome', and 'BP_Position'. Every marker in the genotypic matrix must be listed. If positions are unavailable a numeric dummy (1..m) may be used.
#' @param chromosome A character string matching a value in the 'Chromosome' column of `marker_info`.
#' @param ncor_markers Number of top correlated markers to retain for training. Default 50.
#' @param n_neighbors Number of neighbors to consider in KNN. Default 50.
#' @param percent_testing Proportion of data reserved for validation, strictly between 0 and 1. Default 0.20.
#' @param percent_training Proportion of data used for training, strictly between 0 and 1. Default 0.80.
#' @param include_hets Logical; keep heterozygous calls. Default FALSE.
#' @param include_models Logical; keep the trained models in each permutation (large). Default FALSE.
#' @param verbose Logical; print per-permutation progress. Default FALSE.
#' @param parallel Logical; run permutations in parallel. Default FALSE. When TRUE, textual/graphical feedback is suppressed.
#' @param n_cores Number of cores for parallel processing. If NULL and `parallel = TRUE`, uses all available cores minus one.
#' @param het_label Optional character vector of `Call` values to treat as heterozygous. When NULL (default), calls containing the prefix "het_" are used.
#'
#' @return A list with `Overall_Parameters`, `By_Class_Parameters`,
#'   `Overall_Summary`, `By_Class_Summary`, and `Raw_Permutation_Info`.
#'
#' @export
#'
#' @examples
#'
#' #read in the genotypic data matrix
#' data("geno_mat")
#'
#' #read in the marker information
#' data("marker_info")
#'
#' #read in the gene compendium file
#' data("gene_comp")
#'
#' #run permutational analysis - commented out for package specifications
#' #to run, copy and paste without '#' into the console
#'
#' #fit<-locus_perm_cv(n_perms = 10, #the number of permutations
#' #                   geno_mat=geno_mat, #the genotypic matrix
#' #                   gene_file=gene_comp, #the gene compendium file
#' #                   gene_name="sst1_solid_stem", #the name of the gene
#' #                   marker_info=marker_info, #the marker information file
#' #                   chromosome="3B", #name of the chromosome
#' #                   ncor_markers= 25, #number of markers to retain
#' #                   n_neighbors = 25, #number of nearest-neighbors
#' #                   percent_testing=0.2, #percentage of genotypes in the validation set
#' #                   percent_training=0.8, #percentage of genotypes in the training set
#' #                   include_hets=FALSE, #excludes hets in the model
#' #                   include_models=FALSE, #excludes models in results object
#' #                   verbose = FALSE) #excludes text
#'
#' @importFrom randomForest randomForest
#' @importFrom ggplot2 ggplot
#' @importFrom caret train knn3 createFolds trainControl confusionMatrix
#' @importFrom foreach %dopar% foreach

locus_perm_cv <- function(n_perms = 30, geno_mat, gene_file, gene_name, marker_info,
                          chromosome, ncor_markers = 50, n_neighbors = 50,
                          percent_testing = 0.2, percent_training = 0.8,
                          include_hets = FALSE, include_models = FALSE,
                          verbose = FALSE, parallel = FALSE, n_cores = NULL,
                          het_label = NULL) {

  # validate the logical switches
  .hc_assert_logical(parallel, "parallel")
  .hc_assert_logical(include_hets, "include_hets")
  .hc_assert_logical(include_models, "include_models")
  .hc_assert_logical(verbose, "verbose")
  # validate the permutation count and any custom het labels
  .hc_assert_count(n_perms, "n_perms")
  .hc_assert_labels(het_label, "het_label")
  # n_cores must be NULL or a positive number
  if (!is.null(n_cores)) .hc_assert_count(n_cores, "n_cores")

  # a single argument list reused for every permutation keeps the calls DRY
  cv_args <- list(geno_mat = geno_mat, gene_file = gene_file, gene_name = gene_name,
                  marker_info = marker_info, chromosome = chromosome,
                  ncor_markers = ncor_markers, n_neighbors = n_neighbors,
                  percent_testing = percent_testing, percent_training = percent_training,
                  include_hets = include_hets, include_models = include_models,
                  verbose = FALSE, graph = FALSE, het_label = het_label)

  if (parallel) {

    # default to all cores but one when unset
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    # spin up a cluster and register it as the parallel backend
    cluster <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cluster)
    # always tear the cluster down, even on error
    on.exit({ parallel::stopCluster(cluster); foreach::registerDoSEQ() }, add = TRUE)
    # tell the user parallel mode suppresses feedback
    if (verbose) message("Running permutations in parallel; text output suppressed.")

    # export locus_cv plus all internal helpers to the workers (works whether
    # the package is installed or the single source file was sourced)
    helper_env <- environment(.hc_check_inputs)
    exports <- c("locus_cv", "cv_args", grep("^\\.hc_", ls(helper_env, all.names = TRUE), value = TRUE))
    # run the permutations across the cluster
    results <- foreach::foreach(i = 1:n_perms,
                                .packages = c("caret", "randomForest"),
                                .export = exports) %dopar% {
      do.call(locus_cv, cv_args)
    }

  } else {

    # announce sequential progress
    if (verbose) message("Conducting permutational cross validation for ", gene_name, "...")
    # pre-allocate the results list
    results <- vector("list", n_perms)
    # run each permutation in turn
    for (i in seq_len(n_perms)) {
      results[[i]] <- do.call(locus_cv, cv_args)
      # report percent complete when verbose
      if (verbose) message(round(i / n_perms * 100), "% complete")
    }

  }

  # collapse the raw permutation results into the five-element summary object
  .hc_summarize_perm(results, gene_name)
}
