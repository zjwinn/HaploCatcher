#' Haplotype Prediction: Cross Validation of KNN and RF Models
#'
#' Performs one round of the cross-validation featured in Winn et al. (2022):
#' a random partition of the training data trains KNN and RF models, and a
#' reserved test partition validates them. This is a single permutation; use
#' [locus_perm_cv()] to repeat it.
#'
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
#' @param include_models Logical; keep the trained models in the result (large). Default FALSE.
#' @param verbose Logical; print progress and tables. Default TRUE.
#' @param graph Logical; draw the marker-correlation diagnostic. Default FALSE.
#' @param het_label Optional character vector of `Call` values to treat as heterozygous. When NULL (default), calls containing the prefix "het_" are used.
#'
#' @return A list with `data_frames` (training and test frames), `test_predictions`
#'   (per-model prediction frames), `confusion_matrices` (per-model confusion
#'   objects), and, when `include_models = TRUE`, `trained_models`.
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
#' #run the function without hets for a very limited number of markers and neighbors
#' #due to requirements by cran, this must be commented out
#' #to run, place this code in the console and remove comments
#' #fit<-locus_cv(geno_mat=geno_mat, #the genotypic matrix
#' #             gene_file=gene_comp, #the gene compendium file
#' #             gene_name="sst1_solid_stem", #the name of the gene
#' #             marker_info=marker_info, #the marker information file
#' #             chromosome="3B", #name of the chromosome
#' #             ncor_markers=2, #number of markers to retain
#' #             n_neighbors=1, #number of neighbors
#' #             percent_testing=0.2, #percentage of genotypes in the validation set
#' #             percent_training=0.8, #percentage of genotypes in the training set
#' #             include_hets=FALSE, #include hets in the model
#' #             include_models=TRUE, #include models in the final results
#' #             verbose=TRUE, #allows text output
#' #             graph=TRUE) #allows graph output
#'
#' @importFrom randomForest randomForest
#' @importFrom ggplot2 ggplot
#' @importFrom caret train knn3 createFolds trainControl confusionMatrix

locus_cv <- function(geno_mat, gene_file, gene_name, marker_info, chromosome,
                     ncor_markers = 50, n_neighbors = 50,
                     percent_testing = 0.2, percent_training = 0.8,
                     include_hets = FALSE, include_models = FALSE,
                     verbose = TRUE, graph = FALSE, het_label = NULL) {

  # validate data objects and identifiers in one place
  .hc_check_inputs(geno_mat, gene_file, marker_info, gene_name, chromosome)
  # validate logical and proportion arguments
  .hc_assert_logical(include_hets, "include_hets")
  .hc_assert_logical(include_models, "include_models")
  .hc_assert_prop(percent_training, "percent_training")
  .hc_assert_prop(percent_testing, "percent_testing")
  .hc_assert_labels(het_label, "het_label")
  # the two partitions cannot claim more than all of the data
  if (percent_training + percent_testing > 1)
    stop("'percent_training' + 'percent_testing' cannot exceed 1!", call. = FALSE)

  # build the gene-specific, het-filtered classification table
  classification <- .hc_classification(gene_file, gene_name, include_hets, verbose, het_label)
  # remember the call levels so factors stay consistent across partitions
  call_levels <- unique(classification$Call)

  # randomly choose the training individuals
  train_ids <- sample(classification$FullSampleName,
                      size = floor(percent_training * nrow(classification)))
  # split the classification into training and testing partitions
  training_cls <- classification[classification$FullSampleName %in% train_ids, , drop = FALSE]
  testing_cls  <- classification[!classification$FullSampleName %in% train_ids, , drop = FALSE]

  # pull the chromosome marker block for all involved individuals
  gc <- .hc_geno_chr(geno_mat, marker_info, chromosome, classification)

  # genotype block restricted to the training individuals, ordered by name
  train_geno <- gc$geno[rownames(gc$geno) %in% training_cls$FullSampleName, , drop = FALSE]
  train_geno <- train_geno[order(rownames(train_geno)), , drop = FALSE]
  # training classification ordered to match the genotype block
  ord_train_cls <- training_cls[order(training_cls$FullSampleName), , drop = FALSE]
  # select the top correlated markers using only the training data
  corr <- .hc_top_markers(train_geno, ord_train_cls$Call, gc$markers,
                          ncor_markers, graph, chromosome, gene_name)

  # assemble model-ready training and test frames on the selected markers
  training <- .hc_assemble(training_cls, gc$geno, corr$marker, call_levels)
  testing  <- .hc_assemble(testing_cls,  gc$geno, corr$marker, call_levels)

  # optionally report call frequencies in each partition
  if (verbose) {
    .hc_freq_table(training$Call, "Frequency of Calls in Training")
    .hc_freq_table(testing$Call,  "Frequency of Calls in Test")
  }

  # fit both models with the shared control/grids
  fits <- .hc_fit_models(training, n_neighbors, "all", verbose)
  # predict the test partition and build confusion matrices
  ev <- .hc_eval_models(fits, testing, gene_name, call_levels, verbose)

  # assemble the result, optionally carrying the (large) trained models
  out <- list(data_frames        = list(training = training, test = testing),
              test_predictions   = ev$preds,
              confusion_matrices = ev$confu)
  if (include_models) out <- c(out["data_frames"],
                               list(trained_models = fits),
                               out[c("test_predictions", "confusion_matrices")])
  # return the cross-validation result
  out
}
