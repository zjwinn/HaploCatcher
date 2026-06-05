#' Haplotype Prediction: Training Models for Forward Prediction
#'
#' Trains KNN and/or RF models on the full training data for use in forward
#' prediction of lines that have no locus record. Shares all data preparation
#' and model-fitting logic with [locus_cv()].
#'
#' @param geno_mat An imputed, number-coded genotypic matrix with n rows of individuals and m columns of markers. Row names are genotype IDs; column names are marker IDs. Missing data are not allowed. Numeric coding may vary as long as it is consistent across markers.
#' @param gene_file A data frame with at least the columns 'Gene', 'FullSampleName', and 'Call'. 'Gene' is the gene each observation belongs to, 'FullSampleName' matches a column name in the genotypic matrix, and 'Call' is the marker call for that genotype.
#' @param gene_name A character string matching a value in the 'Gene' column of `gene_file`.
#' @param marker_info A data frame with the columns 'Marker', 'Chromosome', and 'BP_Position'. Every marker in the genotypic matrix must be listed. If positions are unavailable a numeric dummy (1..m) may be used.
#' @param chromosome A character string matching a value in the 'Chromosome' column of `marker_info`.
#' @param ncor_markers Number of top correlated markers to retain for training. Default 50.
#' @param n_neighbors Number of neighbors to consider in KNN. Default 50.
#' @param include_hets Logical; keep heterozygous calls. Default FALSE.
#' @param verbose Logical; print progress and tables. Default FALSE.
#' @param set_seed Numeric seed for reproducibility, or NULL. Default NULL.
#' @param models_request Which models to train: "knn", "rf", or "all". Default "all".
#' @param graph Logical; draw the marker-correlation diagnostic. Default FALSE.
#' @param het_label Optional character vector of `Call` values to treat as heterozygous. When NULL (default), calls containing the prefix "het_" are used.
#'
#' @return A list with `seed`, `models_request`, `trained_models`, and `data`
#'   (the training frame). `trained_models` is a single caret model when one
#'   model was requested, or a list with `knn` and `rf` when "all".
#'
#' @export
#'
#' @examples
#'
#' #set seed for reproducible sampling
#' set.seed(022294)
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
#' #Note: in practice you would have something like a gene file
#' #that does not contain any lines you are trying to predict.
#' #However, this is for illustrative purposes on how to run the function
#'
#' #sample data in the gene_comp file to make a traning population
#' train<-gene_comp[gene_comp$FullSampleName %in%
#'                    sample(gene_comp$FullSampleName,
#'                           round(length(gene_comp$FullSampleName)*0.8),0),]
#'
#' #pull vector of names, not in the train, for forward prediction
#' test<-gene_comp[!gene_comp$FullSampleName
#'                 %in% train$FullSampleName,
#'                 "FullSampleName"]
#'
#' #run the function with hets
#' fit<-locus_train(geno_mat=geno_mat, #the genotypic matrix
#'                  gene_file=train, #the gene compendium file
#'                  gene_name="sst1_solid_stem", #the name of the gene
#'                  marker_info=marker_info, #the marker information file
#'                  chromosome="3B", #name of the chromosome
#'                  ncor_markers=2, #number of markers to retain
#'                  n_neighbors=3, #number of neighbors
#'                  include_hets=FALSE, #include hets in the model
#'                  verbose = FALSE, #allows for text and graph output
#'                  set_seed = 022294, #sets a seed for reproduction of results
#'                  models_request = "knn") #sets what models are requested
#'
#' #predict the lines in the test population
#' pred<-locus_pred(locus_train_results=fit,
#'                  geno_mat=geno_mat,
#'                  genotypes_to_predict=test)
#'
#' #see predictions
#' head(pred)
#'
#' @importFrom randomForest randomForest
#' @importFrom ggplot2 ggplot
#' @importFrom caret train knn3 createFolds trainControl confusionMatrix

locus_train <- function(geno_mat, gene_file, gene_name, marker_info, chromosome,
                        ncor_markers = 50, n_neighbors = 50,
                        include_hets = FALSE, verbose = FALSE,
                        set_seed = NULL, models_request = "all", graph = FALSE,
                        het_label = NULL) {

  # seed the RNG (NULL clears any prior seed) for reproducible model fitting
  set.seed(set_seed)
  # validate data objects and identifiers
  .hc_check_inputs(geno_mat, gene_file, marker_info, gene_name, chromosome)
  # validate the model request up front
  if (!models_request %in% c("knn", "rf", "all"))
    stop("'models_request' must be one of 'knn', 'rf', or 'all'!", call. = FALSE)
  # validate the het flag and any custom het labels
  .hc_assert_logical(include_hets, "include_hets")
  .hc_assert_labels(het_label, "het_label")

  # build the gene-specific, het-filtered classification table
  classification <- .hc_classification(gene_file, gene_name, include_hets, verbose, het_label)
  # call levels used for consistent factor handling
  call_levels <- unique(classification$Call)

  # pull the chromosome marker block for the training individuals
  gc <- .hc_geno_chr(geno_mat, marker_info, chromosome, classification)
  # genotype block ordered by individual name
  train_geno <- gc$geno[order(rownames(gc$geno)), , drop = FALSE]
  # classification ordered to match the genotype block
  ord_cls <- classification[order(classification$FullSampleName), , drop = FALSE]
  # select the top correlated markers across all training data
  corr <- .hc_top_markers(train_geno, ord_cls$Call, gc$markers,
                          ncor_markers, graph, chromosome, gene_name)

  # assemble the model-ready training frame
  training <- .hc_assemble(classification, gc$geno, corr$marker, call_levels)
  # optionally report call frequencies
  if (verbose) .hc_freq_table(training$Call, "Frequency of Calls in Training")

  # fit the requested model set
  fits <- .hc_fit_models(training, n_neighbors, models_request, verbose)
  # for a single-model request, return the bare model (not a one-element list)
  trained_models <- if (models_request == "all") fits else fits[[1]]

  # record the seed actually used (or a sentinel when none was set)
  seed <- if (is.null(set_seed)) "no_seed_set" else set_seed
  # return the training result
  list(seed = seed, models_request = models_request,
       trained_models = trained_models, data = training)
}
