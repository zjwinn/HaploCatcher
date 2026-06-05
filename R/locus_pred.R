#' Haplotype Prediction: Using Trained Models to Make Predictions
#'
#' Applies the models from [locus_train()] to forward-predict the haplotype of
#' genotypes that have genome-wide marker data but no locus record.
#'
#' @param locus_train_results The list returned by [locus_train()].
#' @param geno_mat A genotypic matrix containing the genotypes to predict. The genome-wide markers must be shared with the training population.
#' @param genotypes_to_predict A character vector of genotype names (rows of `geno_mat`) to predict. Names that were in the training data are dropped to avoid bias.
#'
#' @return A data frame with `FullSampleName` and one prediction column per
#'   trained model (`Prediction_KNN` and/or `Prediction_RF`).
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
#' @importFrom caret train knn3

locus_pred <- function(locus_train_results, geno_mat, genotypes_to_predict) {

  # the model request tells us which models the training object holds
  req <- locus_train_results$models_request
  # guard against a malformed training object
  if (is.null(req) || !req %in% c("knn", "rf", "all"))
    stop("'locus_train_results' is malformed; was it produced by locus_train()?", call. = FALSE)

  # the markers used in training are columns 3..n of the stored training frame
  marker_cols <- colnames(locus_train_results$data)[-(1:2)]
  # subset the genotype matrix to the requested individuals and training markers
  x <- geno_mat[rownames(geno_mat) %in% genotypes_to_predict,
                colnames(geno_mat) %in% marker_cols, drop = FALSE]
  # nothing to predict if no requested individuals are present
  if (nrow(x) == 0L)
    stop("None of 'genotypes_to_predict' were found in 'geno_mat'!", call. = FALSE)

  # start the prediction frame with the individual names
  out <- data.frame(FullSampleName = rownames(x), stringsAsFactors = FALSE)
  # pull KNN predictions when a KNN model is available
  if (req %in% c("all", "knn")) {
    knn <- if (req == "all") locus_train_results$trained_models$knn else locus_train_results$trained_models
    out$Prediction_KNN <- stats::predict(knn, x)
  }
  # pull RF predictions when an RF model is available
  if (req %in% c("all", "rf")) {
    rf <- if (req == "all") locus_train_results$trained_models$rf else locus_train_results$trained_models
    out$Prediction_RF <- stats::predict(rf, x)
  }
  # return the prediction frame
  out
}
