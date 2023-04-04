#' Haplotype Prediction: Using Trained Models to Make Predictions
#'
#' @param locus_train_results This object is a the results object from the "locus_train" function.
#' @param geno_mat This is a genotypic matrix with genotypes of individuals you have genotyped and characterized for the locus/QTL/gene of interest and individuals that have only been genotyped with a genome wide marker panel. It is important to note that the markers in the genome wide panel *must* be shared between training population and the population you wish to forward predict. In the case of genotyping-by-sequencing markers, the training and test populations should be discovered and produced together into a genotype file. All marker platforms, however, are compatable, as long as the training and forward prediction population share the same markers genome-wide.
#' @param genotypes_to_predict This is a character vector of genotypes in the geno_mat which the user wishes to predict. If this object contains names in the training population, they will be omitted in the results to avoid bias.
#'
#' @return a data frame with two columns: genotype names and predictions.
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
#'                  models = "knn") #sets what models are requested
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
#' @importFrom lattice qq
#' @importFrom ggplot2 ggplot
#' @importFrom caret train
#' @importFrom caret knn3
#' @importFrom foreach %dopar%

locus_pred<-function(locus_train_results,
                     geno_mat,
                     genotypes_to_predict){

  a<-locus_train_results
  b<-geno_mat
  c<-genotypes_to_predict

  if(a$models_request=="all"){

    #pull models
    knn<-a$trained_models$knn
    rf<-a$trained_models$rf

    #subset to only have the indiviudals who you want to predict
    b<-b[base::rownames(b) %in% c, base::colnames(b) %in% base::colnames(a$data)[3:base::ncol(a$data)]]

    #predict
    pred_1=stats::predict(knn, b)
    pred_2=stats::predict(rf, b)

    #make a prediction object
    pred<-base::data.frame(FullSampleName=rownames(b),
                           Prediction_KNN=pred_1,
                           Prediction_RF=pred_2)

    #make results
    results<-pred

    #return
    return(results)

  }else if(a$models_request=="knn"){

    #pull models
    knn<-a$trained_models

    #subset to only have the indiviudals who you want to predict
    b<-b[base::rownames(b) %in% c, base::colnames(b) %in% base::colnames(a$data)[3:base::ncol(a$data)]]

    #predict
    pred_1=stats::predict(knn, b)

    #make a prediction object
    pred<-base::data.frame(FullSampleName=rownames(b),
                           Prediction_KNN=pred_1)

    #make results
    results<-pred

    #return
    return(results)

  }else if(a$models_request=="rf"){

    #pull models
    rf<-a$trained_models

    #subset to only have the indiviudals who you want to predict
    b<-b[base::rownames(b) %in% c, base::colnames(b) %in% base::colnames(a$data)[3:base::ncol(a$data)]]

    #predict
    pred_2=stats::predict(rf, b)

    #make a prediction object
    pred<-base::data.frame(FullSampleName=rownames(b),
                           Prediction_RF=pred_2)

    #make results
    results<-pred

    #return
    return(results)

  }else{

    base::stop("Error: Is the 'geno_train_results' correct? Check your inputs!")

  }

}
